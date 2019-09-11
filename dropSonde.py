import numpy as np
from netCDF4 import Dataset
import scipy
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import glob
import fnmatch
import math
import pandas as pd
import pdb
import datetime
import warnings

warnings.filterwarnings("ignore")

# Globals
variablesCore = ["Time", "PS_RVSM", "ALT_GIN", "IAS_RVSM", "TAT_DI_R", "TAS_RVSM","Q_RVSM"]
variablesSonde = ["time","time_offset","pres","alt","tdry"]
index0 = -1

def func(h):
	global P0
	global H0
	global T0
	# h0 = 0.0
	# T0 = 288.15
	# P0 = 1013.2500
	L0 = -0.0065
	go = 9.80665
	M = 0.0289644
	R = 8.31432	
	# Calulate pressure from standard atmosphere
	return(P0*(T0/ (T0 + L0 * (h - H0)))**(go * M / (R * L0)))

class Flight:
	def __init__(self, filePath):
		self.filePath = filePath
		self.functions = []
		self.errorData = pd.DataFrame(columns=["error","altitude","pressure","year","flight"])
		print(self.filePath)

	# Get both core and drop sonde data as individual data frames
	def getData(self):
		self.dropTimes = []
		self.sondes = []

		global variablesSonde

		path = [os.path.join(dirpath, f)
			for dirpath, dirnames, files in os.walk(self.filePath)
			for f in fnmatch.filter(files, "faam-dropsonde*proc.nc")]
		
		rev = 0
		appended = 1
		timeStamps = []
		pathFinal = []
		while appended > 0:
			appended = 0
			for i in path:
				#Check if i has revision number and if its not in timeStamps
				if("r{}".format(rev) in i) and i.split("_")[3] not in timeStamps: 
					pathFinal.append(i)
					timeStamps.append(i.split("_")[3])
					appended += 1
				#If it does have revision number but is in timeStamps
				elif("r{}".format(rev) in i):
					#Loop over time stamps and replace file in pathFianl
					for j in range(len(timeStamps)):
						if i.split("_")[3] == timeStamps[j]:
							pathFinal[j] = i
							appended += 1
			rev+= 1
		path = pathFinal
		# Retrive sonde data and drop times
		for i, path in enumerate(path):
			dfSonde = pd.DataFrame(columns = variablesSonde)
			fh = Dataset(path, "r")
			rawTime = fh.variables["base_time"].string.split(" ")[3]
			self.dropTimes.append(int(rawTime[6:8]) + 60*( int(rawTime[3:5]) + 60*(int(rawTime[0:2]))))
			for j in variablesSonde:
				var = fh.variables[j][:]
				var = var.ravel()
				dfSonde[j] = var
			dfSonde["FLAG"] = i
			self.sondes.append(dfSonde)

		
		# Retrive core data
		global variablesCore
		path = [os.path.join(dirpath, f)
				for dirpath, dirnames, files in os.walk(self.filePath)
				for f in fnmatch.filter(files, "core_faam*[!1hz].nc")]
		for i in path:
			fh = Dataset(i, "r")
		dfCore = pd.DataFrame(columns = variablesCore)
		dfdropTimes = pd.DataFrame(columns = variablesCore)
		
		for i in variablesCore:
			var = fh.variables[i][:]
			var = var.ravel()
			if i != "Time":
				var = var[0::32]
			dfCore[i] = var
			
		dfCore = dfCore.rename(columns = {"Time":"TIME"})
		dfdropTimes = dfdropTimes.rename(columns = {"Time":"TIME"})
		
		# Only keep core data around drop time of sonde 
		for i in self.dropTimes:
			dfdropTimes = dfdropTimes.append(dfCore[(dfCore["TIME"] >= i-100) & (dfCore["TIME"] <= i+100)], sort = True)
		self.coreData = dfdropTimes
		
	
	# Stamdardize the times by adding the offset (time since drop) to the drop time
	def standardizeTime(self):
		for i, sond in enumerate(self.sondes):
			sond["TIME"] = sond["time_offset"].add(self.dropTimes[i]) 
	
	# Merge data into one data frame 
	def mergeData(self):
		self.sondeData = pd.concat(self.sondes)
		self.coreData["TIME"] = self.coreData["TIME"].astype(float) 
		self.coreData = pd.merge(self.coreData, self.sondeData, on="TIME", how="outer")
		self.coreData.loc[self.coreData["pres"] == -999.0, "FLAG"] = -1
		self.coreData = self.coreData.replace(-999.0,np.nan)
		self.coreData = self.coreData.drop_duplicates()
		self.coreData.to_csv("coreData")

	# Extrapolate the ISA from the first recorded sonde point
	def extrapolate(self):
		self.functionsAlt = []
		self.errors = np.array([])
		global index0
		global P0
		global H0
		global T0
		for i, dropTime in enumerate(self.dropTimes):
			# Get data set
			data = self.coreData[self.coreData["FLAG"] == i].sort_values(by = "alt").reset_index(drop=True)
			# Remove weird data point (maybe find true fix later)
			data = data[0:len(data)-1]
			data = data.iloc[index0]
			# Get inital values for ISA curve
			P0 = data["pres"] #hPa
			H0 = data["alt"] #m
			T0 = data["tdry"] + 273.15 #K
			#Find expected and observed and save the difference to an array
			exp = func(np.average(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"]))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			self.errors = np.append(self.errors, obs-exp)

	# Calculate indicated mach number
	def calc_mach(self):
		self.coreData['MACH'] = self.coreData['IAS_RVSM'] / (340.294 * np.sqrt(self.coreData['PS_RVSM'] / 1013.25))

	# Collect and generate the data set
	def generateDataSet(self):
		for i, dropTime in enumerate(self.dropTimes):
			#print(np.average(self.coreData[self.coreData["TIME"] == dropTime]["MACH"]))
			self.errorData = self.errorData.append(
				{"error" :self.errors[i],
				"altitude" :np.average(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"]),
				"pressure" :np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"]),
				"mach"     :np.average(self.coreData[self.coreData["TIME"] == dropTime]["MACH"]),
				"dpressure":np.average(self.coreData[self.coreData["TIME"] == dropTime]["Q_RVSM"]),
				"flight"   :self.filePath.split("/")[6][0:4],
				"year"     :self.filePath.split("/")[5]}, ignore_index = True)
			self.errorData = self.errorData.reset_index(drop = True)
			
	# Graph pressure error as a fucntion of altitude
	def plotData(self):
		global index0
		global P0
		global H0
		global T0		
		for i, dropTime in enumerate(self.dropTimes):
			# Get data
			data = self.coreData[self.coreData["FLAG"] == i].sort_values(by = "alt")
			data = self.coreData[self.coreData["FLAG"] == i].sort_values(by = "alt").reset_index(drop=True)
			# Remove weird data point (maybe find true fix later)
			data = data[0:len(data)-1]
			data = data.iloc[index0]
			# Get inital values for ISA curve
			P0 = data["pres"] #hPa
			H0 = data["alt"] #m
			T0 = data["tdry"] + 273.15 #K
			# Plot data
			plt.cla()
			# Initital pressure at drop
			plt.axvline(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"].values, c = "black", alpha = 0.4)
			plt.axhline(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"].values, c = "black", alpha = 0.5)
			plt.scatter(data["alt"], data["pres"], c="green", label = "Drop sonde", zorder = 10)
			plt.scatter(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"].values,
						self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"].values, c="red", label = "PS_RVSM", zorder = 11)
			plt.plot(np.arange(math.floor(data["alt"]), math.ceil(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"])),
				func(np.arange(math.floor(data["alt"]), math.ceil(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"]))),
				ls = "--",c="blue", label="Extrapolation")
			plt.legend(loc = "upper right")
			plt.ylabel("Static Pressure [hPa]")
			plt.xlabel("Altitude above geodesic [m]")
			plt.title("Pressure Altitude extrapolation from a Drop-Sonde ejected from the ARA")
			plt.show()

