import numpy as np
from netCDF4 import Dataset
import scipy
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
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
variablesSonde =["time","time_offset","pres","alt","tdry"]


# def func(x,a,b,c,d):
# 		return a * ((b * (1/c+x)) **d)

# def func(x,a):
# 		return -a*x

def func(x,a,b,c):
	return a * np.exp(-b * x) + c

class Flight:
	def __init__(self, filePath):
		self.filePath = filePath
		self.functions = []
		self.errorData = pd.DataFrame(columns=["subError","fullError","altError","altitude","pressure","year","flight"])
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

	# Extrapolate sonde data to relase point using first 200 seconds of data 
	def subExtrapolate(self):
		self.functionsSub = []
		self.subErrors = np.array([])
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[(self.coreData["FLAG"] == i) & (self.coreData["TIME"] <= dropTime + 200)].sort_values("TIME")
			# Save a list of derrived functions
			self.functionsSub.append(np.polyfit(data["TIME"],data["pres"], 2))
			
			# Find expected and observed
			f = np.poly1d(np.polyfit(data["TIME"],data["pres"], 2))
			exp = np.average(f(self.coreData[self.coreData["TIME"] == dropTime]["TIME"]))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			
			# Update data frame
			self.subErrors = np.append(self.subErrors, exp-obs)
														
	# Extrapolate sonde data using to release point using entire dataset		
	def fullExtrapolate(self):
		self.functionsFull = []
		self.fullErrors = np.array([])
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[(self.coreData["FLAG"] == i)].sort_values("TIME")
			# Save a list of derrived functions
			self.functionsFull.append(np.polyfit(data["TIME"],data["pres"], 2))
			
			# Find expected and observed
			f = np.poly1d(np.polyfit(data["TIME"],data["pres"], 2))
			exp = np.average(f(self.coreData[self.coreData["TIME"] == dropTime]["TIME"]))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])

			# Update data frame
			self.fullErrors = np.append(self.fullErrors,exp-obs) 

	# Extrapolate sonde data using pressure altitude curve 
	def altitudeExtrapolate(self):
		self.functionsAlt = []
		self.altErrors = np.array([])
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[self.coreData["FLAG"] == i].sort_values(by = "alt")
			
			# Save a list of derrived functions
			popt, pcov = curve_fit(func, data["alt"].values, data["pres"].values)
			self.functionsAlt.append(popt)
			
			# Find expected and observed
			exp = np.average(func(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"], *popt))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
						
			# If error is anomolous return NaN
			if (exp-obs) < 10 and (exp-obs) >-10:
				self.altErrors = np.append(self.altErrors,exp-obs)
			else:
				self.altErrors = np.append(self.altErrors,np.nan)

	# Calculate indicated mach number
	def calc_mach(self):
		self.coreData['MACH'] = self.coreData['IAS_RVSM'] / (340.294 * np.sqrt(self.coreData['PS_RVSM'] / 1013.25))

	# Collect and generate the data set
	def generateDataSet(self):
		for i, dropTime in enumerate(self.dropTimes):
			#print(np.average(self.coreData[self.coreData["TIME"] == dropTime]["MACH"]))
			self.errorData = self.errorData.append(
				{"subError":self.subErrors[i],
				"fullError":self.fullErrors[i],
				"altError" :self.altErrors[i],
				"altitude" :np.average(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"]),
				"pressure" :np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"]),
				"mach"     :np.average(self.coreData[self.coreData["TIME"] == dropTime]["MACH"]),
				"dpressure":np.average(self.coreData[self.coreData["TIME"] == dropTime]["Q_RVSM"]),
				"flight"   :self.filePath.split("/")[6][0:4],
				"year"     :self.filePath.split("/")[5]}, ignore_index = True)
			self.errorData = self.errorData.reset_index(drop = True)
	
	# Graph pressure error as a fucntion of time
	def plotDataTime(self):
		# Scatter flight pressure and dropsonde pressure as a function of time
		plt.scatter(self.coreData["TIME"], self.coreData["PS_RVSM"], s=0.5)
		plt.scatter(self.coreData["TIME"], self.coreData["pres"], s=0.5, c = self.coreData["FLAG"])
		# For all dropTimes get fit and extrapolate a line plot from the point of release to first point of data recording
		for i, dropTime in enumerate(self.dropTimes):
			# Initial drop
			plt.axvline(dropTime, c = "red", alpha=0.5)
			# Get data
			data = self.coreData[self.coreData["FLAG"] == i].sort_values(by = "TIME")
			f = np.poly1d(self.functionsFull[i])
			# Plot Data
			extrapolationData = self.coreData[(self.coreData["TIME"] >= dropTime - 50) & (self.coreData["TIME"] <= dropTime + 100)].sort_values("TIME")
			extrapolationData["F(x)"] = f(extrapolationData["TIME"])
			plt.plot(extrapolationData["TIME"], extrapolationData["F(x)"], c="black")
		plt.show()

	# Graph pressure error as a fucntion of altitude
	def plotDataAlt(self):
		plt.scatter(self.coreData["alt"], self.coreData["pres"], s = 1.5, alpha = 0.5)
		for i, dropTime in enumerate(self.dropTimes):
			# Initital pressure at drop
			plt.axvline(np.average(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"]), c = "red")
			plt.axhline(np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"]), c = "red", alpha = 0.5)

			# Get data
			data = self.coreData[self.coreData["FLAG"] == i].sort_values(by = "alt")
			popt = self.functionsAlt[i]
				
			# Plot data
			plt.plot(data["alt"], func(data["alt"], *popt))
		plt.show()

