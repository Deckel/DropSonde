import numpy as np
from netCDF4 import Dataset
import scipy
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import os
import glob
import fnmatch
import math
import pandas as pd
import pdb
import datetime


# Globals
variablesCore = ["Time", "PS_RVSM", "ALT_GIN", "IAS_RVSM", "TAT_DI_R", "TAS_RVSM","Q_RVSM"]
variablesSonde =["time","time_offset","pres","alt","tdry"]

class Flight:
	def __init__(self, filePath):
		self.filePath = filePath
		self.functions = []
		self.errorData = pd.DataFrame(columns=["subError","fullError","Altitude","Pressure"])
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
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[(self.coreData["FLAG"] == i) & (self.coreData["TIME"] <= dropTime + 200)].sort_values("TIME")
			z = np.polyfit(data["TIME"],data["pres"], 2)
			self.functions.append(z)
			f = np.poly1d(z)
			exp = np.average(f(self.coreData[self.coreData["TIME"] == dropTime]["TIME"]))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			self.errorData = self.errorData.append({"subError":exp-obs,
													"Altitude":np.average(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"]),
													"Pressure":np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"]),
													"flight":self.filePath.split("/")[6][0:4],
													"year":self.filePath.split("/")[5]
													}, ignore_index = True)
	
	# Extrapolate sonde data using to release point using entire dataset		
	def fullExtrapolate(self):
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[(self.coreData["FLAG"] == i)].sort_values("TIME")
			z = np.polyfit(data["TIME"],data["pres"], 2)
			self.functions.append(z)
			f = np.poly1d(z)
			exp = np.average(f(self.coreData[self.coreData["TIME"] == dropTime]["TIME"]))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			self.errorData["fullError"] = exp-obs


	# Graph the data (temporary)
	def plotData(self):
		plt.scatter(self.coreData["TIME"], self.coreData["PS_RVSM"], s=0.5)
		plt.scatter(self.coreData["TIME"], self.coreData["pres"], s=0.5, c = self.coreData["FLAG"])
		#print(self.dropTimes)
		for i, dropTime in enumerate(self.dropTimes):
			plt.axvline(dropTime, c = "red")
			f = np.poly1d(self.functions[i])
			x = self.coreData[self.coreData["FLAG"] == i]
			x = x.sort_values(by = "TIME")
			
			xExtrapol = self.coreData[(self.coreData["TIME"] >= dropTime - 50) & (self.coreData["TIME"] <= dropTime + 100)].sort_values("TIME")
			xExtrapol["f"] = f(xExtrapol["TIME"])
			xExtrapol =  xExtrapol.drop_duplicates("TIME")
			#print(i, dropTime)

			plt.plot(xExtrapol["TIME"], f(xExtrapol["TIME"]), c="black")
			#plt.plot(x["TIME"], f(x["TIME"]), c = "black")
		#plt.show()



