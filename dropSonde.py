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
	# Get both core and drop sonde data as individual data frames
	def getData(self):
		self.dropTimes = []
		self.sondes = []
		global variablesSonde
		path = [os.path.join(dirpath, f)
			for dirpath, dirnames, files in os.walk("/".join(self.filePath.split("/")[:-1]))
			for f in fnmatch.filter(files, "faam-dropsonde*proc.nc")]
		
		# Retrive sonde data and drop times
		for i, path in enumerate(path):
			print(path)
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
				for dirpath, dirnames, files in os.walk("/".join(self.filePath.split("/")[:-1]))
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
		#print("Drop times /UTC after midnight: {}".format(self.dropTimes))
		# print(self.coreData)
		# print(self.sondes)
	
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
		#self.coreData.to_csv("Data/coreData")

	# Extrapolate sonde data to relase point
	def extrapolate(self):
		print("WHAT IS HAPPENIGN?")
		self.functions = []
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[(self.coreData["FLAG"] == i) & (self.coreData["TIME"] <= dropTime + 200)].sort_values("TIME")
			data.to_csv("Data/DataSet{}".format(i))
			z = np.polyfit(data["TIME"],data["pres"], 2)
			self.functions.append(z)
			f = np.poly1d(z)
			exp = f(self.coreData[self.coreData["TIME"] == dropTime]["TIME"])
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			error = error.append(exp-obs)
			alt = np.average(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"])

	# Graph the data (temporary)
	def plotData(self):
		plt.scatter(self.coreData["TIME"], self.coreData["PS_RVSM"], s=0.5)
		plt.scatter(self.coreData["TIME"], self.coreData["pres"], s=0.5, c = self.coreData["FLAG"])
		
		for i, dropTime in enumerate(self.dropTimes):
			plt.axvline(dropTime, c = "red")
			f = np.poly1d(self.functions[i])
			x = self.coreData[self.coreData["FLAG"] == i].sort_values("TIME")
			xExtrapol = self.coreData[(self.coreData["TIME"] >= dropTime - 50) & (self.coreData["TIME"] <= dropTime + 100)].sort_values("TIME")
			plt.plot(xExtrapol["TIME"], f(xExtrapol["TIME"]), c="black")
			plt.plot(x["TIME"], f(x["TIME"]), c = "black")
		plt.show()



