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
		self.errorData = pd.DataFrame(columns=["subError","fullError","altitude","pressure"])
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
		self.functionsSubExrapolation = []
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[(self.coreData["FLAG"] == i) & (self.coreData["TIME"] <= dropTime + 200)].sort_values("TIME")
			# Save list of functions derrived
			self.functionsSubExrapolation.append(np.polyfit(data["TIME"],data["pres"], 2))
			
			# Find expected and observed
			f = np.poly1d(z)
			exp = np.average(f(self.coreData[self.coreData["TIME"] == dropTime]["TIME"]))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			
			# update data frame
			self.errorData["subError"] = exp-obs
			self.errorData["altitude"] = np.average(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"])
			self.errorData["pressure"] = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			self.errorData["flight"] = self.filePath.split("/")[6][0:4]
			self.errorData["year"] = :self.filePath.split("/")[5]
			self.errorData = self.errorData.reset_index(drop = True)
	
	# Extrapolate sonde data using to release point using entire dataset		
	def fullExtrapolate(self):
		self.functions = []
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[(self.coreData["FLAG"] == i)].sort_values("TIME")
			self.functions.append(np.polyfit(data["TIME"],data["pres"], 2))
			f = np.poly1d(np.polyfit(data["TIME"],data["pres"], 2))
			exp = np.average(f(self.coreData[self.coreData["TIME"] == dropTime]["TIME"]))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			self.errorData["fullError"] = exp-obs

	def altitudeExtrapolate(self):
		self.functions = []
		for i, dropTime in enumerate(self.dropTimes):
			data = self.coreData[self.coreData["FLAG"] == i].sort_values(by = "alt")
			popt, pcov = curve_fit(func, data["alt"].values, data["pres"].values)
			self.functions.append(popt)
			exp = np.average(func(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"], *popt))
			obs = np.average(self.coreData[self.coreData["TIME"] == dropTime]["PS_RVSM"])
			error = exp - obs
			if error < 10 and error >-10:
				self.errorData["altError"] = error
			else:
				self.errorData["altError"] = np.nan


	# Graph the data (temporary)
	def plotDataTime(self):
		plt.scatter(self.coreData["TIME"], self.coreData["PS_RVSM"], s=0.5)
		plt.scatter(self.coreData["TIME"], self.coreData["pres"], s=0.5, c = self.coreData["FLAG"])
		for i, dropTime in enumerate(self.dropTimes):
			plt.axvline(dropTime, c = "red")
			f = np.poly1d(self.functions[i])
			x = self.coreData[self.coreData["FLAG"] == i]
			x = x.sort_values(by = "TIME")
			
			xExtrapol = self.coreData[(self.coreData["TIME"] >= dropTime - 50) & (self.coreData["TIME"] <= dropTime + 100)].sort_values("TIME")
			xExtrapol["f"] = f(xExtrapol["TIME"])
			xExtrapol =  xExtrapol.drop_duplicates("TIME")

			plt.plot(xExtrapol["TIME"], f(xExtrapol["TIME"]), c="black")
			#plt.plot(x["TIME"], f(x["TIME"]), c = "black")
		plt.show()

	def plotDataAlt(self):
		for i, dropTime in enumerate(self.dropTimes):
			popt = self.functions[i]			
			x = self.coreData[self.coreData["FLAG"] == i]
			x = x.sort_values(by = "alt")

			# print(x["alt"].reset_index(drop=True))
			# xExtrapol = np.arange(x["alt"][-1],self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"],1)
			

			plt.scatter(x["alt"], x["pres"],s = 1.5, alpha = 0.5)
			plt.plot(x["alt"], func(x["alt"], *popt))

			#plt.plot(xExtrapol, func(xExtrapol), c="Green")

			plt.axvline(np.average(self.coreData[self.coreData["TIME"] == dropTime]["ALT_GIN"]))
		plt.show()

