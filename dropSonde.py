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

class Sonde:
	def __init__(self, filePath):
		self.filePath = filePath

	# Get both core and drop sonde data as two data frames
	def getData(self):
		self.dropTimes = []
		self.sondes = []
		global variablesSonde
		path = [os.path.join(dirpath, f)
			for dirpath, dirnames, files in os.walk("/".join(self.filePath.split("/")[:-1]))
			for f in fnmatch.filter(files, "faam-dropsonde*proc.nc")]
		
		for n, i in enumerate(path):
			dfSonde = pd.DataFrame(columns = variablesSonde)
			fh = Dataset(i, "r")
			rawTime = fh.variables["base_time"].string.split(" ")[3]
			self.dropTimes.append(int(rawTime[6:8]) + 60*( int(rawTime[3:5]) + 60*(int(rawTime[0:2]))))
			for i in variablesSonde:
				var = fh.variables[i][:]
				var = var.ravel()
				dfSonde[i] = var
			self.sondes.append(dfSonde)
		

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

		for i in self.dropTimes:
			dfdropTimes = dfdropTimes.append(dfCore[(dfCore["TIME"] >= i-100) & (dfCore["TIME"] <= i+100)], sort = True)
		self.coreData = dfdropTimes

		# print(self.filePath)
		# print("Drop times /UTC after midnight: {}".format(self.dropTimes))
		# print(self.coreData)
		# print(self.sondes)

	def standardizeTime(self):
		for i, sond in enumerate(self.sondes):
			sond["TIME"] = sond["time_offset"].add(self.dropTimes[i]) 

	def mergeData(self):
		self.sondeData = pd.concat(self.sondes)
		self.coreData["TIME"] = self.coreData["TIME"].astype(float) 
		self.coreData = pd.merge(self.coreData, self.sondeData, on="TIME", how="outer")
		self.coreData = self.coreData[self.coreData["pres"] != -999.0]

	def plotData(self):
		plt.plot(self.coreData["TIME"], self.coreData["PS_RVSM"])
		plt.scatter(self.coreData["TIME"], self.coreData["pres"], s=0.5)
		for i in self.dropTimes:
			plt.axvline(i, c = "red")
		plt.show()

# Get flights where dropsondes were depolyed from 2015-2019
def sondeFilePaths():
	fnames = []
	for i in range(5,10):
		path = "/media/faamarchive/badcMirror/data/201{}/".format(i)
		fnames += [os.path.join(dirpath, f)
			for dirpath, dirnames, files in os.walk(path)
			for f in fnmatch.filter(files, "faam-dropsonde*[!raw].nc")]
	return(fnames)

fnames = sondeFilePaths()


sonde = Sonde("/media/faamarchive/badcMirror/data/2019/c153-mar-13/core_processed/faam-dropsonde_faam_20190313115953_r0_c153_proc.nc")
sonde.getData()
sonde.standardizeTime()
sonde.mergeData()
sonde.plotData()

#print(sonde.filePath)
#print(sonde.dropTimes)
#print(sonde.coreData)
# drop.getcoreData()
# print(drop.coreFiles)
# print(drop.dropSondeFiles)

