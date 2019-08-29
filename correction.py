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
import time
from dropSonde import Flight

# Get flights where dropsondes were depolyed from 2015-2019
def sondeFilePaths():
	fnames = []
	for i in range(4,10):
		path = "/media/faamarchive/badcMirror/data/201{}/".format(i)
		fnames += [os.path.join(dirpath, f)
			for dirpath, dirnames, files in os.walk(path)
			for f in fnmatch.filter(files, "faam-dropsonde*[!raw].nc")]
	for i in range(len(fnames)):
		fnames[i] = "/".join(fnames[i].split("/")[:-1])
	fnames = np.unique(fnames)
	return(fnames)

# Process data
def processing(flight, plotAlt = False):
	flight.getData()
	flight.standardizeTime()
	flight.mergeData()
	flight.calc_mach()
	flight.extrapolate()
	flight.generateDataSet()
	if plotAlt == True:
		flight.plotDataTime()

# Plot data
def plotData(errorData):
	fit = np.polyfit(errorData["mach"], errorData["CPI"], 1)
	print("Fitting paramaters:\ny = {}x {}".format(fit[0],fit[1]))
	fig, ax = plt.subplots()
	fn = np.poly1d(fit)
	plt.scatter(errorData["mach"],errorData["CPI"], c = errorData["year"], label ="Pressure-Altitude Extrapolation")
	plt.plot(errorData["mach"], fn(errorData["mach"]))
	ax.set_ylabel("Cpi PS_RVSM-P/Q_RVSM []")
	ax.set_xlabel("Indicated Mach number []")
	ax.grid(True)
	ax.legend(loc = "upper left")
	plt.show()

# Main
def main():
	errorData = pd.DataFrame()
	limit = 1000
	for i, filePath in enumerate(sondeFilePaths()):
		if limit > i:
			try:
				flight = Flight(filePath)
				processing(flight)
				errorData = errorData.append(flight.errorData)
			except Exception as e:
				print("Error occured during processing:\n{}".format(e))
				pass	
	errorData["CPI"] = errorData["error"]/errorData["dpressure"]
	errorData = errorData.dropna()
	errorData = errorData[errorData["mach"] > 0.45]
	errorData.to_csv("dropsonde.csv")
	return(errorData)


plotData(main())

