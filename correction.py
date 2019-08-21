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
	for i in range(5,10):
		path = "/media/faamarchive/badcMirror/data/201{}/".format(i)
		fnames += [os.path.join(dirpath, f)
			for dirpath, dirnames, files in os.walk(path)
			for f in fnmatch.filter(files, "faam-dropsonde*[!raw].nc")]
	for i in range(len(fnames)):
		fnames[i] = "/".join(fnames[i].split("/")[:-1])
	fnames = np.unique(fnames)
	return(fnames)
def processing(flight):
	flight.getData()
	flight.standardizeTime()
	flight.mergeData()
	flight.subExtrapolate()
	flight.fullExtrapolate()

# Main
errorData = pd.DataFrame()
limit = 0
for i in sondeFilePaths():
	limit += 1
	if limit > 50:
		break
	try:
		flight = Flight(i)
		processing(flight)
		errorData = errorData.append(flight.errorData)
	except Exception as e:
		print(e)
		pass	

print(errorData.reset_index(drop=True))
plt.scatter(errorData["Altitude"],errorData["subError"], c = "Blue")
plt.scatter(errorData["Altitude"],errorData["fullError"], c = "Red")
plt.show()