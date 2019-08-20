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
	return(fnames)


def processing(flight):
	flight.getData()
	flight.standardizeTime()
	flight.mergeData()
	flight.extrapolate()


errorData = pd.DataFrame(columns = ["error", "ALT_GIN"])

# Main
for i in sondeFilePaths():
	try:
		flight = Flight(i)
		processing(flight)
	except:
		pass
