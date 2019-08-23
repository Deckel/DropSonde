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
	#flight.altitudeExtrapolate()
	#flight.plotDataTime()
	#flight.plotDataAlt()

# Main

errorData = pd.DataFrame()
limit = 0
for i in sondeFilePaths():
	limit += 1
	if limit > 10:
		break
	try:
		flight = Flight(i)
		processing(flight)
		errorData = errorData.append(flight.errorData)
	except Exception as e:
		print(e)
		pass	
errorData = errorData.reset_index(drop=True)

fig, ax = plt.subplots()

plt.scatter(errorData["Altitude"],errorData["subError"], c = "Blue", label = "Sub-sample Temporal Extrapolation")
sc = plt.scatter(errorData["Altitude"],errorData["fullError"], c = "Red", label= "Full Temporal Extrapolation")
#plt.scatter(errorData["Altitude"],errorData["altError"], c = "Green", label ="Pressure-Altitude Extrapolation")
#ax.set_ylabel("Error /hPa Expected - Observed")
#ax.set_xlabel("Altitude /m")
#ax.grid(True)
annot = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points")
annot.set_visible(False)

def update_annot(ind):
	pos = sc.get_offsets()[ind["ind"][0]]
	annot.xy = pos
	text = "{}".format(" ".join([errorData["flight"][n] for n in ind["ind"]]))

	annot.set_text(text)

def hover(event):
	vis = annot.get_visible()
	if event.inaxes == ax:
		cont, ind = sc.contains(event)
		if cont:
			update_annot(ind)
			annot.set_visible(True)
			fig.canvas.draw_idle()
		else:
			if vis:
				annot.set_visible(False)
				fig.canvas.draw_idle()

# for i in errorData:
# 	plt.annotate(errorData["flight"], (errorData["Altitude"],errorData["subError"]))
fig.canvas.mpl_connect("motion_notify_event", hover)
plt.show()