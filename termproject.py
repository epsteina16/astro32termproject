import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.io import fits
import math

hdulist = fits.open("sdss.fits")

hdulist.info()

uData = hdulist[1].data['u']
gData = hdulist[1].data['g']
zData = hdulist[1].data['z']


uBars = np.zeros(9)
gBars = np.zeros(9)


for i in range(len(uData) - 1):
	if uData[i] == -9999.0:
		uData = np.delete(uData, i)
		zData = np.delete(zData, i)
		gData = np.delete(gData, i)
		i -= 1
	else:
		bottom = math.floor(uData[i])
		if (bottom < 22 and bottom > 12):
			uBars[bottom - 13] += 1

for i in range(len(gData) - 1):
	if gData[i] == -9999.0:
		gData = np.delete(gData, i)
		zData = np.delete(zData, i)
		uData = np.delete(uData, i)
		i -= 1
	else:
		bottom = math.floor(gData[i])
		if (bottom < 22 and bottom > 12):
			gBars[bottom - 13] += 1
i = 0
while i < len(zData):
	#limit redshifts to .3 - z << 1 - "arbitrary as fuck"
	if zData[i] <= 0 or zData[i] > .3:
		gData = np.delete(gData, i)
		zData = np.delete(zData, i)
		uData = np.delete(uData, i)
		i -= 1
	i += 1


def bar1():
	plt.bar(np.arange(12.8,21.8, 1), uBars, width=0.25, align="edge", label="u", color="b", hatch="//")
	plt.bar(np.arange(13.2,22.2, 1), gBars, width=0.25, label="g", color="g")
	plt.xticks(np.arange(13,22))
	plt.xlabel("Apparent Magnitude")
	plt.ylabel("Frequency")
	plt.legend()

	plt.title("Distribution of U-G color in sdss data")
	#color is 18.5

	plt.savefig("histogram.png")

def redshiftDistance(x):
	#H0 used is 73.8
	H0 = 73.8
	y = [((299792 * z)/73.8) * 1000000 for z in x]
	return y

def absMags(x, u):
	absMagnitudes = []
	j = 0
	while j < len(x):
		result = u[j] - ((5 * math.log(x[j])) - 5)
		absMagnitudes.append(result)
		j += 1
	return absMagnitudes

def distributeMagnitudes(mags):
	magBins = np.zeros(58)
	for mag in mags:
		magBins[math.floor(mag) + 24] += 1
	return magBins

def calcEffectiveVolume():


def rawLum():
	#RA range 145 to 236 degrees
	#Dec range -2 to 2 degrees
	#r is 1219.5122 Mpc
	#volume calculated for survey is 2338967.98142
	distances = redshiftDistance(zData)
	maxDistance = max(distances)
	absoluteMags = absMags(distances, uData)

	#need volume to get number density (#/volume in each bin)
	magBins = distributeMagnitudes(absoluteMags)
	vol = 2338967.98142
	lumFunc = [x / vol for x in magBins]

	#correct malmquist bias
	#using h=300pc, Omega = 1
	


	plt.plot(np.arange(-82, -24, 1), lumFunc, "r")
	plt.xlabel("Absolute Magnitude")
	plt.ylabel("#/volume (Mpc^3)")
	plt.title("Raw Luminosity Function for SDSS data")
	plt.savefig("rawlum.png")

rawLum()

hdulist.close()