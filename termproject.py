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


# Step 1, 2: Clean the data
# Remove entries with null values (-9999)
# Remove entries with redshifts greater than 0.3 (too far)

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
	if zData[i] <= 0 or zData[i] > .3:
		gData = np.delete(gData, i)
		zData = np.delete(zData, i)
		uData = np.delete(uData, i)
		i -= 1
	i += 1


# Step 3: u-g color distribution plot

def ugDistPlot():
	plt.bar(np.arange(12.8,21.8, 1), uBars, width=0.25, align="edge", label="u", color="b", hatch="//")
	plt.bar(np.arange(13.2,22.2, 1), gBars, width=0.25, label="g", color="g")
	plt.xticks(np.arange(13,22))
	plt.xlabel("Apparent Magnitude")
	plt.ylabel("Frequency")
	plt.legend()

	plt.title("Distribution of U-G color in sdss data")
	#color is 18.5

	plt.savefig("histogram.png")

# Step 4: Ram Luminosity Function

# compute distance via redshift (d = cz/H0)
def redshiftDistance(x):
	#H0 used is 73.8
	H0 = 73.8
	y = [((299792 * z)/H0) * 1000000 for z in x]
	return y

# compute absolute magnitude via distance and apparent magnitude (distance modulus)
def absMags(x, u):
	absMagnitudes = []
	j = 0
	while j < len(x):
		result = u[j] - ((5 * math.log(x[j])) - 5)
		absMagnitudes.append(result)
		j += 1
	return absMagnitudes


# separate absolute magnitudes into bins
def distributeMagnitudes(mags):
	magBins = np.zeros(58)
	for mag in mags:
		magBins[math.floor(mag) + 24] += 1
	return magBins


def calcEffectiveVolume():
	#from range -82 to -24
	#lim mag is 17.6
	limMag = 17.6
	h = .3
	Omega = 1

	effectiveVols = []
	Vmaxs = []
	for y in range(-82, -23):
		dmax = pow(10, (limMag - y + 5)/5)
		dmax = dmax / 1000000
		x = dmax/h
		effV = Omega * (h**3) * ((.5 * pow(x,2)) + ((x + 1)*math.exp(-1 * x)) - 1)
		effectiveVols.append(effV)

		vol = math.pi * pow((h / 2), 2) * (dmax / 3)
		Vmaxs.append(vol)
	return effectiveVols, Vmaxs

def calcVratios(Veffs, Vmaxs):
	result = []
	for i in range(len(Veffs)):
		result.append(Veffs[i]/Vmaxs[i])
	return result


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
	print(magBins)
	vol = 2338967.98142
	lumFunc = [x / vol for x in magBins]

	#correct malmquist bias
	#using h=300pc, Omega = 1
	effectiveVols, Vmaxs = calcEffectiveVolume()
	Vratios = calcVratios(effectiveVols, Vmaxs)
	#print(effectiveVols)
	correctedLumFunc = []
	for x in range(len(magBins)):
		correctedLumFunc.append(magBins[x]/Vratios[x])


	#plt.plot(np.arange(-82, -24, 1), lumFunc, "r")
	plt.plot(np.arange(-82, -24, 1), correctedLumFunc, "b")
	plt.xlabel("Absolute Magnitude")
	plt.ylabel("#/volume (Mpc^3)")
	plt.title("Raw Luminosity Function for SDSS data")
	plt.savefig("rawlum.png")

rawLum()

hdulist.close()