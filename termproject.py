import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.io import fits
import math

zmax = 0.1
zmin = 1 * 1.e-2
H0 = 73.8 # km /s / Mpc
c = 299792. #km /s
Mega = 1000000. # Mpc to pc
solidAngle = 0.00387
maxR = (c * zmax) / H0
limMag = 17.6
h = 3.0/Mega # Mpc
colorDifference = 1.565

hdulist = fits.open("sdss.fits")

hdulist.info()

uData = hdulist[1].data['u']
gData = hdulist[1].data['g']
zData = hdulist[1].data['z']
rData = hdulist[1].data['r']

#list to hold u-g color
uMinusg = []


# Step 1, 2: Clean the data
# Remove entries with null values (-9999)
# Remove entries with redshifts greater than 0.3 (too far)

for i in range(len(rData) - 1):
	if rData[i] == -9999.0:
		gData = np.delete(gData, i)
		zData = np.delete(zData, i)
		uData = np.delete(uData, i)
		rData = np.delete(rData, i)
		i -= 1
i = 0
while i < len(zData):
	if zData[i] <= zmin or zData[i] > zmax:
		gData = np.delete(gData, i)
		zData = np.delete(zData, i)
		uData = np.delete(uData, i)
		rData = np.delete(rData, i)
		i -= 1
	i += 1

i = 0
while i < len(uData):
	if uData[i] == -9999.0 or gData[i] == -9999.0:
		uData = np.delete(uData, i)
		zData = np.delete(zData, i)
		gData = np.delete(gData, i)
		rData = np.delete(rData, i)
		i -= 1
	else:
		color = uData[i] - gData[i]
		if color < 3:
			uMinusg.append(uData[i] - gData[i])
		else:
			uData = np.delete(uData, i)
			zData = np.delete(zData, i)
			gData = np.delete(gData, i)
			rData = np.delete(rData, i)
			i -= 1
	i += 1

# Step 3: u-g color distribution plot

def ugDistPlot():
	n, bins, patches = plt.hist(uMinusg, bins = 100)
	plt.axis([0,3,0,2000])
	plt.plot(bins)
	plt.axvline(colorDifference, c='r', label="Red Blue Galaxy Split")
	plt.legend()
	plt.xlabel("Color")
	plt.ylabel("Frequency (Hz)")
	plt.title("Distribution of U-G color in sdss data")
	plt.savefig("histogram.png")

# Step 4: Raw Luminosity Function

# compute distance via redshift (d = cz/H0)
def redshiftDistance(x):
	y = [((c * z)/H0) * Mega for z in x] # 1000000 to get pc
	return y

# compute absolute magnitude via distance and apparent magnitude (distance modulus)
# dists is list of distance
# mags is list  of apparent magnitudes
def absMags(dists, mags):
	absMagnitudes = []
	j = 0
	while j < len(dists):
		result = mags[j] - ((5 * np.log10(dists[j])) - 5)
		absMagnitudes.append(result)
		j += 1
	return absMagnitudes

# Returns two lists of absolute magnitudes, one for red galaxies, one for blue galaxies
# compute absolute magnitude via distance and apparent magnitude (distance modulus)
# dists is a list of distance
# mags is a list of apparent magnitude
def sortAbsMags(dists, mags):
	red = []
	blue = []
	j = 0
	while j < len(dists):
		magnitude = mags[j] - ((5 * np.log10(dists[j])) - 5)
		if uMinusg[j] > colorDifference:
			red.append(magnitude)
		else:
			blue.append(magnitude)
		j += 1
	return red, blue


# separate absolute magnitudes into bins
def distributeMagnitudes(mags, length, minV):
	magBins = np.zeros(length + 1)
	for mag in mags:
		magBins[math.floor(mag) + minV] += 1
	return magBins


def effectiveVolume(mags):
	print("m", mags)
	effectiveVols = []
	Vmaxs = []
	for y in range(len(mags)):
		dmax = pow(10, (limMag - mags[y] + 5)/5)
		dmax = dmax / Mega
		x = dmax/h
		Omega = 1.e6 #scalar to shift the graph
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

# raw and corrected lminosity function plot for all galaxies
def allGalaxiesPlot(x, lumFunc, correctedLumFunc):
	plt.plot(x, lumFunc, "r", label="Raw Luminosity Function")
	plt.plot(x, correctedLumFunc, "b", label="Corrected Luminosity Function")
	plt.legend(prop={'size': 7})
	plt.yscale('log')
	plt.xlabel("Absolute Magnitude")
	plt.ylabel("#/volume (Mpc^3)")
	plt.title("Raw Luminosity Function for SDSS data")
	plt.savefig("luminosity_function.png")

# corrected luminosity function plots for all galaxies, red and blue separate
def redBluePlot(x, correctedLumFunc, x1, redLumFunc, x2, blueLumFunc):
	plt.plot(x, correctedLumFunc, "g", label="Corrected Luminosity Function for all Galaxies")
	plt.plot(x1, redLumFunc, "r", label="Corrected Luminosity Function for Red Galaxies")
	plt.plot(x2, blueLumFunc, "b", label="Corrected Luminosity Function for Blue Galaxies")
	plt.legend(prop={'size': 7})
	plt.yscale('log')
	plt.savefig("red_blue_all_luminosity.png")

# malmquist bias correction
# returns correct luminosity function
def malmquistCorrection(absoluteMags, magBins):
	absMagnitudesList = []
	absMagsRange = range(math.floor(max(absoluteMags)), math.floor(min(absoluteMags)), -1)
	for i in absMagsRange:
		absMagnitudesList.insert(0, i)
	
	effectiveVols, Vmaxs = effectiveVolume(absMagnitudesList)
	Vratios = calcVratios(effectiveVols, Vmaxs)
	
	correctedLumFunc = []
	for x in range(len(magBins)):
		correctedLumFunc.append(magBins[x]/effectiveVols[x])

	return correctedLumFunc

# raw and corrected luminosity functions
def luminosityFunction():
	distances = redshiftDistance(zData)
	absoluteMags = absMags(distances, rData)

	#need volume to get number density (#/volume in each bin)
	magBins = distributeMagnitudes(absoluteMags, math.floor(max(absoluteMags) - min(absoluteMags)), abs(math.floor(min(absoluteMags))) - 1)
	
	#solid angle * r^3
	vol = pow(maxR,3) * solidAngle #mpc^3
	lumFunc = [x / vol for x in magBins]

	#correct malmquist bias
	correctedLumFunc = malmquistCorrection(absoluteMags, magBins)
	
	x = np.arange(min(absoluteMags), max(absoluteMags), 1)

	#plot for all galaxies
	#allGalaxiesPlot(x, lumFunc, correctedLumFunc)

	#separated by red and blue galaxies
	redMags, blueMags = sortAbsMags(distances, rData)
	redBins = distributeMagnitudes(redMags, math.floor(max(redMags) - min(redMags)), abs(math.floor(min(redMags))) - 1)
	blueBins = distributeMagnitudes(blueMags, math.floor(max(blueMags) - min(blueMags)), abs(math.floor(min(blueMags))) - 1)
	
	redLumFunc = malmquistCorrection(redMags, redBins)
	blueLumFunc = malmquistCorrection(blueMags, blueBins)

	x1 = np.arange(min(redMags), max(redMags), 1)
	x2 = np.arange(min(blueMags), max(blueMags), 1)

	#plot for red, blue, all
	redBluePlot(x, correctedLumFunc, x1, redLumFunc, x2, blueLumFunc)
	
luminosityFunction()
#ugDistPlot()

hdulist.close()