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
stepsPerMag = 5

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
	print("sort abs mags")
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


# # separate absolute magnitudes into bins
# def distributeMagnitudes(mags, length, maxM, minM, bins):
# 	print("dist mags")
# 	magBins = np.zeros(bins)
# 	magRange = maxM - minM
# 	binSize = magRange/bins
# 	for mag in mags:
# 		magDif = mag-minM
# 		magBins[math.floor(magDif/binSize)] += 1
# 	return magBins


def effectiveVolume(mags):
	print("effective vol")
	effectiveVols = []
	Vmaxs = []
	for y in range(len(mags)):
		dmax = pow(10, (limMag - mags[y] + 5)/5)
		dmax = dmax / Mega # convert to Mpc
		x = dmax/h
		Omega = 1.e5 #scalar to shift the graph
		effV = Omega * (h**3) * ((.5 * pow(x,2)) + ((x + 1)*math.exp(-1 * x)) - 1)
		effectiveVols.append(effV)

		vol = math.pi * pow((h / 2), 2) * (dmax / 3)
		Vmaxs.append(vol)
	return effectiveVols, Vmaxs

def calcVratios(Veffs, Vmaxs):
	print("vratios")
	result = []
	for i in range(len(Veffs)):
		result.append(Veffs[i]/Vmaxs[i])
	return result

# raw and corrected lminosity function plot for all galaxies
def allGalaxiesPlot(x, lumFunc, correctedLumFunc):
	plt.clf()
	plt.plot(x, lumFunc, "c--", label="Raw Luminosity Function")
	plt.plot(x, correctedLumFunc, "g", label="Corrected Luminosity Function")
	plt.legend(prop={'size': 7})
	plt.yscale('log')
	plt.xlabel("Absolute Magnitude")
	plt.ylabel("#/volume (Mpc^3)")
	plt.title("Luminosity Function for SDSS data")
	plt.savefig("luminosity_function.png")

# corrected luminosity function plots for all galaxies, red and blue separate
def redBluePlot(x, correctedLumFunc, x1, redLumFunc, x2, blueLumFunc):
	print("redblue")
	plt.clf()
	plt.plot(x1, redLumFunc, "r.", label="Corrected Luminosity Function for Red Galaxies")
	plt.plot(x2, blueLumFunc, "b--", label="Corrected Luminosity Function for Blue Galaxies")
	plt.plot(x, correctedLumFunc, "g", label="Corrected Luminosity Function for all Galaxies")
	plt.legend(prop={'size': 7})
	plt.yscale('log')
	plt.xlabel("Absolute Magnitude")
	plt.ylabel("#/volume (Mpc^3)")
	plt.title("Corrected Luminosity Functions for SDSS data")
	plt.savefig("red_blue_all_luminosity.png")

#schecter function
def schecter(M):
	print("schecter")
	alpha = -1.13
	Mstar = -20.8

	result = []
	for m in M:
		part1 = pow(10, -0.4 * (alpha + 1) * m)
		# print(m)
		power = -10**(0.4 * (Mstar - m))
		part2 = math.exp(power)
		result.append(part1 * part2)
	return result

# schecter function over correct luminosity functions
def redBluePlotSchecter(x, correctedLumFunc, x1, redLumFunc, x2, blueLumFunc):
	print("redblue schecter")
	plt.clf()
	plt.plot(x1, redLumFunc, "r.", label="Corrected Luminosity Function for Red Galaxies")
	plt.plot(x2, blueLumFunc, "b--", label="Corrected Luminosity Function for Blue Galaxies")
	plt.plot(x, correctedLumFunc, "g", label="Corrected Luminosity Function for all Galaxies")
	plt.plot(x, schecter(x), "yo", label="Schecter Function, α=-1.13, M*=-20.8")
	plt.legend(prop={'size': 7})
	plt.yscale('log')
	plt.xlabel("Absolute Magnitude")
	plt.ylabel("#/volume (Mpc^3)")
	plt.title("Schecter Function for SDSS data")
	plt.savefig("schecter.png")

# malmquist bias correction
# returns correct luminosity function
def malmquistCorrection(magHist, magBins):
	print("malmquist")
	magBins = magBins.tolist()
	magBins = [round(i, 2) for i in magBins]

	effectiveVols, Vmaxs = effectiveVolume(magBins)
	Vratios = calcVratios(effectiveVols, Vmaxs)
	
	correctedLumFunc = []
	for x in range(len(magHist)):
		correctedLumFunc.append(magHist[x]/effectiveVols[x])

	return correctedLumFunc

# raw and corrected luminosity functions
def luminosityFunction():
	print("luminosity")
	distances = redshiftDistance(zData)
	absoluteMags = absMags(distances, rData)
	print(max(absoluteMags),min(absoluteMags))
	dif = max(absoluteMags)-min(absoluteMags)
	bins = np.arange(math.floor(min(absoluteMags)),math.ceil(max(absoluteMags)),1/stepsPerMag)
	print(bins)
	#need volume to get number density (#/volume in each bin)
	magHist, magBins, magPatches = plt.hist(absoluteMags, bins, histtype='step')
	# magBins = distributeMagnitudes(absoluteMags, math.floor(max(absoluteMags) - min(absoluteMags)), 
							# math.ceil(max(absoluteMags)), math.floor(min(absoluteMags)), len(absoluteMags)*3)
	plt.clf()
	trailingZeroes = 0
	for x in range(len(magHist)-1,0,-1):
		if magHist[x] == 0.:
			trailingZeroes += 1
		else:
			break
	print("mag trail: ", trailingZeroes)
	for x in range(1,len(magHist)-trailingZeroes):
		if magHist[x] == 0.:
			magHist[x] = magHist[x-1]

	magBins = np.delete(magBins, 0)
	magHist = magHist[0:len(magHist)-trailingZeroes]
	magBins = magBins[0:len(magBins)-trailingZeroes]
	print("maghistlen", len(magHist))
	# print("rdata len: " + str(len(rData)), "magHist len: " + str(len(magHist)), "magBins len: " + str(len(magBins)))
	print(len(magHist), magHist, "\n", len(magBins), magBins)
	#solid angle * r^3
	vol = pow(maxR,3) * solidAngle #mpc^3
	lumFunc = [x / vol for x in magHist]
	#correct malmquist bias
	correctedLumFunc = malmquistCorrection(magHist, magBins)
	
	# x = np.arange(min(absoluteMags), max(absoluteMags)-(magBins[1]-magBins[0])*1, magBins[1]-magBins[0])

	#plot for all galaxies
	allGalaxiesPlot(magBins, lumFunc, correctedLumFunc)
	#separated by red and blue galaxies
	redMags, blueMags = sortAbsMags(distances, rData)

	rbins = np.arange(math.floor(min(redMags)),math.ceil(max(redMags)),1/stepsPerMag)
	bbins = np.arange(math.floor(min(blueMags)),math.ceil(max(blueMags)),1/stepsPerMag)
	
	redHist, redBins, redPatches = plt.hist(redMags, rbins, histtype="step")
	blueHist, blueBins, bluePatches = plt.hist(blueMags, bbins, histtype="step")
	redBins = np.delete(redBins, 0)
	blueBins = np.delete(blueBins, 0)
	redHist = redHist[0:len(redHist)-trailingZeroes]
	redBins = redBins[0:len(redBins)-trailingZeroes]
	blueHist = blueHist[0:len(blueHist)-trailingZeroes]
	blueBins = blueBins[0:len(blueBins)-trailingZeroes]
	print(len(redHist), redHist, "\n")#, len(redBins), redBins)
	print(len(blueHist), blueHist,"\n")#, len(blueBins), blueBins)
	plt.clf()
	# redBins = distributeMagnitudes(redMags, math.floor(max(redMags) - min(redMags)), 
					# math.ceil(max(redMags)), math.floor(min(redMags)), len(redMags)*3)
	# blueBins = distributeMagnitudes(blueMags, math.floor(max(blueMags) - min(blueMags)), 
					# math.ceil(max(blueMags)), math.floor(min(blueMags)), len(blueMags)*3)
	
	trailingZeroes = 0
	for x in range(len(redHist)-1,0,-1):
		if redHist[x] == 0.:
			trailingZeroes += 1
		else:
			break
	print("red trail: ", trailingZeroes)
	for x in range(1,len(redHist)-trailingZeroes):
		if redHist[x] == 0.:
			redHist[x] = redHist[x-1]

	trailingZeroes = 0
	for x in range(len(blueHist)-1,0,-1):
		if blueHist[x] == 0.:
			trailingZeroes += 1
		else:
			break
	print("blue trail: ", trailingZeroes)
	for x in range(1,len(blueHist)-trailingZeroes):
		if blueHist[x] == 0.:
			blueHist[x] = blueHist[x-1]
	redLumFunc = malmquistCorrection(redHist, redBins)
	blueLumFunc = malmquistCorrection(blueHist, blueBins)

	# x1 = np.arange(min(redMags), max(redMags)-(redBins[1]-redBins[0])+.01, redBins[1]-redBins[0])
	# x2 = np.arange(min(blueMags), max(blueMags)-(blueBins[1]-blueBins[0])+.01, blueBins[1]-blueBins[0])
	x1 = np.arange(math.floor(min(redMags)),math.ceil(max(redMags)),1/4)
	x2 = np.arange(math.floor(min(blueMags)),math.ceil(max(blueMags)),1/4)

	#plot for red, blue, all
	redBluePlot(magBins, correctedLumFunc, redBins, redLumFunc, blueBins, blueLumFunc)

	#plot for red, blue, all, schecter
	redBluePlotSchecter(magBins, correctedLumFunc, redBins, redLumFunc, blueBins, blueLumFunc)
	
luminosityFunction()
#ugDistPlot()

hdulist.close()