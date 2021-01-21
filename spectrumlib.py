import numpy, os
import json
import scipy.interpolate
import scipy.integrate

class spectrumObject:
	def __init__(self):
		self.wavelengths = []
		self.flux = []
		self.fluxErrors = []
		self.length = 0 
		self.wavelengthRange = (0, 0)
		self.name = 'unknown'
		self.loadedFromFilename = 'unknown'
		self.wavelengthUnits = 'unknown'
		self.fluxUnits = 'unknown'
		self.objectName = 'unknown'

	def __str__(self):
		retStr = ""
		retStr+= "%s\n"%(self.name)
		retStr+= "Wavelength range: (%4.2f, %4.2f) %s\n"%(self.wavelengthRange[0], self.wavelengthRange[1], self.wavelengthUnits)
		retStr+= "Flux units: %s"%(self.fluxUnits)
		return retStr
		
	def getProperty(self, property):
		try:
			data = getattr(self, property)
			return data
		except AttributeError:
			return None
		  
		
	def setData(self, wavelengths, flux, errors = []):
		if len(wavelengths) != len(flux):
			return -1
		self.wavelengths = []
		self.flux = []
		self.fluxErrors = []
		
		if len(errors) == 0:
			for w, f in zip(wavelengths, flux):
				self.wavelengths.append(w)
				self.flux.append(f)
		else:
			for w, f, e in zip(wavelengths, flux, errors):
				self.wavelengths.append(w)
				self.flux.append(f)
				self.fluxErrors.append(e)
		
		self.length = len(wavelengths)
		self.wavelengthRange = (min(wavelengths), max(wavelengths))
		
		return self.length
		
	
	def fitPoly(self, degree=3, mask = []):
		# Apply the mask
		inputSpectrum = spectrumObject()
		inputSpectrum.setData(self.wavelengths, self.flux, self.fluxErrors)
		for m in mask:
			inputSpectrum.snipWavelengthRange(m[0], m[1])
		x = [w for w in inputSpectrum.wavelengths]
		y = [f for f in inputSpectrum.flux]
		fit = numpy.poly1d(numpy.polyfit(x, y, degree))
		fittedFlux = fit(self.wavelengths)
		return fittedFlux

	def rebin(self, factor = 5):
		wavelengths = self.wavelengths
		fluxes = self.flux
		fluxErrors = self.fluxErrors
		steps = int(len(wavelengths) / factor)
		newLambdas = []
		newFluxes = []
		for n in range(steps):
			startIndex = n*factor
			lambdas = []
			sampleFluxes = []

			for f in range(factor):
				lambdas.append(wavelengths[startIndex+f])
				sampleFluxes.append(fluxes[startIndex+f])
			newLambda = numpy.mean(lambdas)
			newFlux = numpy.mean(sampleFluxes)
			newLambdas.append(newLambda)
			newFluxes.append(newFlux)
			print(lambdas, ':', newLambda)
			print(sampleFluxes, ':', newFlux)
		
		self.wavelengths = newLambdas
		self.flux = newFluxes
		print(len(self.wavelengths), len(self.flux))

	def sortData(self):
		wavelengths = self.wavelengths
		fluxes = self.flux
		list1, list2 = zip(*sorted(zip(wavelengths, fluxes)))
		self.wavelengths = list1
		self.flux = list2

	def resample(self, sampleWavelengths):
		startWavelength = min(sampleWavelengths)
		endWavelength = max(sampleWavelengths)
		self.trimWavelengthRange(startWavelength, endWavelength)
		# print "num points:", len(self.wavelengths)
		spline = scipy.interpolate.splrep(self.wavelengths, self.flux, s=0)
		sampleFlux = scipy.interpolate.splev(sampleWavelengths, spline, der=0)
		spline = scipy.interpolate.splrep(self.wavelengths, self.fluxErrors, s=0)
		self.fluxErrors = scipy.interpolate.splev(sampleWavelengths, spline, der=0)
		self.flux = sampleFlux
		self.wavelengths = sampleWavelengths
		return sampleFlux

	def removeNegatives(self):
		newFlux = []
		newFluxErrors = []
		for w, f, fe in zip(self.wavelengths, self.flux, self.fluxErrors):
			if f<0: 
				f = 0
				fe = 0
			newFlux.append(f)
			newFluxErrors.append(fe)
		self.flux = newFlux
		self.fluxErrors = newFluxErrors

	def removeZeros(self):
		newWavelengths = []
		newFlux = []
		newFluxErrors = []
		for w, f, fe in zip(self.wavelengths, self.flux, self.fluxErrors):
			if f==0: continue 
			newWavelengths.append(w)
			newFlux.append(f)
			newFluxErrors.append(fe)
		self.wavelengths = newWavelengths
		self.flux = newFlux
		self.fluxErrors = newFluxErrors



	def writeCSV(self, filename):
		outputfile = open(filename, 'w')
		outputfile.write("wavelength, flux\n")
		for w, f in zip(self.wavelengths, self.flux):
			outputfile.write("%f, %f\n"%(w, f))
		outputfile.close()
		
	def snipWavelengthRange(self, lower, upper):
		""" Removes a section from the spectrum """
		newWavelengths = []
		newFlux = []
		newFluxErrors = []
		if lower>=upper: return self.length
			
		for w, f, fe in zip(self.wavelengths, self.flux, self.fluxErrors):
			if (w<lower) or (w>upper):
				newWavelengths.append(w)
				newFlux.append(f)
				newFluxErrors.append(fe)
		self.length = len(newWavelengths)
		self.wavelengths = newWavelengths
		self.flux = newFlux
		self.fluxErrors = newFluxErrors
		self.wavelengthRange = (min(self.wavelengths), max(self.wavelengths))
			
		return self.length		
	
	def integrate(self, wavelengthrange = (-1, -1)):
		""" Integrates under the spectrum between two wavelength limits. Defaults to all of the spectrum """
		if wavelengthrange[0]!=-1:
			wavelengths, fluxes = self.getSubsetByWavelength(wavelengthrange[0], wavelengthrange[1])
		else: 
			wavelengths, fluxes = (self.wavelengths, self.flux)
		
		total = scipy.integrate.simps(fluxes, wavelengths)
		return total

	def divide(self, constant):
		""" Divides the spectrum by a constant value """
		newFlux = []
		for w, f in zip(self.wavelengths, self.flux):
			newFlux.append(f / constant)
		self.flux = newFlux
		return 

	def divideArray(self, fluxValues):
		""" Divides each value of the flux by the corresponding value in the array """
		if len(self.flux) != len(fluxValues): return
		self.flux = self.flux/fluxValues
		self.fluxErrors = self.fluxErrors/fluxValues


	def subtractSpectrum(self, subtractSpectrum):
		if len(self.wavelengths)!=len(subtractSpectrum.wavelengths):
			print("Can't subtract spectra of different lengths.")
			return
		newFlux = []
		for (aw, af, bw, bf) in zip(self.wavelengths, self.flux, subtractSpectrum.wavelengths, subtractSpectrum.flux):
			f = af - bf
			#print af, bf, f
			newFlux.append(f)
		self.flux = newFlux
		return
		
	def trimWavelengthRange(self, lower, upper):
		""" Trims out the lower and upper portions of the spectrum """ 
		newWavelengths = []
		newFlux = []
		newFluxErrors = []
		if lower>=upper: return self.length
		
		for w, f, fe in zip(self.wavelengths, self.flux, self.fluxErrors):
			if (w>lower) and (w<upper):
				newWavelengths.append(w)
				newFlux.append(f)
				newFluxErrors.append(fe)
		self.length = len(newWavelengths)
		self.wavelengths = newWavelengths
		self.flux = newFlux
		self.fluxErrors = newFluxErrors
		self.wavelengthRange = (min(self.wavelengths), max(self.wavelengths))
			
		return self.length		
		
	def getWavelengths(self):
		return self.wavelengths
		
	def convertFluxes(self):
		print("Current fluxUnits are:", self.fluxUnits)
		print("Converting from mJy to erg/s/cm^2/A")
		c = 3E8
		newFlux = []
		newFluxErrors = []
		for w, f, fe in zip(self.wavelengths, self.flux, self.fluxErrors):
			flambda = f * 1E-16 * c / (w*w)
			error = fe * 1E-16 * c / (w*w)
			newFlux.append(flambda)
			newFluxErrors.append(error)
		self.flux = newFlux 
		self.fluxErrors = newFluxErrors
		self.fluxUnits = "erg/s/cm^2/A"
		
	def getFlux(self):
		return self.flux
	
	def getFluxErrors(self):
		return self.fluxErrors
		
	def getNearestFlux(self, wavelength):
		minDistance = self.wavelengthRange[1]
		for w, f in zip (self.wavelengths, self.flux):
			distance = abs(w - wavelength)
			if distance < minDistance:
				minDistance = distance
				result = f
		return result
		
	def getSubsetByWavelength(self, lower, upper):
		newWavelengths = []
		newFlux = []
		for w, f in zip(self.wavelengths, self.flux):
			if (w>=lower) and (w<=upper):
				newWavelengths.append(w)
				newFlux.append(f)
		return (newWavelengths, newFlux)

	def writeToJSON(self, filename, clobber=True):
		object = {}
		for key in self.__dict__.keys():
			data = getattr(self, key)
			# print key, type(data)
			if type(data)==numpy.float32:
				data = float(data)
			if type(data)==numpy.ndarray:
				data = numpy.array(data).tolist()
			if type(data)==list:
				data = numpy.array(data).tolist()
			object[key] = data
			
		if not clobber:
			if os.path.exists(filename): 
				print("Warning: file %s exists... skipping"%filename)
				return
		outputfile = open(filename, 'w')
		json.dump(object, outputfile, indent=4)
		outputfile.close()

	def writeToMollyText(self, filename):
		outputFile = open(filename, 'wt')
		for w, f, fe in zip(self.wavelengths, self.flux, self.fluxErrors):
			outputFile.write("%f\t%f\t%f\n"%(w, f, fe))
		outputFile.close()
	
	def loadFromJSON(self, filename):
		inputfile = open(filename, "r")
		jsonObject = json.load(inputfile)
		for key in jsonObject.keys():
			keyString = str(key)
			value = jsonObject[key]
			if isinstance(value, (str)): 
				value = str(value)
			if isinstance(value, (list)):
				value = numpy.array(value)
			setattr(self, key, value)
		inputfile.close()
		self.loadedFromFilename = filename
		
	def loadFromSLOAN(self, filename):
		inputfile = open(filename, "rt")
		print("SDSS filename is:", filename)
		self.objectName = filename[:23]
		print("Target name is:", self.objectName)
		
		for line in inputfile:
			fields = line.strip().split()
			self.wavelengths.append(float(fields[0]))
			self.flux.append(float(fields[1]))
			self.fluxErrors.append(float(fields[2]))
			
		inputfile.close() 
		
	def extractFITSHeaders(self, headers):
		for h in headers:
			if h== "OBJECT": 
				self.objectName = headers[h]
				print("Object name: %s"%self.objectName)
			if h== "HJD": 
				self.HJD = headers[h]
				print("HJD: %f"%self.HJD)
			if h== "BUNIT": 
				self.fluxUnits = headers[h]
				print("Flux units: %s"%self.fluxUnits)
			if h== "DATE-OBS": 
				parts = headers[h].split('-')
				try:
					self.year = int(parts[0])
					self.month = int(parts[1])
					self.day = int(parts[2])
				except Exception as e:
					print(e)
				print("Date observed: %04d/%02d/%02d"%(self.year, self.month, self.day))
			if h== "TELESCOP": 
				self.telescope = headers[h]
				print("Telescope: %s"%self.telescope)
			if h== "AIRMASS": 
				self.airmass = float(headers[h])
				print("Airmass: %f"%self.airmass)
			if h== "EXPTIME": 
				self.exptime = float(headers[h])
				print("Exposure time: %f s"%self.exptime)
			if h=="VHELIO":
				self.vhelio = float(headers[h])
				print("Heliocentric velocity: %f km/s"%self.vhelio)

		self.name = "%s-%f"%(self.objectName, self.HJD)

		return
		
	def parseHeaderInfo(self, headers):
		
		try: 
			self.objectName = headers['Object']
			self.telescope = headers['Telescope']
			self.ra = headers['RA']
			self.dec = headers['Dec']
			self.UTC = headers['UTC']
			self.RJD = headers['RJD']
			self.equinox = headers['Equinox']
			self.Vearth = headers['Vearth']
			self.hourangle = headers['Hour angle']
			self.longitude = headers['Longitude']
			self.latitude = headers['Latitude']
			self.siderial = headers['Sidereal time']
			self.site = headers['Site']
			self.day = headers['Day']
			self.month = headers['Month']
			self.year = headers['Year']
			self.dwell = headers['Dwell']
			self.HJD = headers['HJD']
			self.airmass = headers['Airmass']
			self.galLatitude = headers['Gal latitude']
			self.galLongitude = headers['Gal longitude']
			self.extractPosition = headers['Extract position']
			self.comment = str(headers['COMMENT'])
			self.phase = float(header['Orbital phase'])
		except Exception as e:
			print(e)
		
		return self.objectName

	
	def appendDataAtNewWavelengths(self, wavelengths, flux):
		currentWavelengthRange = self.wavelengthRange
		newWavelengthsRange = (min(wavelengths), max(wavelengths))
		print("Current wavelength range:", currentWavelengthRange)
		print("New data wavelength range:", newWavelengthsRange)
		if currentWavelengthRange[1] > newWavelengthsRange[0]:
			print("There *IS* an overlap in wavelengths")
			print("First data loaded takes preference")
			startWavelength = currentWavelengthRange[1]
			print("Starting to add new data from above: ", startWavelength)
			for w, f in zip(wavelengths, flux):
				if w>startWavelength:
					self.wavelengths.append(w)
					self.flux.append(f)
					
		self.length = len(self.wavelengths)
		self.wavelengthRange = (min(wavelengths), max(wavelengths))
		
		return self.length
		
	def appendNewData(self, newSpectrum):
		print("HJD of existing spectrum: %7.7f\nHJD of spectrum to be added: %7.7f"%(self.HJD, newSpectrum.HJD))
		timeDifferenceSeconds = 86400. * (self.HJD - newSpectrum.HJD)
		print("Time difference is %f seconds."%timeDifferenceSeconds)
		
		currentWavelengthRange = self.wavelengthRange
		wavelengths = newSpectrum.getWavelengths()
		flux = newSpectrum.getFlux()
		newWavelengthsRange = (min(wavelengths), max(wavelengths))
		print("Current wavelength range:", currentWavelengthRange)
		print("New data wavelength range:", newWavelengthsRange)
		
		# Append data at the end of the current data
		if currentWavelengthRange[0] < newWavelengthsRange[0]:
			if currentWavelengthRange[1] > newWavelengthsRange[0]:
				print("There *IS* an overlap in wavelengths")
				print("First data loaded takes preference")
				startWavelength = currentWavelengthRange[1]
				print("Starting to add new data from above: ", startWavelength)
				for w, f in zip(wavelengths, flux):
					if w>startWavelength:
						self.wavelengths.append(w)
						self.flux.append(f)
						
		# Add new data to the beginning of the old data
		if currentWavelengthRange[0] > newWavelengthsRange[0]:
			if currentWavelengthRange[1] > newWavelengthsRange[0]:
				print("There *IS* an overlap in wavelengths")
				print("First data loaded takes preference")
				endWavelength = currentWavelengthRange[0]
				newFlux = []
				newWavelengths = []
				print("Starting to add new data until: ", endWavelength)
				for w, f in zip(wavelengths, flux):
					if w<endWavelength:
						newWavelengths.append(w)
						newFlux.append(f)
				for w, f in zip(self.wavelengths, self.flux):
					newWavelengths.append(w)
					newFlux.append(f)
			self.wavelengths = newWavelengths
			self.flux = newFlux
		
						
		self.length = len(self.wavelengths)
		self.wavelengthRange = (min(wavelengths), max(wavelengths))
		
		return self.wavelengthRange

