import numpy
import json
import generallib
import datetimelib
import astropy, datetime
import sys, os

class target():
	def __init__(self, id):
		self.id = id
		self.MJD = []
		self.mag = []
		self.err = []
		self.ephemeris = None
		self.data = []
		self.hasEphemeris = False
		self.telescope = "CSS"
		self.fluxColumn = ""
		self.dateColumn = ""
		self.fluxErrorColumn = ""
		self.debug = False
		
	def appendData(self, dataDict):
		self.data.append(dataDict)
	
	def getColumn(self, columnName):
		return [d[columnName] for d in self.data]	
		
	def loadEphemeris(self, filename="none"):
		# Look in the local directory for a file called 'id'-ephem.dat and load it
		if filename=='none': filename = self.id + "-ephem.dat"
		if os.path.exists(filename):
			self.ephemeris = datetimelib.ephemerisObject()
			self.ephemeris.loadFromFile(filename)
			self.hasEphemeris = True
			return True
		
		return False

	def trimLargeErrors(self, limit = 0.25):
		trimmedData = []
		errorColumn = self.fluxErrorColumn
		for d in self.data:
			print(d)
			if not abs(d[errorColumn]) > limit:
				trimmedData.append(d)
		print("Trimmed %d points out of %d."%(len(self.data) - len(trimmedData), len(self.data)))
		self.data = trimmedData

	def calcPhase(self):
		for d in self.data:
			HJD = d['HJD']
			phase = self.ephemeris.getPhase(HJD)
			d['phase'] = phase

	def phaseBin(self, numBins = 100):
		bins = []
		phases = []
		errors = []
		stddevs = []
		for bin in range(numBins):
			if self.debug: print("bin #: %d"%bin)
			
			# Find all the datapoints that lie in this bin
			phaseStart = bin / numBins
			phaseEnd = (bin+1) / numBins
			midPhase = (phaseStart+phaseEnd)/2.
			self.debug: print("phase range: %f to %f"%(phaseStart, phaseEnd))
			binData = []
			for d in self.data:
				if (d['phase']>=phaseStart) and (d['phase']<phaseEnd):
					binData.append(d)
			if self.debug: print("%d points in this phase bucket"%len(binData))
			binNumbers = [d[self.fluxColumn] for d in binData]
			print(binNumbers)
			if len(binData)==0:
				print("WARNING: No data in phase range: %f to %f"%(phaseStart, phaseEnd))
				mean = None
				error = None
				continue
			if len(binData)==1:
				mean = binData[0][self.fluxColumn]
				error = binData[0][self.fluxErrorColumn]
			if len(binData)>1:
				sum = 0
				errorsum = 0
				N = len(binData)
				for d in binData:
					sum+= d[self.fluxColumn] / d[self.fluxErrorColumn]
					errorsum+= 1 / d[self.fluxErrorColumn]
				mean = sum/errorsum
				squ_sum = 0
				for d in binData:
					squ_sum+= (d[self.fluxColumn]-mean)**2 / d[self.fluxErrorColumn]
				stddev = numpy.sqrt( squ_sum/ N / errorsum)
				error = numpy.sqrt( squ_sum / N / (N-1) / errorsum)
				print("Weighted mean: %f [%f]"%(mean, error))

				straightMean = numpy.mean(binNumbers)
				straightStdDev = numpy.std(binNumbers)
				sterror = straightStdDev / numpy.sqrt(len(binNumbers))
				print("Straight mean: %f [%f] (%f)"%(straightMean, straightStdDev, sterror))	

			bins.append(mean)
			phases.append(midPhase)
			errors.append(error)
			stddevs.append(stddev)
		if self.debug: print(phases, bins)
		self.bins = bins
		self.phases = phases
		self.errors = errors
		self.stddevs = stddevs
		

	def setHJDs(self, MJD, HJD):
		keys = [d['MJD'] for d in self.data]
		dates = zip(MJD, HJD)
		for index, d in enumerate(dates):
			self.data[index]['HJD'] = d[1]

	def computeHJDfromJD(self):
		if self.hasEphemeris:
			print(self.id, self.ephemeris)
			JD = self.getColumn('JD')
			correctHelio = datetimelib.heliocentric()
			correctHelio.setTelescope(self.telescope) 
			correctHelio.setTarget(self.ephemeris.ra, self.ephemeris.dec)
			HJD = correctHelio.convertJDtoHJD(JD)
			for index, d in enumerate(self.data):
				d['HJD'] = HJD[index]
		else:
			print("No ephemeris loaded.")
			
	def computeJDFromMJD(self):
		for d in self.data:
			d['JD'] = d['MJD'] + 2400000.5
			print(d)	
		
		
	def computeHJDs(self):
		if self.hasEphemeris:
			print(self.id, self.ephemeris)
			MJD = self.getColumn('MJD')
			correctHelio = datetimelib.heliocentric()
			correctHelio.setTelescope(self.telescope) 
			correctHelio.setTarget(self.ephemeris.ra, self.ephemeris.dec)
			HJD = correctHelio.convertMJDtoHJD(MJD)
			self.setHJDs(MJD, HJD)

	def writeToJSON(self, filename):
		object = {}
		for key in self.__dict__.keys():
			data = getattr(self, key)
			if key == 'ephemeris': continue
			if type(data)==numpy.float32:
				data = float(data)
			if type(data)==numpy.ndarray:
				data = numpy.array(data).tolist()
			if type(data)==list:
				data = numpy.array(data).tolist()
			object[key] = data
			
		outputfile = open(filename, 'w')
		json.dump(object, outputfile, indent=4)
		outputfile.close()

	def fakeErrors(self):
		for d in self.data:
			d['fluxError'] = 1
		self.fluxErrorColumn = 'fluxError'

	def removeNaNs(self):
		newData = []
		print("Removing NaNs")
		print("\tOld length: %d"%(len(self.data)))
		for d in self.data:
			if numpy.isnan(d[self.fluxColumn]): continue
			newData.append(d)
		print("\tNew length: %d"%(len(newData)))
		self.data = newData

	def trimByErrors(self, error):
		newData = []
		for d in self.data:
			if d[self.fluxErrorColumn] < error:
				newData.append(d)
		self.data = newData

	def sigmaClip(self, size, nsigma):
		flux = self.getColumn(self.fluxColumn)
		fluxError = self.getColumn(self.fluxErrorColumn)
		print("In sigmaClip")
		newData = []
		for i in range(0, len(flux)-size):
			subarray = flux[i : i + size]
			print(subarray)
			mean = numpy.mean(subarray)
			std = numpy.std(subarray)
			nsigma = abs(flux[i]-mean)/std
			print("mean: %f, stddev %f, nsigma %f"%(mean, std, nsigma))
			if(nsigma>2):
				print("Rejecting point!")
				continue
			newData.append(self.data[i])
		for i in range(len(flux) - size, len(flux)):
			subarray = flux[len(flux)-size : ]
			print(i, subarray)
			mean = numpy.mean(subarray)
			std = numpy.std(subarray)
			nsigma = abs(flux[i]-mean)/std
			print("mean: %f, stddev %f, nsigma %f"%(mean, std, nsigma))
			if(nsigma>2):
				print("Rejecting point!")
				continue
			newData.append(self.data[i])
		print("old length %d : new length %d"%(len(self.data), len(newData)))
		self.data = newData
			
		
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
		
		

class loadPhotometry():
	def __init__(self, debug=False):
		self.name = ""
		self.debug = debug

	def loadFromHST(self, filename):
		columns = [ {'name': 'JD',      	  'type':'float'}, 
					{'name': 'counts',        'type':'float'}
		]

		data = []
		HSTFile = open(filename, 'rt')
		for line in HSTFile:
			line = line.strip()
			fields = line.split('\t')
			print(fields)
			d = {}
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				goodline = False
				try: 
					if column['type']=='str': value = str(value)
					if column['type']=='float': value = float(value)
					if column['type']=='int': value = int(value)
					d[column['name']] = value
				except ValueError:
					print("WARNING: Could not interpret the line [%d]: %s"%(lineNumber, line))
					break
				goodline = True
			if goodline: data.append(d)
					
		HSTFile.close()

		object = target("HST data")
		object.source = "HST"
		object.telescope = "HST"
		for d in data:
			object.appendData(d)
		print("%d observations loaded: "%len(object.data))
		object.fluxColumn = "counts"
		object.fluxErrorColumn = "none"
		object.dateColumn = "MJD"
		object.telescope = "HST"
		
		return object
		


	def loadFromW1m(self, filename):
		columns = [ {'name': 'filename',      'type':'string'}, 
					{'name': 'mid-time',      'type':'float'},
					{'name': 'star',          'type':'float'},
					{'name': 'noise',         'type':'float'},
					{'name': 'sky',           'type':'float'},
					{'name': 'x',             'type':'float'},
					{'name': 'y',             'type':'float'},
					{'name': 'FWHM',          'type':'float'} ]
		
		def readFromFile(filename, columnOffset):
			W1mFile = open(filename, 'rt')
			comments = ""
			data = []
			lineNumber = 0
			for line in W1mFile:
				lineNumber+= 1
				line = line.strip()
				if(line[0] == "#"):
					comments+= line + "\n"
					continue
				fields = line.split()
				# Strip out columns from other apertures
				firstFields = fields[0:2]
				nextFields = fields[columnOffset*6 + 2:]
				fields = firstFields + nextFields
				d = {}
				for index, column in enumerate(columns):
					value = fields[index].strip(' \t\n\r')
					goodline = False
					try: 
						if column['type']=='str': value = str(value)
						if column['type']=='float': value = float(value)
						if column['type']=='int': value = int(value)
						d[column['name']] = value
					except ValueError:
						print("WARNING: Could not interpret the line [%d]: %s"%(lineNumber, line))
						break
					goodline = True
				if goodline: data.append(d)
				
				if(self.debug): print(d)
			W1mFile.close()
			return comments, data


		objects = []

		comments, data = readFromFile(filename, 0)

		start = comments.find("FramePattern: ")
		end = start + comments[start:].find("\n")
		line = comments[start : end]
		targetString = generallib.getBetweenChars(line, '^', '(')
		object = target(targetString)
		object.telescope = "W1m"
		for d in data:
			object.appendData(d)
		print("%d observations loaded."%len(object.data))
		object.fluxColumn = "star"
		object.fluxErrorColumn = "noise"
		object.dateColumn = "mid-time"

		print(comments)
		start = comments.find("ReferenceTime: ") + len("ReferenceTime: ")
		end = start + comments[start:].find("\n")
		print(start, end)
		timeString = comments[start : end]
		print(timeString)
		from astropy.time import Time
		reference = Time(timeString, format='iso', scale='utc')
		for d in data:
			offset = datetime.timedelta(seconds = d['mid-time']) 
			d['JD'] = float((reference + offset).jd)
			# print("%f : %f : %f %f[%f]"%(d['mid-time'], offset.seconds, d['JD'], d['star'], d['noise']))
		object.dateColumn = "JD"
		objects.append(object)
		
		# Find the number of comparison stars by counting the columns in the data file.
		start = comments.find('### Filename')
		end = start + comments[start:].find("\n")
		headerString = comments[start:end]
		headers = headerString.split()
		count = 0
		for name in headers:
			if name=="Star": count+=1
		print("%d apertures in data file."%count)
		for comparison in range(1, count):
			print("Comparison number %d"%comparison)
			comments, data = readFromFile(filename, comparison)
			object = target("comparison %d"%comparison)
			object.telescope = "W1m"
			object.source = "W1m"
			for d in data:
				object.appendData(d)
			object.fluxColumn = "star"
			object.fluxErrorColumn = "noise"
			object.dateColumn = "mid-time"
			for d in data:
				offset = datetime.timedelta(seconds = d['mid-time']) 
				d['JD'] = float((reference + offset).jd)
			object.dateColumn = "JD"
			objects.append(object)

		return objects



	def loadFromASASSN(self, filename):
		columns = [ {'name': 'HJD',           'type':'float'}, 
					{'name': 'UTdate',        'type':'string'},
					{'name': 'camera',        'type':'string'},
					{'name': 'FWHM',          'type':'float'},
					{'name': 'limit',         'type':'float'},
					{'name': 'mag',           'type':'float'},
					{'name': 'mag_err',       'type':'float'},
					{'name': 'flux',          'type':'float'},
					{'name': 'flux_err',      'type':'float'},
					{'name': 'filter',        'type':'string'} ]
		
		asassnFile = open(filename, 'rt')
		comments = ""
		data = []
		lineNumber = 0
		for line in asassnFile:
			lineNumber+= 1
			line = line.strip()
			if(line[0] == "\\") or (line[0] =='|'):
				comments+= line + "\n"
				continue
			fields = line.split(',')
			d = {}
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				goodline = False
				try: 
					if column['type']=='str': value = str(value)
					if column['type']=='float': value = float(value)
					if column['type']=='int': value = int(value)
					d[column['name']] = value
				except ValueError:
					print("WARNING: Could not interpret the line [%d]: %s"%(lineNumber, line))
					break
				goodline = True
			if goodline: data.append(d)
			
			if(self.debug): print(d)
			
		object = target("ASASSN_data")
		object.source = "ASAS-SN"
		object.telescope = "ASASSN"
		for d in data:
			object.appendData(d)
		print("%d observations loaded: "%len(object.data))
		object.fluxColumn = "mag"
		object.fluxErrorColumn = "mag_err"
		object.dateColumn = "HJD"

		return object


	def loadFromAAVSO(self, filename):
		columns = [ {'name': 'rowID',         'type':'int'}, 
					{'name': 'target',        'type':'string'},
					{'name': 'JD',            'type':'float'},
					{'name': 'date',          'type':'string'},
					{'name': 'mag',           'type':'float'},
					{'name': 'mag_err',       'type':'float'},
					{'name': 'filter',        'type':'string'},
					{'name': 'observer',      'type':'string'} ]
		
		aavsoFile = open(filename, 'rt')

		self.debug = True
		comments = ""
		data = []
		for line in aavsoFile:
			line = line.strip()
			fields = line.split(',')
			if(line[0] == ","):
				comments+= line + "\n"
				continue
			print(fields)
			d = {}
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				if column['type']=='str': value = str(value)
				if column['type']=='float': value = float(value)
				if column['type']=='int': value = int(value)
				d[column['name']] = value
			if(self.debug): print(d)
			data.append(d)
		aavsoFile.close()

		object = target("AAVSO target")
		object.telescope = "AAVSO network"
		object.source = "AAVSO"
		for d in data:
			object.appendData(d)
		object.fluxColumn = "mag"
		object.fluxErrorColumn = "mag_err"
		object.dateColumn = "JD"

	
		return object
	

	def loadFromSuperWASP(self, filename):
		columns = [ {'name': 'TMID',          'type':'int'}, 
					{'name': 'flux2',         'type':'float'},
					{'name': 'flux2_err',     'type':'float'},
					{'name': 'tamflux2',      'type':'float'},
					{'name': 'tamflux2_err',  'type':'float'},
					{'name': 'imageid',       'type':'string'},
					{'name': 'ccdx',          'type':'int'},
					{'name': 'ccdy',          'type':'int'},
					{'name': 'flag',          'type':'int'},
					{'name': 'HJD',           'type':'float'},
					{'name': 'mag2',          'type':'float'},
					{'name': 'mag2_err',      'type':'float'},
					{'name': 'tammag2',       'type':'float'},
					{'name': 'tammag2_err',   'type':'float'} ]
		
		
		swaspFile = open(filename, 'rt')
		comments = ""
		data = []
		for line in swaspFile:
			line = line.strip()
			if(line[0] == "\\") or (line[0] =='|'):
				comments+= line + "\n"
				continue
			fields = line.split()
			d = {}
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				if column['type']=='str': value = str(value)
				if column['type']=='float': value = float(value)
				if column['type']=='int': value = int(value)
				d[column['name']] = value
			if(self.debug): print(d)
			data.append(d)

		object = target("SWASP_data")
		object.telescope = "SuperWASP"
		object.source = "SWASP"
		for d in data:
			object.appendData(d)
		object.fluxColumn = "tammag2"
		object.fluxErrorColumn = "tammag2_err"
		object.dateColumn = "HJD"

		return object

			
		

	def loadFromCRTS(self, filename):
		columns = [ {'name': 'ID',        'type':'str'}, 
					{'name': 'mag',       'type':'float'},
					{'name': 'err',       'type':'float'},
					{'name': 'ra',        'type':'float'},
					{'name': 'dec',       'type':'float'},
					{'name': 'MJD',       'type':'float'},
					{'name': 'blend',     'type':'int'} ]
		
		columnsLongForm = [ {'name': 'ID',        	'type':'str'}, 
							{'name': 'MasterFrame', 'type':'str'},
							{'name': 'CRTSID',      'type':'float'},
							{'name': 'mag',       	'type':'float'},
							{'name': 'err',      	'type':'float'},
							{'name': 'ra',        	'type':'float'},
							{'name': 'dec',       	'type':'float'},
							{'name': 'FWHM',       	'type':'float'},
							{'name': 'var',       	'type':'float'},
							{'name': 'FrameID',     'type':'float'},
							{'name': 'MJD',       	'type':'float'},
							{'name': 'airmass',     'type':'float'},
							{'name': 'exposure',    'type':'float'},
							{'name': 'X',       	'type':'float'},
							{'name': 'Y',       	'type':'float'},
							{'name': 'Flux',       	'type':'float'},
							{'name': 'Area',       	'type':'float'},
							{'name': 'Flags',       'type':'int'},
							{'name': 'Theta',       'type':'float'},
							{'name': 'Elong',       'type':'float'},
							{'name': 'NuMax',       'type':'float'},
							{'name': 'blend',     	'type':'int'} ]
					
		columnNames = ['ID', 'mag', 'err', 'ra', 'dec', 'MJD', 'blend' ]
		data = []
		catalinaFile = open(filename, 'rt')
		headings = catalinaFile.readline()
		if len(headings.split(',')) == 7:
			longForm = False
			print("Short form of Catalina data. Available columns are:", [c['name'] for c in columns])
		elif len(headings.split(',')) == 22:
			longForm = True
			columns = columnsLongForm
			print("Long form Catalina data. Available columns are:", [c['name'] for c in columns])
		else:
			print("Something is wrong with the input. CRTS data should have 8 cols (short form) or 22 cols (long form).")
			sys.exit()
		for line in catalinaFile:
			fields = line.split(',')
			d = {}
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				if column['type']=='str': value = str(value)
				if column['type']=='float': value = float(value)
				if column['type']=='int': value = int(value)
				d[column['name']] = value
			if(self.debug): print(d)
			data.append(d)
				

		
		# Separate different objects
		objects = []
		ids = []
		for d in data:
			id = d['ID']
			if id not in ids:
				ids.append(id)
		if self.debug: print(ids)
		
		for id in ids:
			o = target(id)
			o.source = "CRTS"
			for d in data:
				if d['ID'] == o.id:
					o.appendData(d)
			o.fluxColumn = "mag"
			o.fluxErrorColumn = "err"
			o.dateColumn = "MJD"
			objects.append(o)
		
		print("%d targets loaded"%len(objects))

		return objects