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

	def setColumn(self, name, data):
		for index, d in enumerate(self.data):
			d[name] = data[index]
		
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
			if self.debug: print("phase range: %f to %f"%(phaseStart, phaseEnd))
			binData = []
			for d in self.data:
				if (d['phase']>=phaseStart) and (d['phase']<phaseEnd):
					binData.append(d)
			if self.debug: print("%d points in this phase bucket"%len(binData))
			binNumbers = [d[self.fluxColumn] for d in binData]
			if self.debug: print(binNumbers)
			if len(binData)==0:
				if self.debug: print("WARNING: No data in phase range: %f to %f"%(phaseStart, phaseEnd))
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
					weight = 1/(d[self.fluxErrorColumn]**2)
					sum+= d[self.fluxColumn] * weight
					errorsum+= weight
				weighted_mean = sum/errorsum
				mean = sum/errorsum
				weighted_error = numpy.sqrt(1/errorsum)
				#squ_sum = 0
				#for d in binData:
				#	squ_sum+= (d[self.fluxColumn]-mean)**2 / d[self.fluxErrorColumn]
				#stddev = numpy.sqrt( squ_sum/ N / errorsum)
				#error = numpy.sqrt( squ_sum / N / (N-1) / errorsum)
				#if self.debug: print("Weighted mean: %f [%f]"%(mean, error))
				if self.debug: print("Weighted mean: %f [%f]"%(weighted_mean, weighted_error))
				mean = weighted_mean
				error = weighted_error
				#straightMean = numpy.mean(binNumbers)
				#straightStdDev = numpy.std(binNumbers)
				#sterror = straightStdDev / numpy.sqrt(len(binNumbers))
				#if self.debug: print("Straight mean: %f [%f] (%f)"%(straightMean, straightStdDev, sterror))	

			bins.append(mean)
			phases.append(midPhase)
			errors.append(error)
			#stddevs.append(stddev)
		#if self.debug: print(phases, bins)
		self.bins = bins
		self.phases = phases
		self.errors = errors
		#self.stddevs = stddevs
		

	def setHJDs(self, MJD, HJD):
		keys = [d['MJD'] for d in self.data]
		dates = zip(MJD, HJD)
		for index, d in enumerate(dates):
			self.data[index]['HJD'] = d[1]

	def computeHJDfromJD(self):
		if self.hasEphemeris:
			#print(self.id, self.ephemeris)
			JD = self.getColumn('JD')
			correctHelio = datetimelib.heliocentric()
			correctHelio.setTelescope(self.telescope) 
			correctHelio.setTarget(self.ephemeris.ra, self.ephemeris.dec)
			HJD = correctHelio.convertJDtoHJD(JD)
			for index, d in enumerate(self.data):
				d['HJD'] = HJD[index]
		else:
			print("No ephemeris loaded.")
		self.dateColumn = 'HJD'
			
	def computeJDfromMJD(self):
		for d in self.data:
			d['JD'] = d['MJD'] + 2400000.5
		self.dateColumn = 'JD'		
		
		
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

	def duplicateColumn(self, old, new):
		for d in self.data:
			d[new] = d[old]

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

	def trimDates(self, startDate, endDate, dateColumn):
		newData = []
		print("Trimming dates in range (%f, %f)."%(startDate, endDate))
		print("\tOld length: %d"%(len(self.data)))
		for d in self.data:
			if d[dateColumn]<startDate or d[dateColumn]>endDate:
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
		self.debug=False
		if self.debug: print("In sigmaClip")
		newData = []
		for i in range(0, len(flux)-size):
			subarray = flux[i : i + size]
			if self.debug: print(subarray)
			mean = numpy.mean(subarray)
			std = numpy.std(subarray)
			nsigma = abs(flux[i]-mean)/std
			if self.debug: print("mean: %f, stddev %f, nsigma %f"%(mean, std, nsigma))
			if(nsigma>2):
				if self.debug: print("Rejecting point!")
				continue
			newData.append(self.data[i])
		for i in range(len(flux) - size, len(flux)):
			subarray = flux[len(flux)-size : ]
			if self.debug: print(i, subarray)
			mean = numpy.mean(subarray)
			std = numpy.std(subarray)
			nsigma = abs(flux[i]-mean)/std
			if self.debug: print("mean: %f, stddev %f, nsigma %f"%(mean, std, nsigma))
			if(nsigma>2):
				if self.debug: print("Rejecting point!")
				continue
			newData.append(self.data[i])
		if self.debug: print("old length %d : new length %d"%(len(self.data), len(newData)))
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


	def loadFromASASSNJSON(self, filename):
		inputfile = open(filename, "r")
		jsonObjectList = json.load(inputfile)
		inputfile.close()
		
		data = []
		for jo in jsonObjectList:
			if jo['mag_err']==99.99: 
				jo['limit'] = 'upper'
				jo['mag_err'] = 0
			else:
				jo['limit'] = 'value'
			print(jo)
			data.append(jo)
		
		
		object = target("ASASSN_data")
		object.source = "ASAS-SN"
		object.telescope = "ASASSN"
		for d in data:
			object.appendData(d)
		print("%d observations loaded: "%len(object.data))
		object.fluxColumn = "mag"
		object.fluxErrorColumn = "mag_err"
		object.dateColumn = "hjd"

		return object

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




	def loadFromAAVSOArto(self, filename):
		columns = [ {'name': 'JD',          'type':'float'}, 
					{'name': 'mag',         'type':'float'},
					{'name': 'fainterthan', 'type':'float'},
					{'name': 'mag_err',     'type':'float'},
					{'name': 'uncertaintyhq', 'type':'string'},
					{'name': 'band', 		'type':'int'} ]

		aavsoFile = open(filename, 'rt')

		self.debug = True
		comments = ""
		data = []
		rejectedLines = []
		for index, line in enumerate(aavsoFile):
			if index==0:
				headings = line.strip()
				continue
			fields = line.strip().split('\t')
			print(fields)
			d = {}
			reject = False
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				try:
					if column['type']=='str': value = str(value)
					if column['type']=='float': value = float(value)
					if column['type']=='int': value = int(value)
					d[column['name']] = value
				except ValueError as e:
					#print(e)
					#print("line number: ", index)
					rejectedLines.append(line.strip())
					reject = True
			if(self.debug): print(d)
			if not reject: data.append(d)
				
		print("%d rejected lines"%len(rejectedLines))
		aavsoFile.close()
		object = target("AAVSO target")
		object.telescope = "Arto"
		object.source = "AAVSO"
		for d in data:
			object.appendData(d)
		object.fluxColumn = "mag"
		object.fluxErrorColumn = "mag_err"
		object.dateColumn = "JD"
		print("%d data points loaded."%len(data))
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

			
	def loadFromCRTSDR3(self, filename):
		columns = [ {'name': 'ID',        'type':'str'}, 
					{'name': 'ra',        'type':'float'},
					{'name': 'dec',       'type':'float'},
					{'name': 'mag',       'type':'float'},
					{'name': 'err',       'type':'float'},
					{'name': 'MJD',       'type':'float'},
					{'name': 'Sep',       'type':'float'} ]
				
		CRTSFile = open(filename, 'rt')
		comments = ""
		data = []
		for index, line in enumerate(CRTSFile):
			line = line.strip()
			if(index==0):
				comments+=line + "\n"
				continue
			if(line[0] == "#"):
				comments+= line + "\n"
				continue
			fields = line.split(',')
			if self.debug: print(fields)
			d = {}
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				if column['type']=='str': value = str(value)
				if column['type']=='float': value = float(value)
				if column['type']=='int': value = int(value)
				d[column['name']] = value
			if(self.debug): print(d)
			data.append(d)

		object = target(data[0]['ID'])	
		object.telescope = "CSS"
		object.source = "CRTS"
		for d in data:
			object.appendData(d)
		object.fluxColumn = "mag"
		object.fluxErrorColumn = "mag_err"
		object.dateColumn = "MJD"

		CRTSFile.close()
		return [object]


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

	def loadFromArto(self, filename):
		columns = [ {'name': 'Object',        'type':'string'}, 
					{'name': 'JD',            'type':'float'},
					{'name': 'mag',           'type':'float'},
					{'name': 'mag_err',     'type':'float'} ]
		
		artoFile = open(filename, 'rt')
		comments = ""
		data = []
		for line in artoFile:
			line = line.strip()
			if(line[0] == "#"):
				comments+= line + "\n"
				continue
			fields = line.split(',')
			if self.debug: print(fields)
			d = {}
			limit = 'value'
			for index, column in enumerate(columns):
				value = fields[index].strip(' \t\n\r')
				if value[0]=='<': 
					value = value[1:]
					limit = 'upper'
				if column['type']=='str': value = str(value)
				if column['type']=='float': value = float(value)
				if column['type']=='int': value = int(value)
				d[column['name']] = value
			d['limit'] = limit
			if(self.debug): print(d)
			data.append(d)

		object = target(data[0]['Object'])	
		object.telescope = "Arto"
		object.source = "Arto"
		for d in data:
			object.appendData(d)
		object.fluxColumn = "mag"
		object.fluxErrorColumn = "mag_err"
		object.dateColumn = "JD"

		artoFile.close()
		return object

	def loadFromUCAM(self, filename):
		""" Loads the photometry data from a log file created by Tom Marsh's ULTRACAM pipeline. """
		ultracam = False
		ultraspec = True
		inputFile = open(filename, 'r')
		objects = []
		xValues = []
		yValues = []
		frameList = []
		headerBlock = ""
		runName = "--unknown--"
		telescope = "--unknown--"
		targetName = "--unknown--"
		filterName = "--unknown--"
		PI = "--unknown--"
		columnCount = 0
		uniqueCCDs = []
		for line in inputFile:
			if line[0] == '#':
				headerBlock+=line
				if ("target" in line) and ("estimated" not in line):
					targetName = generallib.getBetweenChars(line, '=', '/').strip()
					if self.debug: print("Target: %s"%targetName)
				if ("filters" in line):
					filterName = generallib.getBetweenChars(line, '=', '/').strip()
					if self.debug: print("Filters: %s"%filterName)
				if ("Telescope" in line) and ("observing" not in line):
					telescopeName = generallib.getBetweenChars(line, '=', '/').strip()
					if self.debug: print("Telescope name: %s"%telescopeName)
				if (" pi " in line):
					PI = generallib.getBetweenChars(line, '=', '/').strip()
					if self.debug: print("PI: %s"%PI)
				if (" Data file name " in line):
					runName = generallib.getBetweenChars(line, '=', '\n').strip()
					if self.debug: print("run data file: %s"%runName)
				if (" Server file name " in line):
					runName = generallib.getBetweenChars(line, '=', '\n').strip()
					if self.debug: print("run data file: %s"%runName)
					
			if line[0] != '#':
				params = line.split()
				# print params
				frameIndex = int(params[0])
				CCD = int(params[4])
				if CCD not in uniqueCCDs: uniqueCCDs.append(CCD)
				frameList.append(frameIndex)
				columnCount = len(params)
		firstFrame = frameList[0]
		
		numApertures = int( ((columnCount-7)/14) )
		if self.debug: print("ColumnCount: ", columnCount, "which means %d apertures."%numApertures)
		# frameList = generalUtils.removeDuplicatesFromList(frameList)
		if self.debug: print("The run in file %s contains %d frames. Start frame: %d End frame: %d"%(filename, len(frameList), min(frameList), max(frameList)))
		if len(uniqueCCDs) == 3:
			if self.debug: print("This file has 3 CCDs. It is an ULTRACAM file.")
			ultracam = True
			ultraspec = False
		if len(uniqueCCDs) == 1: 
			if self.debug: print("This file has 1 CCD. It is an ULTRASPEC file.")
			ultracam = False
			ultraspec = True

		if (ultracam): CCDs = [1, 2, 3]
		else: CCDs = [1]
		for CCD in CCDs: 
			for aperture in range(1, numApertures+1):
				apertureIndex = 14*(aperture-1) + 7
				if self.debug: print("Reading data for aperture %d, CCD %d"%(aperture, CCD))
				inputFile.seek(0)
				MJDs = []
				counts = []
				skys = []
				sigmas = []
				errors = []
				timeFlags = []
				exposures = []
				FWHMs = []
				betas = []
				xs = []
				ys = []
				lineCounter = 0
				data = []
				for line in inputFile:
					lineCounter+= 1
					if self.debug: 
						sys.stdout.write("\rLine number: %d    "%(lineCounter))
						sys.stdout.flush()
					if line[0] != '#':
						params = line.split()
						# print params
						CCDValue = int(params[4])
						apertureValue = int(params[apertureIndex])
						logEntry = {}
						if CCDValue == CCD: 
							frameIndex = int(params[0])
							MJD = float(params[1])
							exposure = float(params[3])
							FWHM = float(params[5])
							beta = float(params[6])
							x = float(params[apertureIndex + 1])
							y = float(params[apertureIndex + 2])
							counts = float(params[apertureIndex + 7])
							counts_err = float(params[apertureIndex + 8])
							sky = float(params[apertureIndex + 9])
							errorflag = int(params[apertureIndex + 13])
							logEntry = {
								'MJD': MJD,
								'exposure': exposure, 
								'FWHM': FWHM,
								'beta': beta,
								'x': x,
								'y': y,
								'counts' : counts,
								'counts_err' : counts_err,
								'sky' : sky,
								'errorflag' : errorflag
							}
							if errorflag<8: data.append(logEntry)
						
			
				object = target(targetName + "_" + str(apertureValue))
				object.telescope = telescopeName
				if telescopeName=="Thai National Observatory 2.4m": object.telescope = "TNT"
				if ultracam: object.source = "ULTRACAM"
				else: object.source = "ULTRASPEC"
				for d in data:
					object.appendData(d)
				object.fluxColumn = "counts"
				object.fluxErrorColumn = "counts_err"
				object.dateColumn = "MJD"
				objects.append(object)
				if self.debug: print()
		inputFile.close()
		return objects


	def loadFromHCAM(self, filename):
		""" Loads the photometry data from a log file created by Tom Marsh's HIPERCAM pipeline. """
		ultracam = False
		ultraspec = True
		hipercam = False
		inputFile = open(filename, 'r')
		objects = []
		xValues = []
		yValues = []
		frameList = []
		headerBlock = ""
		runName = "--unknown--"
		telescopeName = "--unknown--"
		targetName = "--unknown--"
		filterName = "--unknown--"
		PI = "--unknown--"
		columnCount = 0
		uniqueCCDs = []
		for line in inputFile:
			if line[0] == '#':
				headerBlock+=line
				if ("target" in line) and ("estimated" not in line):
					targetName = generallib.getBetweenChars(line, '=', '/').strip()
					if self.debug: print("Target: %s"%targetName)
				if ("log" in line):
					targetName = line.split('=')[-1].strip()
					if self.debug: print("Log: %s"%filterName)
				if ("filters" in line):
					filterName = generallib.getBetweenChars(line, '=', '/').strip()
					if self.debug: print("Filters: %s"%filterName)
				if ("Telescope" in line) and ("observing" not in line):
					telescopeName = generallib.getBetweenChars(line, '=', '/').strip()
					if self.debug: print("Telescope name: %s"%telescopeName)
				if (" pi " in line):
					PI = generallib.getBetweenChars(line, '=', '/').strip()
					if self.debug: print("PI: %s"%PI)
				if (" Data file name " in line):
					runName = generallib.getBetweenChars(line, '=', '\n').strip()
					if self.debug: print("run data file: %s"%runName)
				if (" Server file name " in line):
					runName = generallib.getBetweenChars(line, '=', '\n').strip()
					if self.debug: print("run data file: %s"%runName)
					
			if line[0] != '#':
				params = line.split()
				frameIndex = int(params[1])
				CCD = int(params[0])
				if CCD not in uniqueCCDs: uniqueCCDs.append(CCD)
				frameList.append(frameIndex)
				columnCount = len(params)
		
		numApertures = int( ((columnCount-7)/15) )
		if self.debug: print("ColumnCount: ", columnCount, "which means %d apertures."%numApertures)
		frameList = generallib.removeDuplicatesFromList(frameList)
		if self.debug: 
			print("The run in file %s contains %d frames. Start frame: %d End frame: %d"%(filename, len(frameList), min(frameList), max(frameList)))
			print("CCDS:", uniqueCCDs)
		ultracam = False
		hipercam = False
		ultraspec = False
		if len(uniqueCCDs) == 3:
			if self.debug: print("This file has 3 CCDs. It is an ULTRACAM file.")
			ultracam = True
		if len(uniqueCCDs) > 3:
			if self.debug: print("This file has mare than 3 CCDs. It is a HIPERCAM file.")
			hipercam = True
		if len(uniqueCCDs) == 1: 
			if self.debug: print("This file has 1 CCD. It is an ULTRASPEC file.")
			ultraspec = True

		if (ultracam): CCDs = [1, 2, 3]
		if (ultraspec): CCDs = [1]
		if (hipercam): CCDs = uniqueCCDs
		for CCD in CCDs: 
			for aperture in range(1, numApertures+1):
				apertureIndex = 15*(aperture-1) + 7
				if self.debug: print("Reading data for aperture %d, CCD %d"%(aperture, CCD))
				inputFile.seek(0)
				MJDs = []
				counts = []
				skys = []
				sigmas = []
				errors = []
				timeFlags = []
				exposures = []
				FWHMs = []
				betas = []
				xs = []
				ys = []
				lineCounter = 0
				data = []
				for line in inputFile:
					lineCounter+= 1
					if self.debug: 
						sys.stdout.write("\rLine number: %d    "%(lineCounter))
						sys.stdout.flush()
					if line[0] != '#':
						params = line.split()
						# print params
						CCDValue = int(params[0])
						#apertureValue = int(params[apertureIndex])
						logEntry = {}
						if CCDValue == CCD: 
							frameIndex = int(params[1])
							MJD = float(params[2])
							exposure = float(params[4])
							FWHM = float(params[5])
							beta = float(params[6])
							x = float(params[apertureIndex])
							y = float(params[apertureIndex + 2])
							counts = float(params[apertureIndex + 8])
							counts_err = float(params[apertureIndex + 9])
							sky = float(params[apertureIndex + 10])
							errorflag = int(params[apertureIndex + 14])
							logEntry = {
								'MJD': MJD,
								'exposure': exposure, 
								'FWHM': FWHM,
								'beta': beta,
								'x': x,
								'y': y,
								'counts' : counts,
								'counts_err' : counts_err,
								'sky' : sky,
								'errorflag' : errorflag
							}
							if errorflag==0: data.append(logEntry)
							if self.debug: print(logEntry) 
						
			
				object = target("%s_aper_%d_ccd_%d"%(targetName, aperture, CCD))
				object.telescope = "NTT"				
				object.telescopeName = "New Technology Telescope"
				if ultracam: object.source = "ULTRACAM"
				else: object.source = "ULTRASPEC"
				for d in data:
					object.appendData(d)
				object.fluxColumn = "counts"
				object.fluxErrorColumn = "counts_err"
				object.dateColumn = "MJD"
				objects.append(object)
				if self.debug: print()
		inputFile.close()
		return objects