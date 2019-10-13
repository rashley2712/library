
class configClass():
	def __init__(self, debug = False):
		self.debug = debug
		self.comments = ""
		self.settings = []
		self.defaultSettings()
		
	def show(self):
		for s in self.settings:
			print(s)
		
	def addSetting(self, name, type):
		self.settings[name] = "test"

	def loadSetting(self, setting, original):
		fields = original.split()
		if self.debug: print("loading setting: ", setting['name'])
		try:
			if setting['type'] == "string":
				setattr(self, setting['name'], original[len(setting['name'])+1:])
			if setting['type'] == "tuple":
				setattr(self, setting['name'], (float(fields[1]), float(fields[2])))
			if setting['type'] == "integer":
				setattr(self, setting['name'], int(fields[1]))
		except ValueError:
			print("Warning: Could not read the setting for: %s"%(setting['name']))

		return setting

	def read(self, name):
		for s in self.settings:
			if s['name'] == name: return s

	def load(self, filename="plot.cfg"):
		try:
			fileHandle = open(filename, 'rt')
		except FileNotFoundError:
			print("Could not open config file %s"%filename)
			return
		for line in fileHandle:
			line = line.strip()
			if line[0]=='#': 
				self.comments+= line[1:] + "\n"
				continue
			moreComments = line.split('#')
			if len(moreComments)>1: 
				self.comments+= moreComments[1] + "\n"
			 
			original = moreComments[0]
			line = line.strip().split()
			if self.debug: print(line)
			command = line[0]
			for setting in self.settings:
				if (setting['name']==command):
					setting = self.loadSetting(setting, original)
		
		fileHandle.close()

		if self.debug: print("Comments", self.comments)


	def defaultSettings(self):
		setting = { 'name' : 'wavelengthRange', 'type' : 'tuple', 'description': 'The upper and lower wavelengths of the plot.'}
		self.settings.append(setting)
		setting = { 'name' : 'plotDimensions', 'type' : 'tuple', 'description': 'Width and height of the plot (inches).'}
		self.settings.append(setting)
		setting = { 'name' : 'number', 'type' : 'integer', 'description': 'Plot number in molly file.'}
		self.settings.append(setting)
		setting = { 'name' : 'xlabel', 'type' : 'string', 'description': 'X-axis label.', 'default': 'Wavelength (\AA)'}
		self.settings.append(setting)
		setting = { 'name' : 'ylabel', 'type' : 'string', 'description': 'Y-axis label.', 'default':'Flux'}
		self.settings.append(setting)
		setting = {
			'name': "title", 
			'type': "string",
			'description': 'Title of the plot.',
		}
		self.settings.append(setting)
		
		for setting in self.settings:
			try:
				setattr(self, setting['name'], setting['default'])
			except KeyError:
				continue

		