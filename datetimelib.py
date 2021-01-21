import numpy
import json
import generallib
import astropy
from astropy.coordinates import SkyCoord, EarthLocation
import sys

class heliocentric:
	telescopes = [
		{ 'name': 'CSS', 'longitude': -110.73167 , 'latitude': 32.41667, 'elev': 2510. },
		{ 'name': 'SSS', 'longitude': -149.1 , 'latitude': -31.3, 'elev': 1150. },
		{ 'name': 'MLS', 'longitude': -110.7889 , 'latitude': 32.4433, 'elev': 2790. },
		{ 'name': 'SuperWASP', 'longitude': -17.8816 , 'latitude': 28.7606, 'elev': 2326. },
		{ 'name': 'W1m', 'longitude': -17.8816 , 'latitude': 28.7606, 'elev': 2326. },
		{ 'name': 'WHT', 'longitude': -17.8816 , 'latitude': 28.7606, 'elev': 2326. },
		{ 'name': 'HST', 'longitude': -17.8816 , 'latitude': 28.7606, 'elev': 2326. },
		{ 'name': 'Arto', 'longitude': 26.6 , 'latitude': 62.15, 'elev': 150. },		
		{ 'name': 'TNT', 'longitude': 98.48 , 'latitude': 18.57, 'elev': 2457. },		
		{ 'name': 'NTT', 'longitude': 289.27 , 'latitude': -29.2567, 'elev': 2347. }		
		]
	
	def __init__(self):
		self.telescope = None
		self.target = None
		self.debug = False

	def showTelescopes(self):
		for t in heliocentric.telescopes:
			print("%s (%.2f, %.2f, %.0f)"%(t['name'], t['longitude'], t['latitude'], t['elev']))

	def setTelescope(self, telescopeName):
		for t in heliocentric.telescopes:
			if t['name'] == telescopeName:
				self.telescope = t
				return True
		return False

	def setTarget(self, ra, dec):
		self.target = { 'ra': ra, 'dec': dec}

	def convertMJDtoHJD(self, MJD):
		if self.telescope is None:
			print("We don't know the location of the telescope. Exiting")
			return

		obsLocation = astropy.coordinates.EarthLocation(lon = self.telescope['longitude'], lat = self.telescope['latitude'], height = self.telescope['elev'])

		if self.target is None:
			print("We don't know the coordinates of the target. Exiting.")

		targetRADEC = generallib.toSexagesimal((self.target['ra'], self.target['dec']))
		if self.debug: print("Target position: %s (%f, %f)"%(targetRADEC, self.target['ra'], self.target['dec']))
		if self.debug: print('Observatory: ', self.telescope)

		targetCoords = astropy.coordinates.SkyCoord(self.target['ra'], self.target['dec'], unit='deg')

		from astropy import time, coordinates as coord, units as u
		times = time.Time(MJD, format='mjd', scale='utc', location=obsLocation)
		ltt_bary = times.light_travel_time(targetCoords)
		# print(ltt_bary) 
		time_barycentre = times.tdb + ltt_bary
		#time_barycentre = times.tdb 
		corrected_times = time_barycentre + 2400000.5
		hjd = [t.mjd for t in corrected_times]
		return hjd



	def convertJDtoHJD(self, JD):
		if self.telescope is None:
			print("We don't know the location of the telescope. Exiting")
			return

		obsLocation = astropy.coordinates.EarthLocation(lon = self.telescope['longitude'], lat = self.telescope['latitude'], height = self.telescope['elev'])

		if self.target is None:
			print("We don't know the coordinates of the target. Exiting.")

		targetRADEC = generallib.toSexagesimal((self.target['ra'], self.target['dec']))
		if self.debug: 
			print("Target position: %s (%f, %f)"%(targetRADEC, self.target['ra'], self.target['dec']))
			print('Observatory: ', self.telescope)

		targetCoords = astropy.coordinates.SkyCoord(self.target['ra'], self.target['dec'], unit='deg')

		from astropy import time, coordinates as coord, units as u
		times = time.Time(JD, format='jd', scale='utc', location=obsLocation)
		ltt_bary = times.light_travel_time(targetCoords)
		# print(ltt_bary) 
		time_barycentre = times.tdb + ltt_bary
		#time_barycentre = times.tdb 
		corrected_times = time_barycentre
		hjd = [t.jd for t in corrected_times]
		return hjd

	def convertHJDtoJD(self, HJD):
		if self.telescope is None:
			print("We don't know the location of the telescope. Exiting")
			return

		obsLocation = astropy.coordinates.EarthLocation(lon = self.telescope['longitude'], lat = self.telescope['latitude'], height = self.telescope['elev'])

		if self.target is None:
			print("We don't know the coordinates of the target. Exiting.")

		targetRADEC = generallib.toSexagesimal((self.target['ra'], self.target['dec']))
		if self.debug: 
			print("Target position: %s (%f, %f)"%(targetRADEC, self.target['ra'], self.target['dec']))
			print('Observatory: ', self.telescope)

		targetCoords = astropy.coordinates.SkyCoord(self.target['ra'], self.target['dec'], unit='deg')

		from astropy import time, coordinates as coord, units as u
		times = time.Time(HJD, format='jd', scale='utc', location=obsLocation)
		ltt_bary = times.light_travel_time(targetCoords)
		# print(ltt_bary) 
		time_barycentre = times.tdb - ltt_bary
		#time_barycentre = times.tdb 
		corrected_times = time_barycentre
		jd = corrected_times.jd
		return jd


class ephemerisObject:
	def __init__(self):
		self.T0 = 0
		self.T0_error = 0
		self.Period = 0
		self.Period_error = 0
		self.ra = 0
		self.dec = 0
		self.gamma = None
		self.gamma_error = None
		self.K2 = None
		self.K2_error = None
		self.orbit = False
		self.debug = False

	def setCoords(self, ra, dec):
		self.ra = ra
		self.dec = dec


	def getPhase(self, HJD):
		HJD_difference = HJD - self.T0
		#norbits = int( HJD_difference / self.Period)
		phase = (HJD_difference % self.Period) / self.Period
		#if phase>0.5: phase = phase - 1.0
		return phase

	def getRV(self, HJD):
		if self.K2 is None: 
			print("No RV solution loaded, so cannot determine RV")			
			return
		phase = self.getPhase(HJD)
		RV = self.gamma + self.K2 * numpy.sin(numpy.pi*2.0*phase)
		return RV

	def getNearestCyclenumber(self, HJD):
		HJD_difference = HJD - self.T0
		norbits = round(HJD_difference / self.Period)
		return norbits		

	def getOrbits(self, HJD):
		HJD_difference = HJD - self.T0
		norbits = int(HJD_difference / self.Period)
		upperorbit = int(round(HJD_difference / self.Period))
		return norbits, upperorbit

	def getOrbitsSince(self, start, now):
		HJD_difference = now - start
		norbits = HJD_difference / self.Period
		return norbits

	def getOffsetOrbits(self, HJD):
		HJD_difference = HJD - self.T0 - self.Period/2.0
		norbits = int( HJD_difference / self.Period)
		return norbits

	def loadFromFile(self, filename):
		file = open(filename, 'r')
		for line in file:
			if line[0] == '#': continue
			tokens = line.split()
			if len(tokens)>1:
				if tokens[0] == 'T0': self.T0 = float(tokens[1])
				if tokens[0] == 'T0_error': self.T0_error = float(tokens[1])
				if tokens[0] == 'E': self.Period = float(tokens[1])
				if tokens[0] == 'E_error': self.Period_error = float(tokens[1])
				if tokens[0] == 'J2000':
					coords = tokens[1:]
					self.ra, self.dec = self.parseCoords(coords)
				if tokens[0] == "Gamma": self.gamma = float(tokens[1])
				if tokens[0] == "Gamma_error": self.gamma_error = float(tokens[1])
				if tokens[0] == "K2": self.K2 = float(tokens[1])
				if tokens[0] == "K2_error": self.K2_error = float(tokens[1])
				if tokens[0] == "N": self.N = int(tokens[1])
				if tokens[0] == "Chi2_red":
					self.chi2red = float(tokens[1])
		if self.Period != 0: self.orbit = True
		file.close()
				

	def parseCoords(self, coords):
		if self.debug: print("Given coords:", coords)
		raHours = int(coords[0])
		raMinutes = int(coords[1])
		raSeconds = float(coords[2])
		raFraction = float(raSeconds)/3600. + float(raMinutes)/60.
		raTotal = float(raHours) + raFraction
		raDegrees = raTotal * 15.0
		decDegreesString = str(coords[3])
		sign = decDegreesString[0]
		decDegrees = abs(float(coords[3]))
		decMinutes = float(coords[4])
		decSeconds = float(coords[5])
		decCalc = decDegrees + decMinutes/60. + decSeconds/3600.
		if sign == '-':
			decCalc = -1.0 * decCalc
		if self.debug: print("RA:", raDegrees, "DEC:", decCalc)
		return raDegrees, decCalc

	def __str__(self):
		outString = "T0: %7.8f [%7.8f] + E x %7.10f [%7.10f]"%(self.T0, self.T0_error, self.Period, self.Period_error)
		return outString
