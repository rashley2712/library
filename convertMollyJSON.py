#!/usr/bin/env python3
import sys
import numpy, math
import argparse
import trm.molly
import spectrumlib, datetimelib

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Loads a series of spectra that were saved from Molly converts them to JSON format.')
	parser.add_argument('mollyfile', type=str, help='Molly file containing the spectra')
	parser.add_argument('-f', '--onefile', type=str, help="Save all the spectra into one file with this filename.")
	parser.add_argument('--suffix', type=str, help='Suffix to add to the end of the filenames.')
	parser.add_argument('-e', type=str, help='[Optional] Ephemeris data file')
	arg = parser.parse_args()
	print(arg)

	if arg.onefile is None:
		onefile = False
	else: onefile = True

	if arg.e!=None:
		# Load the ephemeris file
		hasEphemeris = True
		ephemeris = datetimelib.ephemerisObject()
		ephemeris.loadFromFile(arg.e)
		print(ephemeris)
	else:
		hasEphemeris = False


	spectra = [] 
	
	mollyFilename = arg.mollyfile
	mollyFile = trm.molly.rmolly(mollyFilename)
		
	for index, r in enumerate(mollyFile):
		wavelengths = []
		flux = []
		fluxErrors = []	
		for f, fe, w in zip(r.f, r.fe, r.wave):
			# print(w, f, fe)
			wavelengths.append(w)
			flux.append(f)
			fluxErrors.append(fe)
 
		head = r.head
		
		spectrum = spectrumlib.spectrumObject()
		npoints = spectrum.setData(wavelengths, flux, fluxErrors)
		targetName = spectrum.parseHeaderInfo(head)
		spectrum.wavelengthUnits = "\\A"
		spectrum.fluxLabel = r.label
		spectrum.fluxUnits = r.units
		print("Flux units:", spectrum.fluxUnits)
		# spectrum.fluxUnits = "relative counts"
		
		print("Parsed headers of %s for HJD: %f"%(targetName, spectrum.HJD))
		spectrum.name = "%s-%f"%(spectrum.objectName, spectrum.HJD)
		spectra.append(spectrum)
		
	numSpectra = len(spectra)

	print("%d spectra loaded."%numSpectra)

	count = 0
	spectrumDB  = spectrumlib.spectrumArray()
	for s in spectra:
		outname = "%s_%f.json"%(s.objectName, s.HJD)
		if hasEphemeris:
			if arg.suffix!=None:
				outname = "%s_%f_%s.json"%(s.objectName, ephemeris.getPhase(s.HJD), arg.suffix)
			else: 
				outname = "%s_%f.json"%(s.objectName, ephemeris.getPhase(s.HJD))
			
		print("Writing to %s"%outname)
		if not onefile: s.writeToJSON(outname, clobber=False)
		spectrumDB.addSpectrum(s)
		count+=1

	if onefile:
		spectrumDB.filename = arg.onefile
		spectrumDB.writeDB()
		print("Written %d spectra to %s"%(count, arg.onefile))
	else:
		print("Written %d spectra to disk."%count)

