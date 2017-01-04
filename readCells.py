"""Helper set of functions to read in and parse multiple types of input sources."""

import sys
import os
import datetime
from bs4 import BeautifulSoup
from shapely.geometry.polygon import Polygon

def read(ftype, inDir, inSuffix, startTime, endTime):
	"""
	Determines the user-specified file type and parses it accordingly
	
	Parameters
	----------
	ftype : string
		The file type to process: segmotion, probsevere, or ryan
	inDir : string
		The input directory
	inSuffix : string
		The lowest subdirectory of the input directory
	startTime : string
		The earliest time to process
	endTime : string
		The latest time to process
		
	Returns
	-------
	list
		[stormCells, totNumCells, numTrackTimes, dates] - List containing dicts of all 
		storm cells, the number of cells, the number of files, and all valid dates
		
	"""

	if ftype == 'ryan': return readRyan(inDir, inSuffix, startTime, endTime)	
	elif ftype == 'segmotion': return readSegmotion(inDir, inSuffix, startTime, endTime)
	elif ftype == 'probsevere': return readProbSevere(inDir, inSuffix, startTime, endTime)


def readRyan(inDir, inSuffix, startTime, endTime):
	"""
	Parses post-processed segmotion files (.data) from Ryan's original code
	
	Parameters
	----------
	inDir : string
		The input directory
	inSuffix : string
		The lowest subdirectory of the input directory
	startTime : string
		The earliest time to process
	endTime : string
		The latest time to process
		
	Returns
	-------
	list
		[stormCells, totNumCells, numTrackTimes, dates] - List containing dicts of all 
		storm cells, the number of cells, the number of files, and all valid dates
		
	"""
	
	numTrackTimes = 0
	totNumCells = 0
	stormCells = {} 
	dates = []
	
	# Read in Ryan files
	for root, dirs, files in os.walk(inDir):
		if inSuffix != '' and not (files and not dirs and os.path.split(root)[-1] == inSuffix): continue
		for trackFile in files:
			if trackFile.endswith('.data'):
				
				# Skip hidden files
				if trackFile.startswith('._'): continue
				
				# Check if file falls in date range
				try:
					fileDate = datetime.datetime.strptime(str(trackFile).split('_')[0], '%Y-%m-%d-%H%M%S')
				except ValueError:
					print 'File ' + str(trackFile) + ' has an invalid name.  Expected format YYYY-MM-DD-hhmmss_...'
					continue
				if not startTime <= fileDate < endTime:
					continue
				if fileDate.date() not in dates: dates.append(fileDate.date())
				
				# Open file
				f = open(root + '/' + trackFile)
				lines = f.readlines()
				f.close()
				
				# Skip probSevere files
				if int(lines[28].split()[0]) == 1:
					print '\nWARNING: Unable to process storm objects from probSevere in Ryan format.  Use "-t probsevere" instead.'
					print str(trackFile) + ' will be skipped.\n'
					continue
				
				print trackFile
				numTrackTimes += 1
				
				# Get Individual cell metadata
				cells = lines[32::5]
				numCells = len(cells)
				
				for cell in cells:
					cell = cell.split()
					cellID = totNumCells
					stormCells[cellID] = {'time':fileDate, 'lat':float(cell[0]), 'lon':float(cell[1]), 'latr':float(cell[3]), 
											'lonr':float(cell[4]), 'orientation':float(cell[8]), 'track':str(cell[9]) + '_' + str(fileDate.date())} 
					totNumCells += 1
						
	return [stormCells, totNumCells, numTrackTimes, dates]
	


def readSegmotion(inDir, inSuffix, startTime, endTime):
	"""
	Parses raw segmotion .xml files
	
	Parameters
	----------
	inDir : string
		The input directory
	inSuffix : string
		The lowest subdirectory of the input directory
	startTime : string
		The earliest time to process
	endTime : string
		The latest time to process
		
	Returns
	-------
	list
		[stormCells, totNumCells, numTrackTimes, dates] - List containing dicts of all 
		storm cells, the number of cells, the number of files, and all valid dates
		
	"""
	
	numTrackTimes = 0
	totNumCells = 0
	stormCells = {}
	dates = [] 
	
	# Read in Segmotion files
	for root, dirs, files in os.walk(inDir):
		if inSuffix != '' and not (files and not dirs and os.path.split(root)[-1] == inSuffix): continue
		for trackFile in files:
			if trackFile.endswith('.xml'):
				
				# Skip hidden files
				if trackFile.startswith('._'): continue
				
				# Check if file falls in date range
				try:
					fileDate = datetime.datetime.strptime(str(trackFile).split('.')[0], '%Y%m%d-%H%M%S')
				except ValueError:
					print 'File ' + str(trackFile) + ' has an invalid name.  Expected format YYYYMMDD-hhmmss.xml...'
					continue
				if not startTime <= fileDate < endTime:
					continue
				if fileDate.date() not in dates: dates.append(fileDate.date())
				
				# Open file
				f = open(root + '/' + trackFile)
				lines = BeautifulSoup(f, 'html.parser').find_all('datacolumn')
				f.close()
				
				print trackFile
				numTrackTimes += 1
				
				numCells = len(lines[2].find_all('item'))
				
				for i in range(0, numCells):
					time = fileDate
					latr = float(str(lines[4].find_all('item')[i]).split('"')[1])
					lat = float(str(lines[5].find_all('item')[i]).split('"')[1])
					lonr = float(str(lines[6].find_all('item')[i]).split('"')[1])
					lon = float(str(lines[7].find_all('item')[i]).split('"')[1])
					orientation = float(str(lines[12].find_all('item')[i]).split('"')[1])
					track = str(lines[13].find_all('item')[i]).split('"')[1]
					
					cellID = totNumCells
					stormCells[cellID] = {'time': time, 'latr': latr, 'lat': lat, 'lonr': lonr, 'lon': lon, 
											'orientation': orientation, 'track': track + '_' + str(fileDate.date())}
					totNumCells += 1
														
						
	return [stormCells, totNumCells, numTrackTimes, dates]
	


def readProbSevere(inDir, inSuffix, startTime, endTime):
	"""
	Parses probSevere .ascii files
	
	Parameters
	----------
	inDir : string
		The input directory
	inSuffix : string
		The lowest subdirectory of the input directory
	startTime : string
		The earliest time to process
	endTime : string
		The latest time to process
		
	Returns
	-------
	list
		[stormCells, totNumCells, numTrackTimes, dates] - List containing dicts of all 
		storm cells, the number of cells, the number of files, and all valid dates
		
	"""
	
	numTrackTimes = 0
	totNumCells = 0
	stormCells = {} 
	dates = []
	
	# Read in Segmotion files
	for root, dirs, files in os.walk(inDir):
		if inSuffix != '' and not (files and not dirs and os.path.split(root)[-1] == inSuffix): continue
		for trackFile in files:
			if trackFile.endswith('.ascii'):
				
				# Skip hidden files
				if trackFile.startswith('._'): continue
				
				# Check if file falls in date range
				try:
					date = str(trackFile).split('.')[0].split('_')[3]
					time = str(trackFile).split('.')[0].split('_')[4]
					fileDate = datetime.datetime.strptime(date + '_' + time, '%Y%m%d_%H%M%S')
				except ValueError:
					print 'File ' + str(trackFile) + ' has an invalid name.  Expected format SSEC_AWIPS_PROBSEVERE_YYYYMMDD_hhmmss.ascii...'
					continue
				if not startTime <= fileDate < endTime:
					continue
					
				if fileDate.date() not in dates: dates.append(fileDate.date())
					
				# Open file
				f = open(root + '/' + trackFile)
				lines = f.readlines()
				f.close()
				
				print trackFile
				numTrackTimes += 1
				
				for line in lines:
					if line.startswith('Valid:'): continue
					data = str(line).split(':')
					lats = map(float, data[7].split(',')[0::2])
					lons = map(float, data[7].split(',')[1::2])
					track = data[8]						
					
					latr = (max(lats) - min(lats)) / 2.
					lonr = abs(max(lons) - min(lons)) / 2.
					
					# Calculate centroid
					points = []
					for i in range(0, len(lats)):
						points.append((lons[i], lats[i]))
					poly = Polygon(points)
					
					lon = poly.centroid.x
					lat = poly.centroid.y
					
					cellID = totNumCells
					stormCells[cellID] = {'time': fileDate, 'latr': latr, 'lat': lat, 'lonr': lonr, 'lon': lon, 
											'orientation': 'NaN', 'track': track + '_' + str(fileDate.date())}
					totNumCells += 1
						
						
	return [stormCells, totNumCells, numTrackTimes, dates]
	
	
	
	
	
	
	
