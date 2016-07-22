import sys
import os
import datetime
from bs4 import BeautifulSoup

def read(ftype, inDir, inSuffix, startTime, endTime):
	if ftype == 'ryan': return readRyan(inDir, inSuffix, startTime, endTime)	
	elif ftype == 'segmotion': return readSegmotion(inDir, inSuffix, startTime, endTime)
	
def readRyan(inDir, inSuffix, startTime, endTime):
	numTrackTimes = 0
	totNumCells = 0
	stormCells = {} 
	
	# Read in Ryan files
	for root, dirs, files in os.walk(inDir):
		if files and not dirs and os.path.split(root)[-1] == inSuffix:
			for trackFile in files:
				if trackFile.endswith('.data'):
					
					# Check if file falls in date range
					try:
						fileDate = datetime.datetime.strptime(str(trackFile).split('_')[0], '%Y-%m-%d-%H%M%S')
					except ValueError:
						print 'File ' + str(trackFile) + ' has an invalid name.  Expected format YYYY-MM-DD-hhmmss_...'
						continue
					if not startTime <= fileDate < endTime:
						continue
					
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
												'lonr':float(cell[4]), 'orientation':float(cell[8]), 'track':str(cell[9]) + '_' + str(fileDate.date()), 
												'refl':float(cell[5])} 
						totNumCells += 1
						
	return [stormCells, totNumCells, numTrackTimes]
	
	
def readSegmotion(inDir, inSuffix, startTime, endTime):
	numTrackTimes = 0
	totNumCells = 0
	stormCells = {} 
	
	# Read in Segmotion files
	for root, dirs, files in os.walk(inDir):
		if files and not dirs and os.path.split(root)[-1] == inSuffix:
			for trackFile in files:
				if trackFile.endswith('.xml'):
					
					# Check if file falls in date range
					try:
						fileDate = datetime.datetime.strptime(str(trackFile).split('.')[0], '%Y%m%d-%H%M%S')
					except ValueError:
						print 'File ' + str(trackFile) + ' has an invalid name.  Expected format YYYYMMDD-hhmmss.xml...'
						continue
					if not startTime <= fileDate < endTime:
						continue
					
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
						
										
						
	return [stormCells, totNumCells, numTrackTimes]
	
	
	
	
	
	
	
	
	
