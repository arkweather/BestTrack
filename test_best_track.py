import argparse
import socket
import sys
import os
import copy
import datetime
from datetime import timedelta
import time
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.stats.mstats as stats
from collections import defaultdict

## @name Mapping constants
## @{
MIN_LAT = 20
MAX_LAT = 55
MIN_LON = 230
MAX_LON = 300
NUM_XTICKS = 8
NUM_YTICKS = 10

BORDER_WIDTH = 1
BEFORE_WIDTH = 3
AFTER_WIDTH = 1
FONT_SIZE = 12
# TITLE_FONT_SIZE = 8
BORDER_COLOR = [0/255., 0/255., 0/255.]
WATER_COLOR = [0/255., 0/255., 254/255.]
BEFORE_COLOR = [100/255., 100/255., 100/255.]

IMAGE_RES = 600
IMAGE_RES_STR = str(IMAGE_RES)
## @}

## @name Other constants
## @{
TOLERANCE = 1e-9
MAX_MISSING = 10
DASHES = '\n' + '-' * 80 + '\n\n'
STARS = '\n' + '*' * 80 + '\n\n'
## @}



## Retrieve the user-speficified command line arguments
def getOptions():
	
	# Define legal command arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('start_time', type = str, metavar = 'start_time', help = 'Start time in yyyy-mm-dd-hhmmss, or yyyy-mm-dd, etc')
	parser.add_argument('end_time', type = str, metavar = 'end_time', help = 'End time in yyyy-mm-dd-hhmmss, or yyyy-mm-dd, etc')
	parser.add_argument('-td', '--track_dir', type = str, metavar = '', default = 'segmotion_files4david', help = 'Location of source files')
	parser.add_argument('-tsf', '--track_suffix', type = str, metavar = '', default = 'smooth02_30dBZ_best', help = 'Name of last subdirectory for source files')
	parser.add_argument('-m', '--map', action = 'store_true', help = 'Toggle map creation')
	parser.add_argument('-lat', '--map_lat', type = float, nargs = '*', metavar = '', default = None, help = 'A list of latitude ranges seperated by a space. 22 25 26 27 would produce ' +
																											'ranges [22, 25], [26, 27]')
	parser.add_argument('-lon', '--map_lon', type = float, nargs = '*', metavar = '', default = None, help = 'A list of longitude (E) ranges seperated by a space. 250 255 260 270 would ' +
																											'produce ranges [250, 255], [260, 270]')
	parser.add_argument('-md', '--map_dir', type = str, metavar = '', default = 'maps', help = 'Output directory for maps')
	
	args = parser.parse_args()
	return args
	
## Check the user-specified command line arguments for errors not handled by argparse.
## Errors will print to console before terminating the script.
## @param args A dictionary of user-specified arguments
def checkArgs(args):
	
	# TODO: Add checks for input directories
	
	startTime = args['start_time']
	endTime = args['end_time']
	inSuffix = args['track_suffix']
	mapResults = args['map']
	lats = args['map_lat']
	lons = args['map_lon']
	mapDir = args['map_dir']
	
	# Time Checks
	stimeDetail = len(startTime.split('-'))
	try:
		if stimeDetail == 1:
			datetime.datetime.strptime(startTime, '%Y')
		elif stimeDetail == 2:
			datetime.datetime.strptime(startTime, '%Y-%m')
		elif stimeDetail == 3:
			datetime.datetime.strptime(startTime, '%Y-%m-%d')
		elif stimeDetail == 4:
			datetime.datetime.strptime(startTime, '%Y-%m-%d-%H%M%S')
		else: raise ValueError
	except ValueError:
		print '\nERROR: Invalid start time! Times must be formatted as YYYY-MM-DD-HHMMSS, YYYY-MM-DD, YYYY-MM, or YYYY\n'
		sys.exit(2)
		
	etimeDetail = len(endTime.split('-'))
	try:
		if etimeDetail == 1:
			datetime.datetime.strptime(endTime, '%Y')
		elif etimeDetail == 2:
			datetime.datetime.strptime(endTime, '%Y-%m')
		elif etimeDetail == 3:
			datetime.datetime.strptime(endTime, '%Y-%m-%d')
		elif etimeDetail == 4:
			datetime.datetime.strptime(endTime, '%Y-%m-%d-%H%M%S')
		else: raise ValueError
	except ValueError:
		print '\nERROR: Invalid end time! Times must be formatted as YYYY-MM-DD-hhmmss, YYYY-MM-DD, YYYY-MM, or YYYY\n'
		sys.exit(2)
	
	# Everything else
	if '\\' in inSuffix or '/' in inSuffix:
		print '\nERROR: Input directory suffix must not contain / or \\.  Instead got: ' + inSuffix + '\n'
		sys.exit(2)
	else: print 'Name of last subdirectory for original tracking files:  ' + inSuffix
	
	
	# Handle map creation variables
	if mapResults:
		
		# Latitude ranges
		if lats == None or lats == []:
			lats = [MIN_LAT, MAX_LAT]
		elif len(lats) % 2 != 0:
			print '\nInvalid number of latitudes.  There must be an even number of latitudes. Instead ' + str(len(lats)) + ' were given.\n'
			sys.exit(2)
		for lat in lats:
			if lat < MIN_LAT or lat > MAX_LAT:
				print '\nERROR: Latitude must be in range [' + str(MIN_LAT) + ', ' + str(MAX_LAT) + '].  Instead got: ' + lat + '\n'
				sys.exit(2)
			elif lats.index(lat) % 2 == 0 and abs(lat - lats[lats.index(lat) + 1]) <= TOLERANCE:
				print '\nERROR: Each set of ranges must contain different values.  One row contains the same value twice.\n'
				sys.exit(2)
				
		print '\nThis will create maps for the following latitude ranges: '
		for i in range (0, len(lats), 2):
			print str(lats[i]) + '\t' + str(lats[i+1])
		
		
		# Longitude ranges
		if lons == None or lons == []:
			lons = [MIN_LON, MAX_LON]
		elif len(lons) % 2 != 0:
			print '\nERROR: Invalid number of longitudes.  There must be an even number of longitudes. Instead ' + str(len(lons)) + ' were given.\n'
			sys.exit(2)
		elif len(lons) != len(lats):
			print '\nERROR: The number of longitudes must match the number of latitudes.\n'
			sys.exit(2)
		for lon in lons:
			if lon < MIN_LON or lon > MAX_LON:
				print '\nERROR: Longitude must be in range [' + str(MIN_LON) + ', ' + str(MAX_LON) + '].  Instead got: ' + lon + '\n'
				sys.exit(2)
			elif lons.index(lon) % 2 == 0 and abs(lon - lons[lons.index(lon) + 1]) <= TOLERANCE:
				print '\nERROR: Each set of ranges must contain different values.  One row contains the same value twice.\n'
				sys.exit(2)
				
		print '\nThis will create maps for the following longitude ranges: '
		for i in range (0, len(lons), 2):
			print str(lons[i]) + '\t' + str(lons[i+1])
			
		
		# Map output directory
		if not os.path.isdir(mapDir):
			print '\nERROR: ' + mapDir + ' does not exist or is not a directory.'
			sys.exit(2)
		else: print '\nMaps will be saved in ' + mapDir
		
	else: print 'Mapping disabled'
	

## Groups cells with the same track ID into a track dictionary
## @param stormCells Dictionary of storm cells
## @returns stormTracks Dictionary of storm tracks containing storm cells {ID: [cells]}
def find_clusters(stormCells):
	stormTracks = {}
	for cell in stormCells:
		track = stormCells[cell]['track']
		if track in stormTracks:
			stormTracks[track]['cells'].append(stormCells[cell])
		else:
			stormTracks[track] = {'cells':[stormCells[cell]]}			
	
	return stormTracks


if __name__ == '__main__':
	args = vars(getOptions())
	#print args
	
	# Set Hostname
	hostname = socket.gethostname().split('.')[0]
	print '\n\nSetting hostname to ' + hostname
	print 'Current working directory: ' + os.getcwd() + '\n'
	
	# Check user input.  Type casting is handled by argparse.
	checkArgs(args)
	
	## @name Args
	## Assuming the args check out, save their values here
	## @{
	startTime = args['start_time']
	endTime = args['end_time']
	trackDir = args['track_dir']
	inSuffix = args['track_suffix']
	mapResults = args['map']
	lats = args['map_lat']
	lons = args['map_lon']
	mapDir = args['map_dir']
	## @}
	
	# If the times check out, convert to datetime objects
	stimeDetail = len(startTime.split('-'))
	if stimeDetail == 1:
		startTime = datetime.datetime.strptime(startTime, '%Y')
	elif stimeDetail == 2:
		startTime = datetime.datetime.strptime(startTime, '%Y-%m')
	elif stimeDetail == 3:
		startTime = datetime.datetime.strptime(startTime, '%Y-%m-%d')
	elif stimeDetail == 4:
		startTime = datetime.datetime.strptime(startTime, '%Y-%m-%d-%H%M%S')
		
	etimeDetail = len(endTime.split('-'))
	if etimeDetail == 1:
		endTime = datetime.datetime.strptime(endTime, '%Y')
	elif etimeDetail == 2:
		endTime = datetime.datetime.strptime(endTime, '%Y-%m')
	elif etimeDetail == 3:
		endTime = datetime.datetime.strptime(endTime, '%Y-%m-%d')
	elif etimeDetail == 4:
		endTime = datetime.datetime.strptime(endTime, '%Y-%m-%d-%H%M%S')
	
	print '\nRunning for times ' + str(startTime) + ' through ' + str(endTime)
	
	print STARS
	
	# Check for root directory:
	print 'Reading files:'
	if not os.path.isdir(trackDir):
		print '\nERROR: Unable to find source directory "' + trackDir + '". \nIf using a relative path, please check your working directory.\n'
		sys.exit(2)
		
	## @name Storm Track Variables
	## @{
	numTrackTimes = 0
	totNumCells = 0
	stormCells = {} 
	## @}
	
	# Read in files
	for root, dirs, files in os.walk(trackDir):
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
					
					print trackFile
					numTrackTimes += 1
					
					# Get Individual cell metadata
					cells = lines[32::5]
					numCells = len(cells)
					
					for cell in cells:
						cell = cell.split()
						cellID = totNumCells
						stormCells[cellID] = {'time':fileDate, 'lat':float(cell[0]), 'lon':float(cell[1]), 'track':str(cell[9]) + '_' + str(fileDate.date())} 
						totNumCells += 1
					
	print '\nNumber of files: ' + str(numTrackTimes)
	print 'Number of storm cells: ' + str(totNumCells) + '\n'
	
	if numTrackTimes == 0:
		print 'No valid files found for this time period.  Please check the source directory and specified dates.\n'
		sys.exit(0)
		
	if mapResults:
		print 'Preparing to plot maps...'
		
		stormTracks = find_clusters(stormCells)
		
		# Handle empty specifications
		if lats == None or lats == []:
			lats = [MIN_LAT, MAX_LAT]
			
		if lons == None or lons == []:
			lons = [MIN_LON, MAX_LON]
			
		
		# Generate each map
		for i in range(0, len(lats), 2):
			print 'Plotting figure ' + str((i / 2) + 1) + ' of ' + str(len(lats) / 2)
			
			fig = plt.figure((i / 2) + 1)
			
			theseLats = [lats[i], lats[i+1]]
			theseLons = [lons[i], lons[i+1]]
			
			meanLat = np.mean(theseLats)
			meanLon = np.mean(theseLons)
			
			m = Basemap(llcrnrlon = min(theseLons), llcrnrlat = min(theseLats), urcrnrlon = max(theseLons), urcrnrlat = max(theseLats), 
						projection = 'aeqd', lat_0 = meanLat, lon_0 = meanLon)
            
			# Read in shapefiles
			m.readshapefile('States_Shapefiles/s_11au16', name = 'states', drawbounds = True)
			m.readshapefile('province/province', name = 'canada', drawbounds = True)			
            
				
			# Sort cells in each track by time and then get lat lon pairs for each cell
			for track in stormTracks:
				times = []
				finalCellsX = []
				finalCellsY = []
				for cell in stormTracks[track]['cells']:
					times.append(cell['time'])
				#print times
				#break
				times = sorted(times)
				for cellTime in times:
					for cell in stormTracks[track]['cells']:
						if cell['time'] == cellTime:
							finalCellsX.append(m(cell['lon'], cell['lat'])[0])
							finalCellsY.append(m(cell['lon'], cell['lat'])[1])
							break
							
				m.plot(finalCellsX, finalCellsY, color = 'r', linewidth = 1)			
			
			plt.show()


