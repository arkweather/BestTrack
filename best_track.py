## @package best_track
# Concatenates storm tracks from w2segmotionll.
#
# This package is approximately equivalent to w2besttrack with the 
# potential for additional features and greater flexibility.
# This python version was converted from Ryan Lagerquist's MATLAB code
# ryan_best_tracks.m and associated files.
#
# @author David Harrison
# @date July 2016

import argparse
import socket
import sys
import os
import json
import copy
import datetime
from datetime import timedelta
import time
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.stats.mstats as stats
from collections import defaultdict
import readCells

## @name Best-track constants
## @{ 
MAX_BUFFER_DIST = 20     # Buffer distance [km].  0.1 deg in w2besttrack.
MAX_BUFFER_TIME = 21	 # Buffer time [min].  10 min in w2besttrack.
MAX_JOIN_TIME = 21		 # Buffer time for joining Theil-Sen trajectories [min].  15 min in w2besttrack.
MAX_JOIN_DIST = 70		 # Buffer distance for joining Theil-Sen trajectories [km].
MIN_MIN_CELLS = 2        # Min # storm cells per track.
MAX_MIN_CELLS = 12       # Min # storm cells per track.
MIN_ITERS = 3            # Number of outside iterations.
MAX_ITERS = 25           # Number of outside iterations.
MIN_BREAKUP_ITERS = 1    # Number of break-up iterations.
MAX_BREAKUP_ITERS = 5    # Number of break-up iterations.
## @}

## @name Mapping constants
## @{
MIN_LAT = 20
MAX_LAT = 55
MIN_LON = 230
MAX_LON = 300
NUM_XTICKS = 8
NUM_YTICKS = 10

BEFORE_WIDTH = 4
AFTER_WIDTH = 2
FONT_SIZE = 12

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
	parser.add_argument('-i', '--input_dir', type = str, metavar = '', default = 'segmotion_files4david', help = 'Location of source files')
	parser.add_argument('-s', '--dir_suffix', type = str, metavar = '', default = 'smooth02_30dBZ', help = 'Name of last subdirectory for source files')
	parser.add_argument('-t', '--type', type = str, metavar = '', default = 'ryan', help = 'Type of input data: segmotion (.xml), probsevere (.ascii), or ryan (.data)')
	parser.add_argument('-bd', '--buffer_dist', type = float, metavar = '', default = 10., help = 'Buffer distance between storm cell and Theil-Sen trajectory (km)')
	parser.add_argument('-bt', '--buffer_time', type = float, metavar = '', default = 11., help = 'Buffer time for joining two Theil-Sen trajectories (min)')
	parser.add_argument('-jt', '--join_time', type = float, metavar = '', default = 16., help = 'Time threshold to join two or more storm tracks (min)')
	parser.add_argument('-jd', '--join_dist', type = float, metavar = '', default = 50., help = 'Distance threshold to join two or more storm tracks (km)')
	parser.add_argument('-mc', '--min_cells', type = int, metavar = '', default = 3, help = 'Minimum number of storm cells per track')
	parser.add_argument('-mi', '--main_iters', type = int, metavar = '', default = 5, help = 'Number of main iterations')
	parser.add_argument('-bi', '--breakup_iters', type = int, metavar = '', default = 3, help = 'Number of breakup iterations')
	parser.add_argument('-o', '--out_dir', type = str, metavar = '', default = 'tracks', help = 'Name of output directory for new tracking files')
	parser.add_argument('-ts', '--time_step', action = 'store_true', help = 'Toggle file creation for each time step. Default is to combine all times into one file.')
	parser.add_argument('-m', '--map', action = 'store_true', help = 'Toggle map creation')
	parser.add_argument('-lat', '--map_lat', type = float, nargs = '*', metavar = '', default = None, help = 'A list of latitude ranges seperated by a space. 22 25 26 27 would produce ' +
																											'ranges [22, 25], [26, 27] for mapping')
	parser.add_argument('-lon', '--map_lon', type = float, nargs = '*', metavar = '', default = None, help = 'A list of longitude (E) ranges seperated by a space. 250 255 260 270 would ' +
																											'produce ranges [250, 255], [260, 270] for mapping')
	parser.add_argument('-md', '--map_dir', type = str, metavar = '', default = 'maps', help = 'Output directory for maps')
	
	args = parser.parse_args()
	return args
	
## Check the user-specified command line arguments for errors not handled by argparse.
## Errors will print to console before terminating the script.
## @param args A dictionary of user-specified arguments
def checkArgs(args):
	
	startTime = args['start_time']
	endTime = args['end_time']
	inSuffix = args['dir_suffix']
	fType = args['type']
	bufferDist = args['buffer_dist']
	bufferTime = args['buffer_time']
	joinTime = args['join_time']
	joinDist = args['join_dist']
	minCells = args['min_cells']
	mainIters = args['main_iters']
	breakIters = args['breakup_iters']
	outDir = args['out_dir']
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
	
	types = ['segmotion', 'probsevere', 'ryan']
	if fType not in types:
		print 'ERROR: Invalid file type specified. Expected segmotion, probsevere, or ryan.  Instead got: ' + fType + '\n'
	else: print 'Data file type: ' + fType
	
	if bufferDist <= 0 or bufferDist > MAX_BUFFER_DIST:
		print '\nERROR: Buffer distance must be in range (0, ' + str(MAX_BUFFER_DIST) + '].  Instead got: ' + str(bufferDist) + '\n'
		sys.exit(2)
	else: print 'Buffer distance between storm cell and Theil-Sen trajectory:  ' + str(bufferDist) + ' km'
	
	if bufferTime <= 0 or bufferTime > MAX_BUFFER_TIME:
		print '\nERROR: Buffer time must be in range (0, ' + str(MAX_BUFFER_TIME) + '].  Instead got: ' + str(bufferTime) + '\n'
		sys.exit(2)
	else: print 'Buffer time for joining two Theil-Sen trajectories:  ' + str(bufferTime) + ' min'
	
	if joinTime < bufferTime or joinTime > MAX_JOIN_TIME:
		print '\nERROR: Join time must be in range [' + str(bufferTime) + ', ' + str(MAX_JOIN_TIME) + '].  Instead got: ' + str(joinTime) + '\n'
		sys.exit(2)
	else: print 'Join time:  ' + str(joinTime) + ' min'
	
	if joinDist < bufferDist or joinDist > MAX_JOIN_DIST:
		print '\nERROR: Join distance must be in range [' + str(bufferDist) + ', ' + str(MAX_JOIN_DIST) + '].  Instead got: ' + str(joinDist) + '\n'
		sys.exit(2)
	else: print 'Join Distance:  ' + str(joinDist) + ' km'
	
	if minCells < MIN_MIN_CELLS or minCells > MAX_MIN_CELLS:
		print '\nERROR: Min Cells must be in range ['  + str(MIN_MIN_CELLS) + ', ' + str(MAX_MIN_CELLS) + '].  Instead got: ' + str(minCells) + '\n'
		sys.exit(2)
	else: print 'Minimum number of cells per track:  ' + str(minCells)
	
	if mainIters < MIN_ITERS or mainIters > MAX_ITERS:
		print '\nERROR: Number of main iterations must be in range ['  + str(MIN_ITERS) + ', ' + str(MAX_ITERS) + '].  Instead got: ' + str(mainIters) + '\n'
		sys.exit(2)
	else: print 'Number of main iterations:  ' + str(mainIters)
	
	if breakIters < MIN_BREAKUP_ITERS or breakIters > MAX_BREAKUP_ITERS:
		print '\nERROR: Number of breakup iterations must be in range ['  + str(MIN_BREAKUP_ITERS) + ', ' + str(MAX_BREAKUP_ITERS) + '].  Instead got: ' + str(breakIters) + '\n'
		sys.exit(2)
	else: print 'Number of breakup iterations:  ' + str(breakIters)
	
	if not os.path.isdir(outDir):
		print 'Unable to locate output directory. The specified location will be created.'
		os.makedirs(outDir)
	
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

## Computes the Theil-Sen fit for a single storm track.
## See theil_sen_batch() for more detail.
## @param track The value of a single track within the storm track dictionary 
## @returns a storm track dict value with updated items for the provided track
def theil_sen_single(track):
	
	times = []
	x = []
	y = []
	for cell in track['cells']:
		times.append(time.mktime(cell['time'].timetuple())) # Converts datetime object to seconds since epoch time
		x.append(cell['x'])
		y.append(cell['y'])
		
	if len(np.unique(times)) > 1:
		theilSenDataX = stats.theilslopes(x, times)
		theilSenDataY = stats.theilslopes(y, times)
		
		track['u'] = theilSenDataX[0]
		track['v'] = theilSenDataY[0]
		track['t0'] = datetime.datetime.fromtimestamp(min(times))
		track['tend'] = datetime.datetime.fromtimestamp(max(times))
		track['x0'] = theilSenDataX[1] + theilSenDataX[0] * (min(times))
		track['y0'] = theilSenDataY[1] + theilSenDataY[0] * (min(times))
		track['xf'] = theilSenDataX[1] + theilSenDataX[0] * (max(times))
		track['yf'] = theilSenDataY[1] + theilSenDataY[0] * (max(times)) 
		
	else:
		stormTracks[track]['u'] = 'NaN'
		stormTracks[track]['v'] = 'NaN'
		track['t0'] = datetime.datetime.fromtimestamp(min(times))
		track['tend'] = datetime.datetime.fromtimestamp(max(times))
		track['x0'] = track['cells'][times.index(min(times))]['x']
		track['y0'] = track['cells'][times.index(min(times))]['y']
		track['xf'] = track['cells'][times.index(max(times))]['x']
		track['yf'] = track['cells'][times.index(max(times))]['y']
		
	return track

## Computes the Theil-Sen fit for each storm track.
## Sources: http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mstats.theilslopes.html 
##          https://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
## @param stormTracks A dictionary of track IDs each containing associated storm cells with Lat, Lon, X, Y, and datetime params
##                  {'ID':{['x', 'y', 'lat', 'lon', 'times', 'track']}}
## @returns stormTracks Modified with new values 'u', 'v', 't0', 'tend', 'x0', 'y0', 'xf', 'yf'	
def theil_sen_batch(stormTracks):
	
	for track in stormTracks:
		times = []
		x = []
		y = []
		for cell in stormTracks[track]['cells']:
			times.append(time.mktime(cell['time'].timetuple())) # Converts datetime object to seconds since epoch time
			x.append(cell['x'])
			y.append(cell['y'])
			
		#print times
		#print cellTracks[track]['x']
		
		if len(np.unique(times)) > 1:
			theilSenDataX = stats.theilslopes(x, times)
			theilSenDataY = stats.theilslopes(y, times)
			
			stormTracks[track]['u'] = theilSenDataX[0]
			stormTracks[track]['v'] = theilSenDataY[0]
			stormTracks[track]['t0'] = datetime.datetime.fromtimestamp(min(times))
			stormTracks[track]['tend'] = datetime.datetime.fromtimestamp(max(times))
			stormTracks[track]['x0'] = theilSenDataX[1] + theilSenDataX[0] * (min(times))
			stormTracks[track]['y0'] = theilSenDataY[1] + theilSenDataY[0] * (min(times))
			stormTracks[track]['xf'] = theilSenDataX[1] + theilSenDataX[0] * (max(times))
			stormTracks[track]['yf'] = theilSenDataY[1] + theilSenDataY[0] * (max(times)) 
			
		else:
			stormTracks[track]['u'] = 'NaN'
			stormTracks[track]['v'] = 'NaN'
			stormTracks[track]['t0'] = datetime.datetime.fromtimestamp(min(times))
			stormTracks[track]['tend'] = datetime.datetime.fromtimestamp(max(times))
			stormTracks[track]['x0'] = stormTracks[track]['cells'][times.index(min(times))]['x']
			stormTracks[track]['y0'] = stormTracks[track]['cells'][times.index(min(times))]['y']
			stormTracks[track]['xf'] = stormTracks[track]['cells'][times.index(max(times))]['x']
			stormTracks[track]['yf'] = stormTracks[track]['cells'][times.index(max(times))]['y']
	
	return stormTracks
		
	
######################################################################################################################
#                                                                                                                    #
#  Main Method - Handle user input, read in files, then run calculations                                             #
#                                                                                                                    #
######################################################################################################################

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
	inDir = args['input_dir']
	inSuffix = args['dir_suffix']
	fType = args['type']
	bufferDist = args['buffer_dist']
	bufferTime = timedelta(minutes = int(args['buffer_time']))
	joinTime = timedelta(minutes = int(args['join_time']))
	joinDist = args['join_dist']
	minCells = args['min_cells']
	mainIters = args['main_iters']
	breakIters = args['breakup_iters']
	outDir = args['out_dir']
	outType = args['time_step']
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
	
	####################################################################################################################
	#                                                                                                                  #
	#  Read in the files and process data                                                                              #
	#                                                                                                                  #
	####################################################################################################################
	
	# Check for root directory:
	print 'Reading files:'
	if not os.path.isdir(inDir):
		print '\nERROR: Unable to find source directory "' + inDir + '". \nIf using a relative path, please check your working directory.\n'
		sys.exit(2)
		
	data = readCells.read(fType, inDir, inSuffix, startTime, endTime)
	stormCells = data[0]
	totNumCells = data[1]
	numTrackTimes = data[2]
					
	print '\nNumber of files: ' + str(numTrackTimes)
	print 'Number of storm cells: ' + str(totNumCells) + '\n'
	
	if numTrackTimes == 0:
		print 'No valid files found for this time period.  Please check the source directory and specified dates.\n'
		sys.exit(0)
		
	# Project onto equidistant coord system
	print 'Projecting storm cells onto equidistant coordinate system...'
	
	## @name Projection variables
	## @{
	meanLat = np.mean([MIN_LAT, MAX_LAT])
	meanLon = np.mean([MIN_LON, MAX_LON])
	xyDistMax = 0
	llDistMax = 0
	distanceRatio = 0
	
	# Setup equidistant map projection
	m = Basemap(llcrnrlon = MIN_LON, llcrnrlat = MIN_LAT, urcrnrlon = MAX_LON, urcrnrlat = MAX_LAT, 
						projection = 'aeqd', lat_0 = meanLat, lon_0 = meanLon)
    ## @}
    
	for cell in stormCells:
		stormCells[cell]['x'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[0]
		stormCells[cell]['y'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[1]
	
	# Find ratio between x-y distances and lat-lon distances
	xMin, yMin = m(MIN_LON, MIN_LAT)
	xMax, yMax = m(MAX_LON, MAX_LAT)
	
	xyDistMax = np.sqrt((xMin - xMax)**2 + (yMin - yMax)**2)
	
	# Find distance between two lat lon coordinates
	# Source: https://en.wikipedia.org/wiki/Great-circle_distance
	# point1 = [MAX_LON, MIN_LAT]
	# point2 = [MIN_LON, MAX_LAT]
	
	rlat1 = np.radians(MIN_LAT)
	rlat2 = np.radians(MAX_LAT)
	r = 6371 # Mean radius of Earth (km)
	dlon = abs(MAX_LON - MIN_LON)
	dsig = np.arccos(np.sin(rlat1) * np.sin(rlat2) + np.cos(rlat1) * np.cos(rlat2) * np.cos(np.radians(dlon)))
	llDistMax = r * dsig
	
	distanceRatio = llDistMax / xyDistMax
	
	print 'Ratio between x-y distances and lat-lon distances: ' + str(distanceRatio)
	print DASHES
	
	####################################################################################################################
	#                                                                                                                  #
	#  Calculations!                                                                                                   #
	#                                                                                                                  #
	####################################################################################################################
	
	print 'Beginning Calculations...'
	REPORT_EVERY = 1000
	scOrigin = copy.deepcopy(stormCells)
	
	# Main iterations
	for i in range(0, mainIters):
		print '\nMain iteration: ' + str(i + 1)
		
		anyChanges = False
		
		# Breakup iterations
		for j in range(0, breakIters):
			
			print '\nBreakup iteration: ' + str(j + 1)
			print 'Finding clusters...'
			stormTracks = find_clusters(stormCells)
			print 'Number of clusters: ' + str(len(stormTracks))
			
			print 'Computing Theil-Sen fit for each cluster...'	
			stormTracks = theil_sen_batch(stormTracks)
			
			# Assign cells to nearest cluster
			print 'Assigning each cell to nearest cluster...'
			changedCells = 0
			for cell in stormCells:
				cellTime = stormCells[cell]['time']
				cellX = stormCells[cell]['x']
				cellY = stormCells[cell]['y']
				
				# Calculate distances
				minDist = 1e9
				minTrack = stormTracks[min(stormTracks)]
				for track in stormTracks:
					# Only compare to tracks in temporal range
					if not (stormTracks[track]['t0'] - bufferTime <= cellTime <= stormTracks[track]['tend'] + bufferTime):
						continue
						
					# Preference individual cells to join other tracks
					if len(stormTracks[track]['cells']) < 2 and track == stormCells[cell]['track']:
						continue
						
					# TODO test
					#if len(stormTracks[track]['cells']) > 3: continue
					
					if stormTracks[track]['u'] == 'NaN':
						xPoint = stormTracks[track]['x0']
						yPoint = stormTracks[track]['y0']
					else:
						xPoint = stormTracks[track]['x0'] + (stormTracks[track]['u'] * ((cellTime - stormTracks[track]['t0']).total_seconds()))
						yPoint = stormTracks[track]['y0'] + (stormTracks[track]['v'] * ((cellTime - stormTracks[track]['t0']).total_seconds()))
					
					dist = np.sqrt((cellX - xPoint)**2 + (cellY - yPoint)**2)
					dist = dist * distanceRatio # Convert from x,y to km
					#print str(dist)
					
					# Force cells to be assigned to a track (not NaN)
					# If need be they'll be weeded out in the tie break step later
					if dist < minDist and track != 'NaN':
						minDist = dist
						minTrack = track
						
				if minDist <= bufferDist:
					if minTrack != stormCells[cell]['track']: changedCells += 1
					stormCells[cell]['track'] = minTrack
				else:
					stormCells[cell]['track'] = 'NaN'
					
				if cell % REPORT_EVERY == 0:
					print '......' + str(cell) + ' of ' + str(totNumCells) + ' assigned......'
					
			print 'All cells have been assigned!'
			print 'Number of modified cells: ' + str(changedCells)
			
			# TODO: Check for changes
			#newChanges = (not scLast == stormCells)
			#anyChanges = newChanges or anyChanges
			#if not newChanges:
			#	print 'No new changes. Ending the breakup step.'
			#	break
				
		# ------ End of breakup iteration ------ #
			
		# Find new clusters
		print '\nFinding new clusters after breakup...'
		lastNumTracks = len(stormTracks)
		stormTracks = find_clusters(stormCells)
		print 'Original number of clusters: ' + str(lastNumTracks)
		print 'New number of clusters: ' + str(len(stormTracks))
		
		# Get Theil-Sen fit
		print 'Computing Theil-Sen fit for each new cluster...'
		stormTracks = theil_sen_batch(stormTracks)
		
		# Join similar clusters
		print 'Joining similar clusters...'
		removeTracks = []
		tracks = sorted(stormTracks.keys())
		totNumTracks = len(tracks)
		
		for j in range(0, len(tracks)):
			track1 = tracks[j]
			
			# Skip tracks with only 1 cell
			if len(stormTracks[track1]['cells']) < 2:
				if j % REPORT_EVERY == 0: print '......' + str(j) + ' of ' + str(totNumTracks) + ' processed for joining......'
				continue
			if track1 == 'NaN': continue
			
			for k in range(0, j - 1):
				track2 = tracks[k]
				
				#if len(stormTracks[track2]['cells']) < 2: continue
				if track2 in removeTracks: continue
				if track2 == 'NaN': continue
				if len(stormTracks[track2]['cells']) < 2: continue
				
				# Check time gap between tracks
				if stormTracks[track1]['t0'] > stormTracks[track2]['t0']:
					earlyIndex = track2
					lateIndex = track1
				else:
					earlyIndex = track1
					lateIndex = track2
				timeDiff = stormTracks[lateIndex]['t0'] - stormTracks[earlyIndex]['tend']
				if abs(timeDiff.total_seconds()) > joinTime.total_seconds: continue
				
				# Check distance between tracks
				x1 = stormTracks[earlyIndex]['xf']
				y1 = stormTracks[earlyIndex]['yf']
				x2 = stormTracks[lateIndex]['x0']
				y2 = stormTracks[lateIndex]['y0']
				
				dist = np.sqrt((x1-x2)**2 + (y1 - y2)**2)
				dist = dist * distanceRatio
				
				# Limit track join distance
				if dist > joinDist: continue
				
				# Check velocity difference between tracks
				u1 = stormTracks[earlyIndex]['u'] * distanceRatio # Km / s
				v1 = stormTracks[earlyIndex]['v'] * distanceRatio # Km / s
				u2 = stormTracks[lateIndex]['u'] * distanceRatio  # Km / s
				v2 = stormTracks[lateIndex]['v'] * distanceRatio  # Km / s
			
				velocityDiff = np.sqrt((u1 - u2)**2 + (v1 - v2)**2)
				if velocityDiff > float(bufferDist) / bufferTime.total_seconds(): continue
				
				# Check if track predictions are close enough				
				dist = []
				for cell in stormTracks[lateIndex]['cells']:
					xActual = cell['x']
					yActual = cell['y']
				
					cellTime = cell['time']
					xPredict = stormTracks[earlyIndex]['xf'] + (stormTracks[earlyIndex]['u'] * ((cellTime - stormTracks[earlyIndex]['tend']).total_seconds()))
					yPredict = stormTracks[earlyIndex]['yf'] + (stormTracks[earlyIndex]['v'] * ((cellTime - stormTracks[earlyIndex]['tend']).total_seconds()))
				
					dist.append(np.sqrt((xPredict - xActual)**2 + (yPredict - yActual)**2) * distanceRatio)
				
				if np.mean(dist) > bufferDist: continue
				
				# If the two tracks survived the process, join them 'cause clearly they're meant to be together ;-)
				removeTracks.append(track2)
				for cell in stormTracks[track2]['cells']:
					cell['track'] = track1
					stormTracks[track1]['cells'].append(cell)
				
				stormTracks[track1] = theil_sen_single(stormTracks[track1])
					
			if j % REPORT_EVERY == 0:
				print '......' + str(j) + ' of ' + str(totNumTracks) + ' processed for joining......'
				
		print 'All tracks have been joined if necessary!'
		print 'Merged ' + str(len(removeTracks)) + ' tracks\n'
		
		# ------ End of Joining process ------ #
		
		print 'Finding new clusters after joining...'
		lastNumTracks = len(stormTracks)
		stormTracks = find_clusters(stormCells)
		stormTracks = theil_sen_batch(stormTracks)
		print 'Original number of clusters: ' + str(lastNumTracks)
		print 'New number of clusters: ' + str(len(stormTracks))
		
		# Break ties (multiple cells assigned to same cluster at same time step)
		print '\nBreaking ties...'
		scLast = stormCells
		count = 0
		breaks = 0
		
		for track in stormTracks:
			if len(stormTracks[track]['cells']) < 2: 
				if count % REPORT_EVERY == 0: print '......' + str(count) + ' of ' + str(totNumTracks) + ' tracks processed for ties......'
				count += 1
				continue
			
			# Map all cells to their times
			times = {}
			for cell in stormTracks[track]['cells']:
				if cell['time'] in times:
					times[cell['time']].append(cell)
				else:
					times[cell['time']] = [cell]
					
			# Get duplicate times
			for thisTime in times:
				cells = []
				if len(times[thisTime]) > 1:
					cells = times[thisTime]
					
					# Compare each cell and keep the one closest to the track
					dist = []
					for cell in cells:
						cellX = cell['x']
						cellY = cell['y']
						xPredict = stormTracks[track]['x0'] + (stormTracks[track]['u'] * ((thisTime - stormTracks[track]['t0']).total_seconds()))
						yPredict = stormTracks[track]['y0'] + (stormTracks[track]['v'] * ((thisTime - stormTracks[track]['t0']).total_seconds()))
					
						dist.append(np.sqrt((xPredict - cellX)**2 + (yPredict - cellY)**2) * distanceRatio)
					
					minCell = cells[dist.index(min(dist))]
					for cell in cells:
						if cell != minCell:
							stormTracks[track]['cells'].remove(cell)
							stormCells[stormCells.keys()[stormCells.values().index(cell)]]['track'] = 'NaN'
							
					breaks += 1
							
			if count % REPORT_EVERY == 0:
				print '......' + str(count) + ' of ' + str(totNumTracks) + ' tracks processed for ties......'
			count += 1
				
		print 'All tracks have been processed for tie breaks'
		print 'Number of tie breaks: ' + str(breaks)
					
		# ------ End of Main iteration ------ #
		
	# Remove clusters with too few cells
	print '\nRemoving clusters with too few cells...'
	numRemoved = 0
	for track in stormTracks:
		if len(stormTracks[track]['cells']) < int(minCells):
			for cell in stormTracks[track]['cells']:
				cell['track'] = 'NaN'
			numRemoved += 1
	
	print 'Number of removed tracks: ' + str(numRemoved + 1)
	lastNumTracks = len(stormTracks)
	
	print '\nPerforming final cluster identification...'
	stormTracks = find_clusters(stormCells)
	stormTracks = theil_sen_batch(stormTracks)
	stormTracks.pop('NaN', None)
	print 'Original number of clusters: ' + str(lastNumTracks)
	print 'New number of clusters: ' + str(len(stormTracks))
	
	print DASHES
	
	
	####################################################################################################################
	#                                                                                                                  #
	#  Maps!                                                                                                           #
	#                                                                                                                  #
	####################################################################################################################		
	
	if mapResults:
		print 'Preparing to plot maps...'
		
		# Get original storm tracks
		stOrigin = find_clusters(scOrigin)
		stOrigin = theil_sen_batch(stOrigin)
		
		# Handle empty specifications
		if lats == None or lats == []:
			lats = [MIN_LAT, MAX_LAT]
			
		if lons == None or lons == []:
			lons = [MIN_LON, MAX_LON]
			
		
		# Generate each map
		for i in range(0, len(lats), 2):
			print 'Plotting figure ' + str((i / 2) + 1) + ' of ' + str(len(lats) / 2) + '...'
			
			fig = plt.figure((i / 2) + 1, figsize = (1920 / 96., 1080 / 96.), dpi = 96)
			
			theseLats = [lats[i], lats[i+1]]
			theseLons = [lons[i], lons[i+1]]
			
			meanLat = np.mean(theseLats)
			meanLon = np.mean(theseLons)
			
			m = Basemap(llcrnrlon = min(theseLons), llcrnrlat = min(theseLats), urcrnrlon = max(theseLons), urcrnrlat = max(theseLats), 
						projection = 'aeqd', lat_0 = meanLat, lon_0 = meanLon)
            
			# Read in shapefiles
			m.readshapefile('States_Shapefiles/s_11au16', name = 'states', drawbounds = True)
			m.readshapefile('province/province', name = 'canada', drawbounds = True)			
            
			# Sort cells in each original track by time and then get lat lon pairs for each cell
			for track in stOrigin:
				times = []
				originCellsX = []
				originCellsY = []
				
				for cell in stOrigin[track]['cells']:
					times.append(cell['time'])
				times = sorted(times)
				for cellTime in times:
					for cell in stOrigin[track]['cells']:
						if cell['time'] == cellTime:
							originCellsX.append(m(cell['lon'], cell['lat'])[0])
							originCellsY.append(m(cell['lon'], cell['lat'])[1])
							break
				
				if len(originCellsX) < 2: m.scatter(originCellsX, originCellsY, color = 'grey', marker = 'o')		
				else: m.plot(originCellsX, originCellsY, color = 'grey', linewidth = BEFORE_WIDTH)
				
			# Sort cells in each track by time and then get lat lon pairs for each cell
			for track in stormTracks:
				times = []
				refl = []
				finalCellsX = []
				finalCellsY = []
				for cell in stormTracks[track]['cells']:
					times.append(cell['time'])
					try: refl.append(cell['refl'])
					except KeyError: refl.append(50)
				#print times
				#break
				times = sorted(times)
				for cellTime in times:
					for cell in stormTracks[track]['cells']:
						if cell['time'] == cellTime:
							finalCellsX.append(m(cell['lon'], cell['lat'])[0])
							finalCellsY.append(m(cell['lon'], cell['lat'])[1])
							break
				
				if max(refl) >= 50: color = 'r'
				elif max(refl) >= 40 and max(refl) < 50: color = 'y'
				else: color = 'g'
				m.plot(finalCellsX, finalCellsY, color = color, linewidth = AFTER_WIDTH)			
			
			plt.show()
			
			# Save map to file
			print 'Saving figure ' + mapDir + '/' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + str((i/2) + 1) + '.png' + '...'
			plt.savefig(mapDir + '/' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + str((i/2) + 1) + '.png', bbox_inches = 'tight')
		
		print DASHES
	
	####################################################################################################################
	#                                                                                                                  #
	#  Output                                                                                                          #
	#                                                                                                                  #
	####################################################################################################################
	
	# Reset basemap for conversions		
	print 'Preparing output...'
	# Setup equidistant map projection
	m = Basemap(llcrnrlon = MIN_LON, llcrnrlat = MIN_LAT, urcrnrlon = MAX_LON, urcrnrlat = MAX_LAT, 
						projection = 'aeqd', lat_0 = meanLat, lon_0 = meanLon)
	
	# Remove NaN track cells
	print 'Removing unassigned cells...'
	removeCells = []
	for cell in stormCells:
		if stormCells[cell]['track'] == 'NaN': removeCells.append(cell)
	for cell in removeCells: stormCells.pop(cell, None)
	
	print 'Finding new start time, age, and speed for each cell...'
	
	for track in stormTracks:
		# Get startime and age
		# Convert all datetimes to str for JSON
		for cell in stormTracks[track]['cells']:
			cell['start_time'] = stormTracks[track]['t0']
			cell['age'] = (cell['time'] - cell['start_time']).total_seconds()
			
		# Sort cells by time	
		times = []
		cells = []		
		for cell in stormTracks[track]['cells']:
			times.append(cell['time'])
			
		times = sorted(times)
		for cellTime in times:
			for cell in stormTracks[track]['cells']:
				if cell['time'] == cellTime:
					cells.append(cell)
					break
		
		# Calculate speed and component velocities for each cell		
		for cell in cells:		
			index = cells.index(cell)
			if index == 0: 
				cell['motion_east'] = stormTracks[track]['u'] * distanceRatio * 1000 		# m/s
				cell['motion_south'] = -1 * stormTracks[track]['v'] * distanceRatio * 1000 	# m/s
				
			else:
				prevX = cells[index - 1]['x']
				prevY = cells[index - 1]['y']
				prevTime = cells[index - 1]['time']
				
				cell['motion_east'] = (cell['x'] - prevX) / ((cell['time'] - prevTime).total_seconds()) * distanceRatio * 1000			# m/s
				cell['motion_south'] = -1 * (cell['y'] - prevY) / ((cell['time'] - prevTime).total_seconds()) * distanceRatio * 1000	# m/s
				
			cell['speed'] = np.sqrt(cell['motion_east']**2 + cell['motion_south']**2)
		
		# Cleanup for output
		ids = []
		for cell in stormTracks[track]['cells']:
			# Convert times to strings for JSON
			cell['time'] = str(cell['time'])
			cell['start_time'] = str(cell['start_time'])
			
			# Remove data specific to this run
			cell.pop('x', None)
			cell.pop('y', None)
			
			ids.append(stormCells.keys()[stormCells.values().index(cell)])
		
		# Only save cell IDs to storm track to save space
		stormTracks[track]['cells'] = ids			
			
		stormTracks[track]['t0'] = str(stormTracks[track]['t0'])
		stormTracks[track]['tend'] = str(stormTracks[track]['tend'])
			
		# Convert x, y back to lon, lat and km/s to m/s
		stormTracks[track]['lon0'], stormTracks[track]['lat0'] = m(stormTracks[track]['x0'], stormTracks[track]['y0'], inverse = True)
		stormTracks[track]['lonf'], stormTracks[track]['latf'] = m(stormTracks[track]['xf'], stormTracks[track]['yf'], inverse = True)
		stormTracks[track]['u'] = stormTracks[track]['u'] * distanceRatio * 1000 	# m/s
		stormTracks[track]['v'] = stormTracks[track]['v'] * distanceRatio * 1000 	# m/s
		
		# Remove data specific to this run
		stormTracks[track].pop('x0', None)
		stormTracks[track].pop('y0', None)
		stormTracks[track].pop('xf', None)
		stormTracks[track].pop('yf', None)
				
	if outType:
		# Print data for each time step
		
		# Sort cells by time
		times = []		
		for cell in stormCells:
			times.append(datetime.datetime.strptime(stormCells[cell]['time'], '%Y-%m-%d %H:%M:%S'))
			
		times = sorted(np.unique(times))
		
		for cellTime in times:
			cells = {}
			for cell in stormCells:
				if datetime.datetime.strptime(stormCells[cell]['time'], '%Y-%m-%d %H:%M:%S') == cellTime:
					cells[cell] = stormCells[cell]
			
			# Print stormCell data for this time step
			print 'Printing ' + str(cellTime) + '_' + 'cells.data'
			with open(outDir + '/' + str(cellTime) + '_' + 'cells.data', 'w') as outfile:
				json.dump(cells, outfile, sort_keys = True, indent = 0)
	
			outfile.close()
		
	else:
		# Print stormCells to data file
		print 'Printing ' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + 'cells.data'
		with open(outDir + '/' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + 'cells.data', 'w') as outfile:
			json.dump(stormCells, outfile, sort_keys = True, indent = 0)
	
		outfile.close()

		# Print stormTracks to data file
		print 'Printing ' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + 'tracks.data'
		with open(outDir + '/' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + 'tracks.data', 'w') as outfile:
			json.dump(stormTracks, outfile, sort_keys = True, indent = 0)
	
		outfile.close()

	print '\n\nBest Track has completed succesfully!\n\n'
	
	
