"""
 Concatenates storm tracks from w2segmotionll, probSevere, and post-processed .data (Ryan) files.

 This package is approximately equivalent to w2besttrack with the 
 potential for additional features and greater flexibility.
 This python version was converted from Ryan Lagerquist's MATLAB code
 ryan_best_tracks.m and associated files.

 Author :	David Harrison
 Date   :	July 2016

"""

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
import traceback
from multiprocessing import Pool, Manager, Value, Array, Lock
import multiprocessing
import ctypes
from contextlib import closing

# Best-track constants 
MAX_BUFFER_DIST = 20  	# Buffer distance [km].  0.1 deg in w2besttrack.
MAX_BUFFER_TIME = 21  	# Buffer time [min].  10 min in w2besttrack.
MAX_JOIN_TIME = 21  	# Buffer time for joining Theil-Sen trajectories [min].  15 min in w2besttrack.
MAX_JOIN_DIST = 70  	# Buffer distance for joining Theil-Sen trajectories [km].
MIN_MIN_CELLS = 2  		# Min min number storm cells per track.
MAX_MIN_CELLS = 12  	# Max min number storm cells per track.
MIN_ITERS = 3  			# Number of outside iterations.
MAX_ITERS = 25  		# Number of outside iterations.
MIN_BREAKUP_ITERS = 1  	# Number of break-up iterations.
MAX_BREAKUP_ITERS = 5  	# Number of break-up iterations.

# Mapping constants
MIN_LAT = 20
MAX_LAT = 51
MIN_LON = -119
MAX_LON = -62

BEFORE_WIDTH = 4
AFTER_WIDTH = 2
FONT_SIZE = 12

# Other constants
TOLERANCE = 1e-9
MAX_MISSING = 10
DASHES = '\n' + '-' * 80 + '\n\n'
STARS = '\n' + '*' * 80 + '\n\n'


def getOptions():
	"""
	Retrieve the user-speficified command line arguments
	
	Returns
	--------
	Namespace
		Namespace of parsed arguments returned by ArgumentParser.parse_args()
			
	"""

	# Load default values from config file
	try:
		f = open('best_track.config')
		lines = f.readlines()
		f.close()

		inDir = lines[9].split('=')[1].strip()
		suf = lines[11].split(' =')[1].strip()
		ftype = lines[13].split('=')[1].strip()
		bd = float(lines[15].split('=')[1])
		bt = float(lines[17].split('=')[1])
		jt = float(lines[19].split('=')[1])
		jd = float(lines[21].split('=')[1])
		mc = int(lines[23].split('=')[1])
		mi = int(lines[25].split('=')[1])
		bi = int(lines[27].split('=')[1])
		bg = int(lines[29].split('=')[1])
		outDir = lines[31].split('=')[1].strip()
		ts = bool(int(lines[33].split('=')[1]))
		m = bool(int(lines[35].split('=')[1]))

	except IOError:
		print 'Unable to find best_track.config.  Please make sure the file is in the same directory as best_track.py\n\n'
		sys.exit(2)
	
	except IndexError:
		print 'A value is missing in best_track.config or the file has become corrupted.  See the README file for more info.\n\n'
		sys.exit(2)

	# Define legal command arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('start_time', type=str, metavar='start_time',
		                help='Start time in yyyy-mm-dd-hhmmss, or yyyy-mm-dd, etc')
	parser.add_argument('end_time', type=str, metavar='end_time',
		                help='End time in yyyy-mm-dd-hhmmss, or yyyy-mm-dd, etc')
	parser.add_argument('-i', '--input_dir', type=str, metavar='', default=inDir, help='Location of source files')
	parser.add_argument('-s', '--dir_suffix', type=str, metavar='', default=suf,
		                help='Name of last subdirectory for source files')
	parser.add_argument('-t', '--type', type=str, metavar='', default=ftype,
		                help='Type of input data: segmotion (.xml), probsevere (.ascii), or ryan (.data)')
	parser.add_argument('-bd', '--buffer_dist', type=float, metavar='', default=bd,
		                help='Buffer distance between storm cell and Theil-Sen trajectory (km)')
	parser.add_argument('-bt', '--buffer_time', type=float, metavar='', default=bt,
		                help='Buffer time for joining two Theil-Sen trajectories (min)')
	parser.add_argument('-jt', '--join_time', type=float, metavar='', default=jt,
		                help='Time threshold to join two or more storm tracks (min)')
	parser.add_argument('-jd', '--join_dist', type=float, metavar='', default=jd,
		                help='Distance threshold to join two or more storm tracks (km)')
	parser.add_argument('-mc', '--min_cells', type=int, metavar='', default=mc,
		                help='Minimum number of storm cells per track')
	parser.add_argument('-mi', '--main_iters', type=int, metavar='', default=mi, help='Number of main iterations')
	parser.add_argument('-bi', '--breakup_iters', type=int, metavar='', default=bi, help='Number of breakup iterations')
	parser.add_argument('-bg', '--big_thresh', type=int, metavar='', default=bg,
		                help='Number of cells threshold to activate big data mode')
	parser.add_argument('-o', '--out_dir', type=str, metavar='', default=outDir,
		                help='Name of output directory for new tracking files')
	parser.add_argument('-ts', '--time_step', action='store_true', default=ts,
		                help='Toggle file creation for each time step. Default is to combine all times into one file.')
	parser.add_argument('-m', '--map', action='store_true', default=m, help='Toggle map creation')

	args = parser.parse_args()
	return args


def checkArgs(args):
	"""
	Check the user-specified command line arguments for errors not handled by argparse.
	
	Errors will print to console before terminating the script.
	
	Parameters
	----------
	args: Namespace
		Namespace of user-specified arguments returned from getOptions()
				
	"""

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
	bigThreshold = args['big_thresh']

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
		else:
			raise ValueError
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
		else:
		    raise ValueError
	except ValueError:
		print '\nERROR: Invalid end time! Times must be formatted as YYYY-MM-DD-hhmmss, YYYY-MM-DD, YYYY-MM, or YYYY\n'
		sys.exit(2)

	# Everything else
	if '\\' in inSuffix or '/' in inSuffix:
		print '\nERROR: Input directory suffix must not contain / or \\.  Instead got: ' + inSuffix + '\n'
		sys.exit(2)
	else:
		print 'Name of last subdirectory for original tracking files:  ' + inSuffix

	types = ['segmotion', 'probsevere', 'ryan']
	if fType.lower() not in types:
		print 'ERROR: Invalid file type specified. Expected segmotion, probsevere, or ryan.  Instead got: ' + fType + '\n'
	else:
		print 'Data file type: ' + fType

	if bufferDist <= 0 or bufferDist > MAX_BUFFER_DIST:
		print '\nERROR: Buffer distance must be in range (0, ' + str(MAX_BUFFER_DIST) + '].  Instead got: ' + str(
		    bufferDist) + '\n'
		sys.exit(2)
	else:
		print 'Buffer distance between storm cell and Theil-Sen trajectory:  ' + str(bufferDist) + ' km'

	if bufferTime <= 0 or bufferTime > MAX_BUFFER_TIME:
		print '\nERROR: Buffer time must be in range (0, ' + str(MAX_BUFFER_TIME) + '].  Instead got: ' + str(
		    bufferTime) + '\n'
		sys.exit(2)
	else:
		print 'Buffer time for joining two Theil-Sen trajectories:  ' + str(bufferTime) + ' min'

	if joinTime < bufferTime or joinTime > MAX_JOIN_TIME:
		print '\nERROR: Join time must be in range [' + str(bufferTime) + ', ' + str(
		    MAX_JOIN_TIME) + '].  Instead got: ' + str(joinTime) + '\n'
		sys.exit(2)
	else:
		print 'Join time:  ' + str(joinTime) + ' min'

	if joinDist < bufferDist or joinDist > MAX_JOIN_DIST:
		print '\nERROR: Join distance must be in range [' + str(bufferDist) + ', ' + str(
		    MAX_JOIN_DIST) + '].  Instead got: ' + str(joinDist) + '\n'
		sys.exit(2)
	else:
		print 'Join Distance:  ' + str(joinDist) + ' km'

	if minCells < MIN_MIN_CELLS or minCells > MAX_MIN_CELLS:
		print '\nERROR: Min Cells must be in range [' + str(MIN_MIN_CELLS) + ', ' + str(
		    MAX_MIN_CELLS) + '].  Instead got: ' + str(minCells) + '\n'
		sys.exit(2)
	else:
		print 'Minimum number of cells per track:  ' + str(minCells)

	if mainIters < MIN_ITERS or mainIters > MAX_ITERS:
		print '\nERROR: Number of main iterations must be in range [' + str(MIN_ITERS) + ', ' + str(
		    MAX_ITERS) + '].  Instead got: ' + str(mainIters) + '\n'
		sys.exit(2)
	else:
		print 'Number of main iterations:  ' + str(mainIters)

	if breakIters < MIN_BREAKUP_ITERS or breakIters > MAX_BREAKUP_ITERS:
		print '\nERROR: Number of breakup iterations must be in range [' + str(MIN_BREAKUP_ITERS) + ', ' + str(
		    MAX_BREAKUP_ITERS) + '].  Instead got: ' + str(breakIters) + '\n'
		sys.exit(2)
	else:
		print 'Number of breakup iterations:  ' + str(breakIters)

	print 'Big Data Threshold: ' + str(bigThreshold)

	if not os.path.isdir(outDir):
		print 'Unable to locate output directory. The specified location will be created.'
		os.makedirs(outDir)

    # TODO Automatic saving of maps is disabled for now.


#	# Handle map creation variables
#	if mapResults:
#		
#		# Latitude ranges
#		if lats == None or lats == []:
#			lats = [MIN_LAT, MAX_LAT]
#		elif len(lats) % 2 != 0:
#			print '\nInvalid number of latitudes.  There must be an even number of latitudes. Instead ' + str(len(lats)) + ' were given.\n'
#			sys.exit(2)
#		for lat in lats:
#			if lat < MIN_LAT or lat > MAX_LAT:
#				print '\nERROR: Latitude must be in range [' + str(MIN_LAT) + ', ' + str(MAX_LAT) + '].  Instead got: ' + lat + '\n'
#				sys.exit(2)
#			elif lats.index(lat) % 2 == 0 and abs(lat - lats[lats.index(lat) + 1]) <= TOLERANCE:
#				print '\nERROR: Each set of ranges must contain different values.  One row contains the same value twice.\n'
#				sys.exit(2)
#				
#		print '\nThis will create maps for the following latitude ranges: '
#		for i in range (0, len(lats), 2):
#			print str(lats[i]) + '\t' + str(lats[i+1])
#		
#		
#		# Longitude ranges
#		if lons == None or lons == []:
#			lons = [MIN_LON, MAX_LON]
#		elif len(lons) % 2 != 0:
#			print '\nERROR: Invalid number of longitudes.  There must be an even number of longitudes. Instead ' + str(len(lons)) + ' were given.\n'
#			sys.exit(2)
#		elif len(lons) != len(lats):
#			print '\nERROR: The number of longitudes must match the number of latitudes.\n'
#			sys.exit(2)
#		for lon in lons:
#			if lon < MIN_LON or lon > MAX_LON:
#				print '\nERROR: Longitude must be in range [' + str(MIN_LON) + ', ' + str(MAX_LON) + '].  Instead got: ' + lon + '\n'
#				sys.exit(2)
#			elif lons.index(lon) % 2 == 0 and abs(lon - lons[lons.index(lon) + 1]) <= TOLERANCE:
#				print '\nERROR: Each set of ranges must contain different values.  One row contains the same value twice.\n'
#				sys.exit(2)
#				
#		print '\nThis will create maps for the following longitude ranges: '
#		for i in range (0, len(lons), 2):
#			print str(lons[i]) + '\t' + str(lons[i+1])
#			
#		
#		# Map output directory
#		if not os.path.isdir(mapDir):
#			print '\nERROR: ' + mapDir + ' does not exist or is not a directory.'
#			sys.exit(2)
#		else: print '\nMaps will be saved in ' + mapDir
#		
#	else: print 'Mapping disabled'


def find_clusters(stormCells, activeCells):
	"""
	Groups cells with the same track ID into a track dictionary
	
	Parameters
	----------	
	stormCells :  Dictionary 
		Dictionary of storm cells
		
	Returns
	-------
	Dictionary
		Dictionary of storm tracks containing storm cells {ID: [cells]}
	
	"""

	stormTracks = {}
	for cell in activeCells:
		track = stormCells[cell]['track']
		if track in stormTracks:
		    stormTracks[track]['cells'].append(stormCells[cell])
		else:
		    stormTracks[track] = {'cells': [stormCells[cell]]}

	return stormTracks



def theil_sen_single(track):
	"""
	Computes the Theil-Sen fit for a single storm track.
	
	See theil_sen_batch() for more detail.
	
	Parameters
	----------
	track : Dictionary
		The value of a single track within the storm track dictionary
		
	Returns
	-------
	Dictionary
		A storm track dict value with updated items for the provided track
	
	"""

	times = []
	x = []
	y = []
	for cell in track['cells']:
		times.append(time.mktime(cell['time'].timetuple()))  # Converts datetime object to seconds since epoch time
		x.append(cell['x'])
		y.append(cell['y'])

	if len(np.unique(times)) > 1 and len(np.unique(times)) > 1 and len(np.unique(x)) > 1 and len(np.unique(y)) > 1:
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
		track['u'] = 0
		track['v'] = 0
		track['t0'] = datetime.datetime.fromtimestamp(min(times))
		track['tend'] = datetime.datetime.fromtimestamp(max(times))
		track['x0'] = track['cells'][times.index(min(times))]['x']
		track['y0'] = track['cells'][times.index(min(times))]['y']
		track['xf'] = track['cells'][times.index(max(times))]['x']
		track['yf'] = track['cells'][times.index(max(times))]['y']

	return track


def theil_sen_batch(stormTracks):
	"""
	Computes the Theil-Sen fit for each storm track.
	
	Sources: http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mstats.theilslopes.html 
	         https://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
	         
	Parameters
	----------
	stormTracks : Dictionary
		A dictionary of track IDs each containing associated storm cells with Lat, Lon, X, Y, and datetime params
        {'ID':{['x', 'y', 'lat', 'lon', 'times', 'track']}}
        
    Returns
    -------
    Dictionary
    	stormTracks modified with new values 'u', 'v', 't0', 'tend', 'x0', 'y0', 'xf', 'yf'	
    	
	"""
	for track in stormTracks:
		times = []
		x = []
		y = []
		for cell in stormTracks[track]['cells']:
			times.append(time.mktime(cell['time'].timetuple()))  # Converts datetime object to seconds since epoch time
			x.append(cell['x'])
			y.append(cell['y'])

		# print times
		# print cellTracks[track]['x']

		if len(np.unique(times)) > 1 and len(np.unique(x)) > 1 and len(np.unique(y)) > 1:
			try:
				theilSenDataX = stats.theilslopes(x, times)
				theilSenDataY = stats.theilslopes(y, times)
			except ValueError:
				stormTracks[track]['u'] = 0
				stormTracks[track]['v'] = 0
				stormTracks[track]['t0'] = datetime.datetime.fromtimestamp(min(times))
				stormTracks[track]['tend'] = datetime.datetime.fromtimestamp(max(times))
				stormTracks[track]['x0'] = stormTracks[track]['cells'][times.index(min(times))]['x']
				stormTracks[track]['y0'] = stormTracks[track]['cells'][times.index(min(times))]['y']
				stormTracks[track]['xf'] = stormTracks[track]['cells'][times.index(max(times))]['x']
				stormTracks[track]['yf'] = stormTracks[track]['cells'][times.index(max(times))]['y']

			stormTracks[track]['u'] = theilSenDataX[0]
			stormTracks[track]['v'] = theilSenDataY[0]
			stormTracks[track]['t0'] = datetime.datetime.fromtimestamp(min(times))
			stormTracks[track]['tend'] = datetime.datetime.fromtimestamp(max(times))
			stormTracks[track]['x0'] = theilSenDataX[1] + theilSenDataX[0] * (min(times))
			stormTracks[track]['y0'] = theilSenDataY[1] + theilSenDataY[0] * (min(times))
			stormTracks[track]['xf'] = theilSenDataX[1] + theilSenDataX[0] * (max(times))
			stormTracks[track]['yf'] = theilSenDataY[1] + theilSenDataY[0] * (max(times))

		else:
			stormTracks[track]['u'] = 0
			stormTracks[track]['v'] = 0
			stormTracks[track]['t0'] = datetime.datetime.fromtimestamp(min(times))
			stormTracks[track]['tend'] = datetime.datetime.fromtimestamp(max(times))
			stormTracks[track]['x0'] = stormTracks[track]['cells'][times.index(min(times))]['x']
			stormTracks[track]['y0'] = stormTracks[track]['cells'][times.index(min(times))]['y']
			stormTracks[track]['xf'] = stormTracks[track]['cells'][times.index(max(times))]['x']
			stormTracks[track]['yf'] = stormTracks[track]['cells'][times.index(max(times))]['y']

	return stormTracks



def init(l, c):
	"""
	Instantiates global variables for multiprocessing
	
	Paramters
	l : multiprocessing.Lock()
		A multiprocessing.Lock() object to be shared across processes
	c : multiprocessing.Value()
		A multiprocessing.Value('i') object to be shared across processes
		
	"""
	global lock
	global counter
	lock = l
	counter = c


#====================================================================================================================#
#                                                                                                                    #
#  Main Method - Handle user input, read in files, then run calculations                                             #
#                                                                                                                    #
#====================================================================================================================#

if __name__ == '__main__':
	args = vars(getOptions())
	# print args

	# Set Hostname
	hostname = socket.gethostname().split('.')[0]
	print '\n\nSetting hostname to ' + hostname
	print 'Current working directory: ' + os.getcwd()
	print 'Number of processing cores: ' + str(multiprocessing.cpu_count()) + '\n'

	# Check user input.  Type casting is handled by argparse.
	checkArgs(args)

	# If the args check out, save their values here
	startTime = args['start_time']
	endTime = args['end_time']
	inDir = args['input_dir']
	inSuffix = args['dir_suffix']
	fType = args['type'].lower()
	bufferDist = args['buffer_dist']
	bufferTime = timedelta(minutes=int(args['buffer_time']))
	joinTime = timedelta(minutes=int(args['join_time']))
	joinDist = args['join_dist']
	minCells = args['min_cells']
	mainIters = args['main_iters']
	breakIters = args['breakup_iters']
	outDir = args['out_dir']
	outType = args['time_step']
	mapResults = args['map']
	bigThreshold = args['big_thresh']
	bigData = False

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
	
	# TODO:  Everything below here should be in a function (I think)

    #==================================================================================================================#
    #                                                                                                                  #
    #  Read in the files and process data                                                                              #
    #                                                                                                                  #
    #==================================================================================================================#

	# Check for root directory:
	print 'Reading files:'
	if not os.path.isdir(inDir):
		print '\nERROR: Unable to find source directory "' + inDir + '". \nIf using a relative path, please check your working directory.\n'
		sys.exit(2)

	data = readCells.read(fType, inDir, inSuffix, startTime, endTime)
	stormCells = data[0]
	totNumCells = data[1]
	numTrackTimes = data[2]
	dates = sorted(data[3])

	print '\nNumber of files: ' + str(numTrackTimes)
	print 'Total number of storm cells: ' + str(totNumCells)

	if numTrackTimes == 0:
		print 'No valid files found for this time period.  Please check the source directory and specified dates.\n'
		sys.exit(0)

	# Activate bigData mode above a certain threshold (50000 cells)
	if totNumCells >= bigThreshold:
		bigData = True
		print 'Files will be processed in big data mode...'
		if not outType: print 'An output file will be created for each day in the data...'

	# Project onto equidistant coord system
	print '\nProjecting storm cells onto equidistant coordinate system...'

	#Projection variables
	meanLat = np.mean([MIN_LAT, MAX_LAT])
	meanLon = np.mean([MIN_LON, MAX_LON])
	xyDistMax = 0
	llDistMax = 0
	distanceRatio = 0

	# Setup equidistant map projection
	m = Basemap(llcrnrlon=MIN_LON, llcrnrlat=MIN_LAT, urcrnrlon=MAX_LON, urcrnrlat=MAX_LAT,
		        projection='aeqd', lat_0=meanLat, lon_0=meanLon)

	for cell in stormCells:
		stormCells[cell]['x'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[0]
		stormCells[cell]['y'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[1]

	# Find ratio between x-y distances and lat-lon distances
	xMin, yMin = m(MIN_LON, MIN_LAT)
	xMax, yMax = m(MAX_LON, MAX_LAT)

	xyDistMax = np.sqrt((xMin - xMax) ** 2 + (yMin - yMax) ** 2)

	# Find distance between two lat lon coordinates
	# Source: https://en.wikipedia.org/wiki/Great-circle_distance
	# point1 = [MAX_LON, MIN_LAT]
	# point2 = [MIN_LON, MAX_LAT]

	rlat1 = np.radians(MIN_LAT)
	rlat2 = np.radians(MAX_LAT)
	r = 6371  # Mean radius of Earth (km)
	dlon = abs(MAX_LON - MIN_LON)
	dsig = np.arccos(np.sin(rlat1) * np.sin(rlat2) + np.cos(rlat1) * np.cos(rlat2) * np.cos(np.radians(dlon)))
	llDistMax = r * dsig

	distanceRatio = llDistMax / xyDistMax

	print 'Ratio between x-y distances and lat-lon distances: ' + str(distanceRatio)
	print DASHES

    #==================================================================================================================#
    #                                                                                                                  #
    #  Calculations!                                                                                                   #
    #                                                                                                                  #
    #==================================================================================================================#

	print 'Beginning Calculations...'
	REPORT_EVERY = 1000
	if mapResults: scOrigin = copy.deepcopy(stormCells)
	oldCells = []
	total_seconds = datetime.timedelta.total_seconds

	# Run the whole thing for each date
	# Note this will only run once if not bigData (break at end)
	for date in np.unique(dates):
		runstart = datetime.datetime.now()
		dt = datetime.timedelta(hours=6)
		activeCells = []

		# If dealing with a lot of data, only load 1 day (+- 6 hours) at a time
		if bigData:
			print 'Identifying valid cells...'
			for cell in stormCells:
				if stormCells[cell]['time'].date() >= date - dt and stormCells[cell]['time'].date() <= date + dt:
					activeCells.append(cell)
		else:
			activeCells = stormCells.keys()

		if bigData: print 'Processing ' + str(date) + '...'

		# Main iterations
		for i in range(0, mainIters):
			print '\nMain iteration: ' + str(i + 1)

			# Breakup iterations
			for j in range(0, breakIters):

				print '\nBreakup iteration: ' + str(j + 1)
				print 'Finding clusters...'
				stormTracks = find_clusters(stormCells, activeCells)
				print 'Number of clusters: ' + str(len(stormTracks))

				print 'Computing Theil-Sen fit for each cluster...'
				stormTracks = theil_sen_batch(stormTracks)

				# Assign cells to nearest cluster
				print 'Assigning each cell to nearest cluster...'


				# Separated into function for multiprocessing.
				# Note that global variables are not shared between
				# process, so each process returns the modified subset
				# of stormCells to be rejoined later
				def breakupCells(cellSubset):
					
					changedCells = 0

					for cell in cellSubset:
						cellTime = cellSubset[cell]['time']
						cellX = cellSubset[cell]['x']
						cellY = cellSubset[cell]['y']

						# Calculate distances
						minDist = 1e9
						minTrack = stormTracks[min(stormTracks)]
						for track in stormTracks:
							# Only compare to tracks in temporal range
							if not (stormTracks[track]['t0'] - bufferTime <= cellTime <= stormTracks[track]['tend'] + bufferTime):
								continue

							# Preference individual cells to join other tracks
							if len(stormTracks[track]['cells']) < 2 and track == cellSubset[cell]['track']:
								continue

							if stormTracks[track]['u'] == 'NaN':
								xPoint = stormTracks[track]['x0']
								yPoint = stormTracks[track]['y0']
							else:
								xPoint = stormTracks[track]['x0'] + (stormTracks[track]['u'] * (total_seconds(cellTime - stormTracks[track]['t0'])))
								yPoint = stormTracks[track]['y0'] + (stormTracks[track]['v'] * (total_seconds(cellTime - stormTracks[track]['t0'])))

							dist = np.sqrt((cellX - xPoint) ** 2 + (cellY - yPoint) ** 2)
							dist = dist * distanceRatio  # Convert from x,y to km

							# Force cells to be assigned to a track (not NaN)
							# If need be they'll be weeded out in the tie break step later
							if dist < minDist and track != 'NaN':
								minDist = dist
								minTrack = track

						if minDist <= bufferDist:
							if minTrack != cellSubset[cell]['track']: changedCells += 1
							cellSubset[cell]['track'] = minTrack
						else:
							cellSubset[cell]['track'] = 'NaN'

						lock.acquire()
						if counter.value % REPORT_EVERY == 0:
							print '......' + str(counter.value) + ' of ' + str(len(activeCells)) + ' assigned......'
						counter.value += 1
						lock.release()

					return [changedCells, cellSubset]


				# Determine the number of cells per process
				subsets = []
				numPerProc = int(np.ceil(float(len(activeCells)) / multiprocessing.cpu_count()))
				numDone = 0
				for k in xrange(0, len(activeCells), numPerProc):
					temp = {}
					for key in activeCells[k:k + numPerProc]:
						temp[key] = stormCells[key]
					subsets.append(temp)
				del temp

				# Split processing over avilable cores
				l = Lock()
				counter = Value('i', 0)
				with closing(Pool(initializer=init, initargs=(l, counter), processes=20, maxtasksperchild = 1)) as pool:
					results = [pool.apply_async(breakupCells, (subsets[l],)) for l in range(len(subsets))]
					changedCells = sum([result.get()[0] for result in results])
					for result in results:
						for key in result.get()[1]:
							stormCells[key] = result.get()[1][key]

					del results
					del subsets

					pool.close()
					pool.join()
					pool.terminate()

				print 'All cells have been assigned!'
				print 'Number of modified cells: ' + str(changedCells)

				# Stop if there are no changes
				if changedCells == 0: break

            # ------ End of breakup iteration ------ #

			# Find new clusters
			print '\nFinding new clusters after breakup...'
			lastNumTracks = len(stormTracks)
			stormTracks = find_clusters(stormCells, activeCells)
			print 'Original number of clusters: ' + str(lastNumTracks)
			print 'New number of clusters: ' + str(len(stormTracks))

			# Get Theil-Sen fit
			print 'Computing Theil-Sen fit for each new cluster...'
			stormTracks = theil_sen_batch(stormTracks)
			totNumTracks = len(stormTracks)

			# Join similar clusters
			print 'Joining similar clusters...'
			tracks = sorted(stormTracks.keys())
			totNumTracks = len(tracks)
			removeTracks = np.zeros(totNumTracks, dtype = bool)
			merged = 0

			for j in range(0, totNumTracks):
				track1 = tracks[j]

				# Skip tracks with only 1 cell
				if len(stormTracks[track1]['cells']) < 2:
					if j % REPORT_EVERY == 0: print '......' + str(j) + ' of ' + str(totNumTracks) + ' processed for joining......'
					continue
				if track1 == 'NaN': continue

				for k in range(0, j - 1):
					track2 = tracks[k]

					#if len(stormTracks[track2]['cells']) < 2: continue
					if removeTracks[k]: continue
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
					if abs(total_seconds(timeDiff)) > total_seconds(joinTime): continue
			
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
					if velocityDiff > float(bufferDist) / total_seconds(bufferTime): continue
			
					# Check if track predictions are close enough using a subset of 5 cells (or fewer)				
					dist = [None] * len(stormTracks[lateIndex]['cells'][0:6])
					index = 0
					for cell in stormTracks[lateIndex]['cells'][0:6]:
						xActual = cell['x']
						yActual = cell['y']

						cellTime = cell['time']
						xPredict = stormTracks[earlyIndex]['xf'] + (stormTracks[earlyIndex]['u'] * (total_seconds(cellTime - stormTracks[earlyIndex]['tend'])))
						yPredict = stormTracks[earlyIndex]['yf'] + (stormTracks[earlyIndex]['v'] * (total_seconds(cellTime - stormTracks[earlyIndex]['tend'])))

						dist[index] = np.sqrt((xPredict - xActual)**2 + (yPredict - yActual)**2) * distanceRatio
						index += 1
						
					if np.mean(dist) > bufferDist: continue
			
					# If the two tracks survived the process, join them 'cause clearly they're meant to be together ;-)
					removeTracks[k] = True
					merged += 1
					for cell in stormTracks[track2]['cells']:
						cell['track'] = track1
						stormTracks[track1]['cells'].append(cell)

					stormTracks[track1] = theil_sen_single(stormTracks[track1])
		
				if j % REPORT_EVERY == 0:
					print '......' + str(j) + ' of ' + str(totNumTracks) + ' processed for joining......'

			del tracks
			del removeTracks
			
			print 'All tracks have been joined if necessary!'
			print 'Merged ' + str(merged) + ' tracks\n'

			# ------ End of Joining process ------ #

			print 'Finding new clusters after joining...'
			lastNumTracks = len(stormTracks)
			stormTracks = find_clusters(stormCells, activeCells)
			stormTracks = theil_sen_batch(stormTracks)
			totNumTracks = len(stormTracks)
			print 'Original number of clusters: ' + str(lastNumTracks)
			print 'New number of clusters: ' + str(totNumTracks)

			# Break ties (multiple cells assigned to same cluster at same time step)
			print '\nBreaking ties...'
			
			def tieBreak(trackSubset):
				
				breaks = 0
				modifiedCells = {}
				
				for track in trackSubset:
					if len(stormTracks[track]['cells']) < 2:
						lock.acquire()
						if counter.value % REPORT_EVERY == 0: print '......' + str(counter.value) + ' of ' + str(totNumTracks) + ' tracks processed for ties......'
						counter.value += 1
						lock.release()
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
						
						if len(times[thisTime]) > 1:
							cells = times[thisTime]

							# Compare each cell and keep the one closest to the track
							dist = []
							for cell in cells:
								cellX = cell['x']
								cellY = cell['y']
								xPredict = stormTracks[track]['x0'] + (stormTracks[track]['u'] * (total_seconds(thisTime - stormTracks[track]['t0'])))
								yPredict = stormTracks[track]['y0'] + (stormTracks[track]['v'] * (total_seconds(thisTime - stormTracks[track]['t0'])))

								dist.append(np.sqrt((xPredict - cellX) ** 2 + (yPredict - cellY) ** 2) * distanceRatio)

							minCell = cells[dist.index(min(dist))]
							for cell in cells:
								if cell != minCell:
									stormTracks[track]['cells'].remove(cell)
									cell['track'] = 'NaN'
									modifiedCells[stormCells.keys()[stormCells.values().index(cell)]] = cell
								
							breaks += 1

					lock.acquire()
					if counter.value % REPORT_EVERY == 0:
						print '......' + str(counter.value) + ' of ' + str(totNumTracks) + ' tracks processed for ties......'
					counter.value += 1
					lock.release()
				
				return [breaks, modifiedCells]
			
			# Determine the number of cells per process
			subsets = []
			numPerProc = int(np.ceil(float(len(stormTracks.keys())) / multiprocessing.cpu_count()))
			for k in xrange(0, len(stormTracks.keys()), numPerProc):
				temp = {}
				for key in stormTracks.keys()[k:k + numPerProc]:
					temp[key] = stormTracks[key]
				subsets.append(temp)
			del temp

			# Split processing over avilable cores
			l = Lock()
			counter = Value('i', 0)
			with closing(Pool(initializer=init, initargs=(l, counter), processes=20, maxtasksperchild = 1)) as pool:
				results = [pool.apply_async(tieBreak, (subsets[l],)) for l in range(len(subsets))]
				breaks = sum([result.get()[0] for result in results])
				for result in results:
					for key in result.get()[1]:
						stormCells[key] = result.get()[1][key]

				del results
				del subsets

				pool.close()
				pool.join()
				pool.terminate()
			
			print 'All tracks have been processed for tie breaks'
			print 'Number of tie breaks: ' + str(breaks)

        # ------ End of Main iteration ------ #

		print 'Finding new clusters after breakup...'
		lastNumTracks = len(stormTracks)
		stormTracks = find_clusters(stormCells, activeCells)
		stormTracks = theil_sen_batch(stormTracks)
		totNumTracks = len(stormTracks)
		print 'Original number of clusters: ' + str(lastNumTracks)
		print 'New number of clusters: ' + str(totNumTracks)

		# Remove clusters with too few cells
		print '\nRemoving clusters with too few cells...'
		numRemoved = 0
		for track in stormTracks:
			if len(stormTracks[track]['cells']) < int(minCells):
				for cell in stormTracks[track]['cells']:
					cell['track'] = 'NaN'
				numRemoved += 1
		lastNumTracks = len(stormTracks)

		print '\nPerforming final cluster identification...'
		stormTracks = find_clusters(stormCells, activeCells)
		stormTracks = theil_sen_batch(stormTracks)
		stormTracks.pop('NaN', None)

		print 'Number of removed tracks: ' + str(numRemoved + 1)
		print 'Original number of clusters: ' + str(lastNumTracks)
		print 'New number of clusters: ' + str(len(stormTracks))

		print DASHES

		#==================================================================================================================#
		#                                                                                                                  #
		#  Maps!                                                                                                           #
		#                                                                                                                  #
		#==================================================================================================================#

		if mapResults:
			print 'Preparing to plot maps...'

			# Get original storm tracks
			stOrigin = find_clusters(scOrigin, activeCells)
			stOrigin = theil_sen_batch(stOrigin)

			# Handle empty specifications
			lats = [MIN_LAT, MAX_LAT]
			lons = [MIN_LON, MAX_LON]

			# Generate each map
			print 'Plotting figure...'

			fig = plt.figure(1)

			theseLats = lats
			theseLons = lons

			meanLat = np.mean(theseLats)
			meanLon = np.mean(theseLons)

			m = Basemap(llcrnrlon=-119, llcrnrlat=22, urcrnrlon=-64,
						urcrnrlat=49, projection='lcc', lat_1=33, lat_2=45,
						lon_0=-95, resolution='i', area_thresh=10000)

			# Read in shapefiles
			m.readshapefile('counties/c_11au16', name='counties', drawbounds=True, color='#C9CFD1')
			m.readshapefile('States_Shapefiles/s_11au16', name='states', drawbounds=True)
			m.readshapefile('province/province', name='canada', drawbounds=True)

			# Sort cells in each original track by time and then get lat lon pairs for each cell
			for track in stOrigin:
				times = []
				originCellsX = []
				originCellsY = []

				for cell in stOrigin[track]['cells']:
					if bigData and cell['time'].date() != date: continue
					times.append(cell['time'])
				times = sorted(times)
				for cellTime in times:
					for cell in stOrigin[track]['cells']:
						if cell['time'] == cellTime:
							originCellsX.append(m(cell['lon'], cell['lat'])[0])
							originCellsY.append(m(cell['lon'], cell['lat'])[1])
							break

				if len(originCellsX) < 2:
					m.scatter(originCellsX, originCellsY, color='grey', marker='o')
				else:
					m.plot(originCellsX, originCellsY, color='grey', linewidth=BEFORE_WIDTH)

			# Sort cells in each track by time and then get lat lon pairs for each cell
			for track in stormTracks:
				times = []
				finalCellsX = []
				finalCellsY = []
				for cell in stormTracks[track]['cells']:
					if bigData and cell['time'].date() != date: continue
					times.append(cell['time'])

				times = sorted(times)
				for cellTime in times:
					for cell in stormTracks[track]['cells']:
						if cell['time'] == cellTime:
							finalCellsX.append(m(cell['lon'], cell['lat'])[0])
							finalCellsY.append(m(cell['lon'], cell['lat'])[1])
							break

			m.plot(finalCellsX, finalCellsY, color='r', linewidth=AFTER_WIDTH)

			plt.show()

            # Save map to file
            # print 'Saving figure ' + mapDir + '/' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + str((i/2) + 1) + '.png' + '...'
            # plt.savefig(mapDir + '/' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + str((i/2) + 1) + '.png')

			print DASHES

		#==================================================================================================================#
		#                                                                                                                  #
		#  Output                                                                                                          #
		#                                                                                                                  #
		#==================================================================================================================#

		# Reset basemap for conversions
		print 'Preparing output...'
		# Setup equidistant map projection
		m = Basemap(llcrnrlon=MIN_LON, llcrnrlat=MIN_LAT, urcrnrlon=MAX_LON, urcrnrlat=MAX_LAT,
    				projection='aeqd', lat_0=meanLat, lon_0=meanLon)

		# Remove NaN track cells
		print 'Removing unassigned cells...'
		removeCells = []
		for cell in activeCells:
			if stormCells[cell]['track'] == 'NaN': removeCells.append(cell)
		for cell in removeCells: activeCells.pop(activeCells.index(cell))

		print 'Finding new start time, age, and speed for each cell...'

		# Get a smaller dict with the active cell info for efficiency
		# It's easier to use direct access here than to iterate later
		activeStormCells = {}
		for cell in activeCells:
			activeStormCells[cell] = stormCells[cell]

		removeTracks = []
		for track in stormTracks:
			# Remove tracks that aren't part of this day
			if bigData:
				if stormTracks[track]['t0'].date() != date and stormTracks[track]['tend'].date() != date:
					removeTracks.append(track)
					if stormTracks.keys().index(track) % REPORT_EVERY == 0: print '......' + str(stormTracks.keys().index(track)) + ' of ' + str(len(stormTracks)) + ' tracks processed......'
					continue

			# Get start time and age
			# Convert all datetimes to str for JSON
			times = []
			cells = []

			for cell in stormTracks[track]['cells']:
				cell['start_time'] = stormTracks[track]['t0']
				cell['age'] = total_seconds(cell['time'] - cell['start_time'])
				times.append(cell['time'])

			# Sort cells by time
			times = sorted(np.unique(times))
			for cellTime in times:
				for cell in stormTracks[track]['cells']:
					if cell['time'] == cellTime:
						cells.append(cell)
						break

			# Calculate speed and component velocities for each cell
			for cell in cells:
				index = cells.index(cell)
				if index == 0:
					cell['motion_east'] = stormTracks[track]['u'] * distanceRatio * 1000  # m/s
					cell['motion_south'] = -1 * stormTracks[track]['v'] * distanceRatio * 1000  # m/s

				else:
					prevX = cells[index - 1]['x']
					prevY = cells[index - 1]['y']
					prevTime = cells[index - 1]['time']

					cell['motion_east'] = (cell['x'] - prevX) / (total_seconds(cell['time'] - prevTime)) * distanceRatio * 1000  # m/s
					cell['motion_south'] = -1 * (cell['y'] - prevY) / (total_seconds(cell['time'] - prevTime)) * distanceRatio * 1000  # m/s

				cell['speed'] = np.sqrt(cell['motion_east'] ** 2 + cell['motion_south'] ** 2)

			# Cleanup for output
			ids = []
			for cell in stormTracks[track]['cells']:
				# Convert times to strings for JSON
				cell['time'] = str(cell['time'])
				cell['start_time'] = str(cell['start_time'])

				# Remove data specific to this run
				cell.pop('x', None)
				cell.pop('y', None)

				ids.append(activeStormCells.keys()[activeStormCells.values().index(cell)])

			# Only save cell IDs to storm track to save space
			stormTracks[track]['cells'] = ids

			stormTracks[track]['t0'] = str(stormTracks[track]['t0'])
			stormTracks[track]['tend'] = str(stormTracks[track]['tend'])

			# Convert x, y back to lon, lat and km/s to m/s
			stormTracks[track]['lon0'], stormTracks[track]['lat0'] = m(stormTracks[track]['x0'], stormTracks[track]['y0'], inverse=True)
			stormTracks[track]['lonf'], stormTracks[track]['latf'] = m(stormTracks[track]['xf'], stormTracks[track]['yf'], inverse=True)
			stormTracks[track]['u'] = stormTracks[track]['u'] * distanceRatio * 1000  # m/s
			stormTracks[track]['v'] = stormTracks[track]['v'] * distanceRatio * 1000  # m/s

			# Remove data specific to this run
			stormTracks[track].pop('x0', None)
			stormTracks[track].pop('y0', None)
			stormTracks[track].pop('xf', None)
			stormTracks[track].pop('yf', None)

			if stormTracks.keys().index(track) % REPORT_EVERY == 0: print '......' + str(stormTracks.keys().index(track)) + ' of ' + str(totNumTracks) + ' tracks processed......'

		# Remove tracks not part of this date
		if bigData:
			print '\nRemoving ' + str(len(removeTracks)) + ' invalid clusters...'
			for track in removeTracks: 
				stormTracks.pop(track, None)
			print 'New number of clusters: ' + str(len(stormTracks))

		if outType:
			# Print data for each time step

			# Sort cells by time
			times = []
			for cell in activeCells:
				times.append(datetime.datetime.strptime(activeStormCells[cell]['time'], '%Y-%m-%d %H:%M:%S'))

			times = sorted(np.unique(times))

			for cellTime in times:
				cells = {}
				for cell in activeCells:
					if datetime.datetime.strptime(activeStormCells[cell]['time'], '%Y-%m-%d %H:%M:%S') == cellTime:
						cells[cell] = activeStormCells[cell]

				# Print stormCell data for this time step
				filename = (str(cellTime.year) + str(cellTime.month).zfill(2) + str(cellTime.day).zfill(2) + '_' +
	            			str(cellTime.hour).zfill(2) + str(cellTime.minute).zfill(2) + str(cellTime.second).zfill(
		    				2) + '_cells.data')
				print 'Printing ' + filename
				with open(outDir + '/' + filename, 'w') as outfile:
					json.dump(cells, outfile, sort_keys=True, indent=0)

				outfile.close()

		else:
			# Print stormCells to data file
			if bigData:
				filename = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '_cells.data'
				cells = {}
				for track in stormTracks:
					for cell in stormTracks[track]['cells']:
						cells[cell] = activeStormCells[cell]
				print '\nPrinting ' + filename
				with open(outDir + '/' + filename, 'w') as outfile:
					json.dump(cells, outfile, sort_keys=True, indent=0)

			else:
				filename = (str(startTime.year) + str(startTime.month).zfill(2) + str(startTime.day).zfill(2) + '_' +
							str(endTime.year) + str(endTime.month).zfill(2) + str(endTime.day).zfill(2) + '_cells.data')
				print 'Printing ' + filename
				with open(outDir + '/' + filename, 'w') as outfile:
					json.dump(activeStormCells, outfile, sort_keys=True, indent=0)

			outfile.close()

			# Print stormTracks to data file
			if bigData:
				filename = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '_tracks.data'
				print 'Printing ' + filename
				with open(outDir + '/' + filename, 'w') as outfile:
					json.dump(stormTracks, outfile, sort_keys=True, indent=0)
			else:
				filename = (str(startTime.year) + str(startTime.month).zfill(2) + str(startTime.day).zfill(2) + '_' +
			    			str(endTime.year) + str(endTime.month).zfill(2) + str(endTime.day).zfill(2) + '_tracks.data')
				print 'Printing ' + filename
				with open(outDir + '/' + filename, 'w') as outfile:
					json.dump(stormTracks, outfile, sort_keys=True, indent=0)

			outfile.close()

		# Print metadata
		if bigData:
			filename = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '.meta'
		else:
			filename = (str(startTime.year) + str(startTime.month).zfill(2) + str(startTime.day).zfill(2) + '_' +
        				str(endTime.year) + str(endTime.month).zfill(2) + str(endTime.day).zfill(2) + '.meta')
		print 'Printing ' + filename + '\n\n'
		f = open(outDir + '/' + filename, 'w')
		f.write('Run Start: ' + str(runstart) + '\n')
		f.write('Start Time: ' + str(startTime) + '\n')
		f.write('End Time: ' + str(endTime) + '\n')
		f.write('File Type: ' + fType + '\n')
		f.write('Buffer Distance: ' + str(bufferDist) + '\n')
		f.write('Buffer Time: ' + str(bufferTime) + '\n')
		f.write('Join Distance: ' + str(joinDist) + '\n')
		f.write('Join Time: ' + str(joinTime) + '\n')
		f.write('Min Cells per Track: ' + str(minCells) + '\n')
		f.write('Main Iterations: ' + str(mainIters) + '\n')
		f.write('Breakup Iterations: ' + str(breakIters) + '\n')
		f.write('Number of Cells: ' + str(len(activeCells)) + '\n')
		f.write('Completed: ' + str(datetime.datetime.now()))
		f.close()

		# Don't do it again if not bigData
		if not bigData: break

		# Recreate cell values for next iteration
		for cell in activeCells:
			stormCells[cell]['x'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[0]
			stormCells[cell]['y'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[1]
			stormCells[cell]['time'] = datetime.datetime.strptime(str(stormCells[cell]['time']), '%Y-%m-%d %H:%M:%S')

	print '\n\nBest Track has completed succesfully!\n\n\n\n'
