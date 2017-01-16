"""
 Concatenates storm tracks from w2segmotionll, probSevere, and post-processed .data (Ryan) files.

 This package is approximately equivalent to w2besttrack with the 
 potential for additional features and greater flexibility.
 This python version was converted from Ryan Lagerquist's MATLAB code
 ryan_best_tracks.m and associated files.

 Author :	David Harrison
 Date   :	July 2016

"""

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
import traceback
from multiprocessing import Pool, Manager, Value, Array, Lock
import multiprocessing
import ctypes
from contextlib import closing
import inspect

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
REPORT_EVERY = 1000

total_seconds = datetime.timedelta.total_seconds


# These functions must be outside the class definition for multiprocessing
def initmultiprocess(l, c):
	"""
	Instantiates global variables for multiprocessing
	
	Paramters
	---------
	l : multiprocessing.Lock()
		A lock object to be shared across processes
	c : multiprocessing.Value()
		A integer value object to be shared across processes
		
	"""
	global lock
	global counter
	lock = l
	counter = c
	
def callBreakup(obj, cellSubset, stormTracks, bufferTime, bufferDist, distanceRatio, totCells):
	"""Calls breakupCells() for multiprocessing
	
	Parameters
	----------
	obj : btengine
		The current btengine object (self)
	Everything else same as btengine.breakupCells()
	
	"""
	
	return obj.breakupCells(cellSubset, stormTracks, bufferTime, bufferDist, distanceRatio, totCells)
	
def callTieBreak(obj, trackSubset, stormTracks, stormCells, totNumTracks, distanceRatio):
	"""Calls tieBreak() for multiprocessing
	
	Parameters
	----------
	obj : btengine
		The current btengine object (self)
	Everything else same as btengine.tieBreak()
	
	"""
	
	return obj.tieBreak(trackSubset, stormTracks, stormCells, totNumTracks, distanceRatio)
	
	
#==================================================================================================================#
#                                                                                                                  #
#  Class definition starts here                                                                                    #
#                                                                                                                  #
#==================================================================================================================#

class btengine:
	"""Class containing all initial parameters and algorithms of the best_track program"""
	
	def __init__(self, stormCells, mainIters = 5, breakIters = 3, bufferDist = 10, bufferTime = 11, joinTime = 16, 
				   joinDist = 50, minCells = 3, dates = [0], startTime = None, endTime = None, mapResults = False, 
				   bigData = False, output = False, fType = '', outDir = '', outType = False):
		"""
		Constructor for the best track engine.  This instantiates all user parameters 
		required for the class methods.  See the README for more info.
	
		Parameters
		----------
		stormCells : Dictionary
			Full dictionary of all stormCells in the dataset
		mainIters : int
			Default 5
			The number of times the whole process is run on the data
		breakIters : int
			Default 3
			The number of times the breakup process is run per main iteration
		bufferDist : int
			Default 10 (km) 
			The distance threshold to use when associated cells with a track
		bufferTime : int
			Default 11 (minutes)
			The time threshold to use when associated cells with a track
		joinTime : int
			Default 16 (minutes)
			The time threshold to use when joining two tracks
		joinDist : int
			Default 50 (km)
			The distance threshold to use when joining two tracks
		minCells : int
			Default 3
			The minimum number of cells required to be in a single track
		dates : List
			List containing all dates (datetime objects) to be processed.
			If not bigData, this can be left at the default value of [0]
		startTime : datetime
			datetime object with the start of the processed time range
			Only required if not in bigData mode
		endTime : datetime
			datetime object with the end of the processed time range
			Only required if not in bigData mode
		mapResults : Bool
			Set True to plot the results of the BestTrack calculations.
			Requires user interaction!
		bigData : Bool
			Set True if handling very large datasets.  Recommended if more than
			50,000 cells are being processed at once.  See README for more info.
		output : Bool
			Set True to generate output files as specified in the README
		fType : String
			The type of the input data: segmotion (.xml),
			probsevere (.ascii), or ryan (.data).  See README for 
			more details
		outDir : String
			Filepath where the files will be saved (can be specified in args)
		outType : Bool
			Will produce a file for every timestep if set to True.  See README
			for more info.
			
		"""	
		
		# Instantiate class variables
		self.stormCells = stormCells
		self.mainIters = mainIters
		self.breakIters = breakIters
		self.bufferDist = bufferDist
		self.bufferTime = bufferTime
		self.joinTime = joinTime
		self.joinDist = joinDist
		self.minCells = minCells
		self.dates = dates
		self.startTime = startTime
		self.endTime = endTime
		self.mapResults = mapResults
		self.bigData = bigData
		self.output = output
		self.fType = fType
		self.outDir = outDir
		self.outType = outType
		
		self.runstart = None
		self.thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
			
	#==================================================================================================================#
	#                                                                                                                  #
	#  Cluster identification                                                                                          #
	#                                                                                                                  #
	#==================================================================================================================#

	def find_clusters(self, stormCells, activeCells):
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

	#==================================================================================================================#
	#                                                                                                                  #
	#  Theil-sen calculations                                                                                          #
	#                                                                                                                  #
	#==================================================================================================================#

	def theil_sen_single(self, track):
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


	def theil_sen_batch(self, stormTracks):
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
		
	#==================================================================================================================#
	#                                                                                                                  #
	#  Breakup Cells                                                                                                   #
	#                                                                                                                  #
	#==================================================================================================================#

	def breakupCells(self, cellSubset, stormTracks, bufferTime, bufferDist, distanceRatio, totCells):
		"""
		Multiprocessing function used to breakup a subset of cells
	
		All cells in the subset are compared to each track and 
		added to the most appropriate one based on distance and time
	
		Parameters
		----------
		cellSubset : Dictionary
			Dictionary of storCells subset to be processed by this worker
		stormTracks : Dictionary
			Full stormTracks dictionary containing information about the current 
			tracks and the cells contained within them
		bufferTime : int
			The time threshold to use when associated cells with a track
		bufferDist : int
			The distance threshold to use when associated cells with a track
		distanceRatio : float
			The ratio between x-y distances and lat-lon distances
		totCells : int
			The total number of active cells being processed across all workers
		
		Returns
		-------
		List
			List containing the number of changed cells and the modified cell subset
		"""
	
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
				print '......' + str(counter.value) + ' of ' + str(totCells) + ' assigned......'
			counter.value += 1
			lock.release()

		return [changedCells, cellSubset]
	
	#==================================================================================================================#
	#                                                                                                                  #
	#  Tie Break	                                                                                                   #
	#                                                                                                                  #
	#==================================================================================================================#
	
	def tieBreak(self, trackSubset, stormTracks, stormCells, totNumTracks, distanceRatio):
		"""
		Mulitprocess function to resolve multiple cells assigned to same cluster at same time step
	
		Parameters
		----------
		trackSubset : Dictionary
			Subset of stormTracks dictionary to be processed by this worker
		stormTracks : Dictionary
			Full dictionary of all stormTracks in the dataset
		stormCells : Dictionary
			Full dictionary of all stormCells in the dataset
		totNumTracks : int
			Number of tracks being processed by all workers
		distanceRatio : float
			The ratio between x-y distances and lat-lon distances
		
		Returns
		-------
		List
			List containing the number of tie-breaks and the modified storm cells
	
		"""
				
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
		
	#==================================================================================================================#
	#                                                                                                                  #
	#  Calculations!                                                                                                   #
	#                                                                                                                  #
	#==================================================================================================================#

	def calculateBestTrack(self):
		"""
		Takes a dictionary of storm cells and merges them into a series of optimal tracks
	
		The data undergoes 3 processes per iteration.  First cells are broken up into first-guess
		track groupings.  These tracks are then joined with each other as necessary, and finally
		any temporal ties between individual cells are resolved.  The process is repeated as many
		times as specified by the user.  See the README for more info.
		
		Returns
		-------
		List
			List containing the modified stormCells and stormTracks dictionaries
	
		"""
		
		# TODO : Go through and eliminate the duplicate variables
	
		stormCells = self.stormCells
		mainIters = self.mainIters
		breakIters = self.breakIters
		bufferDist = self.bufferDist
		bufferTime = self.bufferTime
		joinTime = self.joinTime
		joinDist = self.joinDist
		minCells = self.minCells
		dates = self.dates
		startTime = self.startTime
		endTime = self.endTime
		mapResults = self.mapResults
		bigData = self.bigData
		output = self.output
		outDir = self.outDir
		outType = self.outType
	
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

		# Begin Calculations!

		print 'Beginning Calculations...'
		if mapResults: scOrigin = copy.deepcopy(stormCells)
		oldCells = []

		# Run the whole thing for each date
		# Note this will only run once if not bigData (break at end)
		for date in np.unique(dates):
			self.runstart = datetime.datetime.now()
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
					stormTracks = self.find_clusters(stormCells, activeCells)
					print 'Number of clusters: ' + str(len(stormTracks))

					print 'Computing Theil-Sen fit for each cluster...'
					stormTracks = self.theil_sen_batch(stormTracks)

					# Assign cells to nearest cluster
					print 'Assigning each cell to nearest cluster...'


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
					with closing(Pool(initializer=initmultiprocess, initargs=(l, counter), processes=20, maxtasksperchild = 1)) as pool:
						results = [pool.apply_async(callBreakup, (self, subsets[l], stormTracks, bufferTime, bufferDist, 
													distanceRatio, len(activeCells),)) for l in range(len(subsets))]
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
				stormTracks = self.find_clusters(stormCells, activeCells)
				print 'Original number of clusters: ' + str(lastNumTracks)
				print 'New number of clusters: ' + str(len(stormTracks))

				# Get Theil-Sen fit
				print 'Computing Theil-Sen fit for each new cluster...'
				stormTracks = self.theil_sen_batch(stormTracks)
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

						stormTracks[track1] = self.theil_sen_single(stormTracks[track1])
		
					if j % REPORT_EVERY == 0:
						print '......' + str(j) + ' of ' + str(totNumTracks) + ' processed for joining......'

				del tracks
				del removeTracks
			
				print 'All tracks have been joined if necessary!'
				print 'Merged ' + str(merged) + ' tracks\n'

				# ------ End of Joining process ------ #

				print 'Finding new clusters after joining...'
				lastNumTracks = len(stormTracks)
				stormTracks = self.find_clusters(stormCells, activeCells)
				stormTracks = self.theil_sen_batch(stormTracks)
				totNumTracks = len(stormTracks)
				print 'Original number of clusters: ' + str(lastNumTracks)
				print 'New number of clusters: ' + str(totNumTracks)

				# Break ties (multiple cells assigned to same cluster at same time step)
				print '\nBreaking ties...'
			
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
				with closing(Pool(initializer=initmultiprocess, initargs=(l, counter), processes=20, maxtasksperchild = 1)) as pool:
					results = [pool.apply_async(callTieBreak, (self, subsets[l], stormTracks, stormCells, 
												totNumTracks, distanceRatio,)) for l in range(len(subsets))]
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
			stormTracks = self.find_clusters(stormCells, activeCells)
			stormTracks = self.theil_sen_batch(stormTracks)
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
			stormTracks = self.find_clusters(stormCells, activeCells)
			stormTracks = self.theil_sen_batch(stormTracks)
			stormTracks.pop('NaN', None)

			print 'Number of removed tracks: ' + str(numRemoved + 1)
			print 'Original number of clusters: ' + str(lastNumTracks)
			print 'New number of clusters: ' + str(len(stormTracks))

			print DASHES

			if mapResults:
				self.generateMap(scOrigin, activeCells, stormTracks, bigData, date)
				print DASHES
			
			# Save output
			if output:
				self.generateOutput(activeCells, stormCells, stormTracks, distanceRatio, outDir, 
								startTime, endTime, date, bigData, outType)

	   		# Don't do it again if not bigData
			if not bigData: break
			
			print DASHES
		
		return [stormCells, stormTracks]
		
		
	#==================================================================================================================#
	#                                                                                                                  #
	#  Maps!                                                                                                           #
	#                                                                                                                  #
	#==================================================================================================================#
			
	def generateMap(self, scOrigin, activeCells, stormTracks, bigData, date):
		"""
		Plots a map showing the new tracks compared to the original dataset
	
		Parameters
		----------
		scOrigin : Dictionary
			The original stormCells dictionary (before the calculations)
		activeCells : List
			List containing the cells currently being processed.
			If not in BigData mode, this will be the same as stormCells.keys()
		stormTracks : Dictionary
			Dictionary containing the modified stormTracks (after the calculations)
		bigData : Bool
			Will run in BigData mode if set to True.  See README for more info.
		date : datetime
			The date currently being processed.  Only required if bigData == True
	
		"""
		

		print 'Preparing to plot maps...'

		# Get original storm tracks
		stOrigin = self.find_clusters(scOrigin, activeCells)
		stOrigin = self.theil_sen_batch(stOrigin)

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
		# Source: http://www.nws.noaa.gov/geodata/
		m.readshapefile(self.thisDir + '/mapdata/counties/c_11au16', name='counties', drawbounds=True, color='#C9CFD1')
		m.readshapefile(self.thisDir + '/mapdata/states/s_11au16', name='states', drawbounds=True)
		m.readshapefile(self.thisDir + '/mapdata/provinces/province', name='canada', drawbounds=True)

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
	
		print "Displaying plot.  Please close the figure to continue"	
		plt.show()
	
		# Save map to file
		# print 'Saving figure ' + mapDir + '/' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + str((i/2) + 1) + '.png' + '...'
		# plt.savefig(mapDir + '/' + str(startTime.date()) + '_' + str(endTime.date()) + '_' + str((i/2) + 1) + '.png')
	   
	#==================================================================================================================#
	#                                                                                                                  #
	#  Output                                                                                                          #
	#                                                                                                                  #
	#==================================================================================================================# 
		
	def generateOutput(self, activeCells, stormCells, stormTracks, distanceRatio, outDir, startTime, endTime, date, 
						bigData = False, outType = False) :
		"""
		Generates output files with the results of the BestTrack calculations.
	
		This function produces at least 2 json-encoded files with information about
		stormCells, stormTracks, and meta data.  If outType is True, this will produce
		a stormCells file for every timestep (can be very large!).  See the README for 
		more information.
	
		Parameters
		----------
		activeCells : List
			List containing the cells currently being processed.
			If not in BigData mode, this will be the same as stormCells.keys()
		stormCells : Dictionary
			Full dictionary of all stormCells in the dataset
		stormTracks : Dictionary
			Dictionary containing the modified stormTracks (after the calculations)
		distanceRatio : float
			The ratio between x-y distances and lat-lon distances
		outDir : String
			Filepath where the files will be saved (can be specified in args)
		startTime : datetime
			datetime object with the start of the processed time range
			Only required if not in bigData mode
		endTime : datetime
			datetime object with the end of the processed time range
			Only required if not in bigData mode
		date : 
			datetime object with the date currently being processed
			Only required in bigData mode
		bigData : Bool
			Default False
			Will run in BigData mode if set to True.  See README for more info.
		outType : Bool
			Default False
			Will produce a file for every timestep if set to True.  See README
			for more info		
		
		""" 
		
		totNumTracks = len(stormTracks)
	
		lats = [MIN_LAT, MAX_LAT]
		lons = [MIN_LON, MAX_LON]

		meanLat = np.mean(lats)
		meanLon = np.mean(lons)
	
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

				#print cell['motion_east']
				#print cell['motion_south']
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
		f.write('Run Start: ' + str(self.runstart) + '\n')
		f.write('Start Time: ' + str(self.startTime) + '\n')
		f.write('End Time: ' + str(self.endTime) + '\n')
		f.write('File Type: ' + self.fType + '\n')
		f.write('Buffer Distance: ' + str(self.bufferDist) + '\n')
		f.write('Buffer Time: ' + str(self.bufferTime) + '\n')
		f.write('Join Distance: ' + str(self.joinDist) + '\n')
		f.write('Join Time: ' + str(self.joinTime) + '\n')
		f.write('Min Cells per Track: ' + str(self.minCells) + '\n')
		f.write('Main Iterations: ' + str(self.mainIters) + '\n')
		f.write('Breakup Iterations: ' + str(self.breakIters) + '\n')
		f.write('Number of Cells: ' + str(len(activeCells)) + '\n')
		f.write('Number of Tracks: ' + str(totNumTracks) + '\n')
		f.write('Completed: ' + str(datetime.datetime.now()))
		f.close()

		# Recreate cell values for next iteration
		if bigData:
			for cell in activeCells:
				stormCells[cell]['x'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[0]
				stormCells[cell]['y'] = m(stormCells[cell]['lon'], stormCells[cell]['lat'])[1]
				stormCells[cell]['time'] = datetime.datetime.strptime(str(stormCells[cell]['time']), '%Y-%m-%d %H:%M:%S')

	
