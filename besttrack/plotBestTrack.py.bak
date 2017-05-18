"""
Script to easily read in and plot the results of BestTrack output

Example Usage:  python plotBestTrack.py 20150408_20150409 -i tracks

"""

import argparse
import datetime
from datetime import timedelta
import time
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.stats.mstats as stats
from collections import defaultdict
import json
import os

MIN_LAT = 20
MAX_LAT = 51
MIN_LON = -119
MAX_LON = -62

BEFORE_WIDTH = 4
AFTER_WIDTH = 2

MIN_MIN_CELLS = 2        # Min min number storm cells per track.
MAX_MIN_CELLS = 12       # Max min number storm cells per track.


def getOptions():
	"""
	Retrieve the user-speficified command line arguments
	
	Returns
	--------
	Namespace
		Namespace of parsed arguments returned by ArgumentParser.parse_args()
			
	"""
	
	# Define legal command arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('file_name', type = str, metavar = 'file_name', help = 'File to read in exluding _tracks or _cells, etc')
	parser.add_argument('-i', '--input_dir', type = str, metavar = '', default = 'tracks', help = 'Location of source files')
	parser.add_argument('-s', '--dir_suffix', type = str, metavar = '', default = '', help = 'Name of last subdirectory for source files')
	parser.add_argument('-mc', '--min_cells', type = int, metavar = '', default = 3, help = 'Minimum number of storm cells per track')
	
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
	
	inSuffix = args['dir_suffix']
	minCells = args['min_cells']
	
	if '\\' in inSuffix or '/' in inSuffix:
		print '\nERROR: Input directory suffix must not contain / or \\.  Instead got: ' + inSuffix + '\n'
		sys.exit(2)
	else: print 'Name of last subdirectory for original tracking files:  ' + inSuffix
	
	if minCells < MIN_MIN_CELLS or minCells > MAX_MIN_CELLS:
		print '\nERROR: Min Cells must be in range ['  + str(MIN_MIN_CELLS) + ', ' + str(MAX_MIN_CELLS) + '].  Instead got: ' + str(minCells) + '\n'
		sys.exit(2)
	else: print 'Minimum number of cells per track:  ' + str(minCells)
	
	
#====================================================================================================================#
#                                                                                                                    #
#  Main Method - Handle user input, read in files, then plot 			                                             #
#                                                                                                                    #
#====================================================================================================================#	
if __name__ == '__main__':
	"""Handle user input, read in files, then plot the tracks"""
	
	args = vars(getOptions())
	checkArgs(args)
	
	stormTracks = {}
	stormCells = {}
	
	# Read in track files
	for root, dirs, files in os.walk(args['input_dir']):
		if args['dir_suffix'] != '' and not (files and not dirs and os.path.split(root)[-1] == args['dir_suffix']): continue
		for trackFile in files:
			if trackFile.startswith(args['file_name']) and trackFile.endswith('_tracks.data'):
				# Load tracks
				f = open(root + '/' + trackFile)
				stormTracks = json.load(f)
				f.close()
			
			elif trackFile.startswith(args['file_name']) and trackFile.endswith('_cells.data'):
				# Load cells
				f = open(root + '/' + trackFile)
				stormCells = json.load(f)
				f.close()

	#print stormCells
	#print stormTracks
	
	# Load dimensions
	lats = [MIN_LAT, MAX_LAT]
	lons = [MIN_LON, MAX_LON]


	# Generate each map
	print 'Plotting figure...'

	fig = plt.figure( 1)

	theseLats = lats
	theseLons = lons

	meanLat = np.mean(theseLats)
	meanLon = np.mean(theseLons)

	m = Basemap(llcrnrlon=-119, llcrnrlat=22, urcrnrlon=-64,
					urcrnrlat=49, projection='lcc', lat_1=33, lat_2=45,
					lon_0=-95, resolution='i', area_thresh=10000)

	# Read in shapefiles
	m.readshapefile('counties/c_11au16', name = 'counties', drawbounds = True, color = '#C9CFD1')
	m.readshapefile('States_Shapefiles/s_11au16', name = 'states', drawbounds = True)
	m.readshapefile('province/province', name = 'canada', drawbounds = True)

	# Sort cells in each track by time and then get lat lon pairs for each cell
	for track in stormTracks:
		times = []
		finalCellsX = []
		finalCellsY = []
		for cell in stormTracks[track]['cells']:
			times.append(stormCells[str(cell)]['time'])

		times = sorted(times)
		for cellTime in times:
			for cell in stormTracks[track]['cells']:
				if stormCells[str(cell)]['time'] == cellTime:
					finalCellsX.append(m(stormCells[str(cell)]['lon'], stormCells[str(cell)]['lat'])[0])
					finalCellsY.append(m(stormCells[str(cell)]['lon'], stormCells[str(cell)]['lat'])[1])
					break

		m.plot(finalCellsX, finalCellsY, color = 'r', linewidth = AFTER_WIDTH)			

	plt.show()
        		
        		
        		
        		
        		
	

