#=================================================================================#
#                                                                                 #
#                           Best Track README                                     #
#                                                                                 #
#                        Author: David Harrison                                   #
#                         V 1.0 July 26, 2016                                     #
#                                                                                 #
#=================================================================================#


Best Track is approximately equivalent to w2besttrack with the potential for 
additional features and greater flexibility. This python version was converted 
from Ryan Lagerquist's MATLAB code ryan_best_tracks.m and associated files.


REQUIREMENTS:

Best Track was built and tested with Python 2.7.11.  No compatability checks were
performed for Python 3.x or older versions of Python 2.7.x.

The following 3rd party libraries are used in Best Track and must be installed
prior to use:

        * matplotlib
        * numpy
        * scipy
        * BeautifulSoup (bs4)
        * shapely

Additionally, Best Track uses the following Python libraries which should be 
installed with your distribution:

        * argparse
        * socket
        * json
        * datetime
        * collections

If you receive errors when running Best Track, ensure these libraries are
installed and available for use by your distribution.

-----------------------------------------------------------------------------------

TERMINAL USAGE:

usage: best_track.py [-h] [-i] [-s] [-t] [-bd] [-bt] [-jt] [-jd] [-mc] [-mi]
                     [-bi] [-bg] [-o] [-ts] [-m]
                     start_time end_time
                     
        -h    --help               Show the help message and exit.
        
        -i    --input_dir          Specify the directory where the files are stored.
        
        -s    --dir_suffix         Specify the lowest level subdirectory where the files 
                                   are stored.  This is useful if files are spread out
                                   in a hierarchy.  e.g. Use 'scale_0' to get all files
                                   in the hierarchy cref/01/scale_0/, cref/02/scale_0/, etc.
        
        -t    --type               The type of the input data: segmotion (.xml),
                                   probsevere (.ascii), or ryan (.data).  See INPUT for 
                                   more details.
                                                
        -bd   --buffer_dist        The distance threshold to use when associating cells with
                                   a Theil-Sen trajectory (km)
                                                
        -bt   --buffer_time        The time threshold to use when associating cells with
                                   a Theil-Sen trajectory (minutes)
                                                
        -jt   --join_time          The time threshold to use when joining two tracks
                                   (minutes)
                                                
        -jd   --join_dist          The distance threshold to use when joining two tracks
                                   (km)
                                                
        -mc   --min_cells          The minimum number of cells required in a single track
        
        -mi   --main_iters         The number of main iterations the code will run.  More
                                   iterations may increase the likelihood of convergence, 
                                   but will also increase runtime.
                                                
        -bi   --breakup_iters      The number of breakup iterations the code will run in a
                                   single main iteration.  More iterations may increase the
                                   likelihood of convergence, but will also increase runtime.

       	-bg   --big_thresh         The number of storm cells used as a threshold to determine
                                   whether to run the code in "big data mode" (see below)

        -o    --out_dir            Specify where the output files will be saved.
        
        -ts   --time_step          (No arguments) Flag to specify if output files will be 
                                   broken up by time step.  If used, a file will be created
                                   for every time stamp in the data (this can get very large
                                   if used for more than just a day or two at a time).
                                   Default is to create one file for the entire dataset.
                                   See OUTPUT for details.
                                                
        -m    --map                (No arguments) Flag to plot the data on a map.  This will
                                   open a graphic in Python where you can manipulate the 
                                   display and save manually.  Note: this requires user 
                                   interaction!
                                                
        start_time                 A required parameter to specify the beginning of the time
                                   range.  Format must be one of the following:
                                        * YYYY-MM-DD-HHmmss
                                        * YYYY-MM-DD
                                        * YYYY-MM
                                        * YYYY
                                   Example: specifying 2014 would set the start date to
                                   2014-01-01-000000
                                                
        end_time                   A required parameter to specify the end of the time
                                   range.  Format must be one of the following:
                                        * YYYY-MM-DD-HHmmss
                                        * YYYY-MM-DD
                                        * YYYY-MM
                                        * YYYY
                                   Example: specifying 2015 would set the end date to
                                   2015-01-01-000000

Example:
        python best_track.py 2013-05-20 2013-05-21 -i /cref -s scale_0 -t segmotion -m
        
        This would load all .xml files from /cref/.../scale_0 that fall between
        2013-05-20-000000 and 2013-05-21-000000.  After processing the data, it would
        then plot the results and await user input.  Note: the plot must be closed for
        the process to produce the output files.


-----------------------------------------------------------------------------------
        
CONFIG:
        
The default values for Best_track are located in best_track.config.
You can modify these values to meet your needs without having to type in the
arguments at runtime.  Note that any arguments specified at runtime will 
bypass the default arguments for that run.

Ensure that the best_track.config file is located in the same directory as
best_track.py and that the structure of the file (blank lines, etc) remains
unchanged.

-----------------------------------------------------------------------------------

INPUT:

Best Track is able to read in and process 3 different file formats: segmotion,
probsevere, and ryan.  Following is a description of each file type and its
attributes:

        Type:               segmotion
        Extension:          .xml
        Naming Convention:  YYYYMMDD-HHmmss.xml        
        
        Segmotion files are the output .xml files from w2segmotionll.  Cells are
        determined by a specified variable (cref, azshear, etc) and their lat,
        lon, lat radius, lon radius, orientation, and track (row name) are stored.
        Best Track is able to read in segmotion files regardless of what 
        variable is used to determine the cells as long as the aforementioned 
        attributes are available.
        
        Type:               probsevere
        Extension:          .ascii
        Naming Convention:  SSEC_AWIPS_PROBSEVERE_YYYYMMDD_HHmmss.ascii
        
        Probsevere files contain lat, lon, and track (row name) attributes for 
        probsevere objects (cells) in an .ascii file.  Best Track uses the provided
        lat and lon coordinates to calculate the object's centroid (used for the 
        output lat, lon) and the lat radius and lon radius.  Cell orientation is 
        output as NaN for probsevere objects.
        
        Type:               ryan
        Extension:          .data
        Naming Convention:  YYYY-MM-DD-HHmmss_...
        
        Ryan files are processed segmotion files used in the original Best Track
        MATLAB code produced by Ryan Lagerquist.  These files contain 32 lines
        of metadata followed by cell attributes, lon, x, lat, y repeating every
        5 lines for each cell.  Best Track continues to support ryan processed
        files to maintain backwards compatability with the original MATLAB code.
        
Note that the naming conventions for each file type MUST be followed for Best
Track to successfully read in the data!

Big Data Mode:

After reading in all files in the input directory in the specified time range,
if there are more storm cells than the specified threshold (default 50000), Best
Track will run in "Big Data Mode".  In this mode, each day will be processed 
individually and output files (See Output) will be created for each valid date 
within the specified time range.  This is intended to decrease total runtime
on very large datasets, but will not necessarily reduce the memory usage of the script.
CAUTION: Reading in a significant number of files or files with a significant
number of storm cells may result in poor performance and/or out of memory errors.

-----------------------------------------------------------------------------------

OUTPUT:

Single File Method:

If the -ts flag is not set at runtime or in the CONFIG file, Best Track will
default to single file output.  Contrary to the name, this will create a total
of 3 files in the output directory: ..._cells.data, ..._tracks.data, and ...meta.
Here we will describe the contents of each file:

        File Name:   YYYYMMDD_YYYYMMDD_cells.data
        Attributes:  cellID               Unique identifier for the cell
                     age                  Age of the cell in seconds
                     lat                  Centroid latitude
                     lon                  Centroid longitude
                     latr                 Latitude radius
                     lonr                 Longitude radius
                     motion_east          Eastward speed in m/s
                     motion_south         Southward speed in m/s
                     orientation          Object orientation in degrees
                     speed                Cell speed in m/s
                     start_time           Start time of associated track
                     time                 Cell valid time
                     track                Associated storm track
        Structure:
                     {
                     cellID:{
                     age:,
                     lat:,
                     latr:,
                     lon:,
                     lonr:,
                     motion_east:,
                     motion_south:,
                     orientation:,
                     speed:,
                     start_time:,
                     time:,
                     track: }
                     ...
                     }
                                
                *        *        *        *        *        *
                                
        File Name:   YYYYMMDD_YYYYMMDD_tracks.data
        Attributes:  trackID          Unique identifier for the track
                     cells            Array of cell IDs associated with the track
                     lat0             Starting lat of the track
                     latf             Ending lat of the track
                     lon0             Starting lon of the track
                     lonf             Ending lon of the track
                     t0               Starting time of the track
                     tend             Ending time of the track
                     u                Mean zonal velocity in m/s
                     v                Mean meridonal velocity in m/s
                                
        Structure:
                     {
                     trackID:{
                     cells:[cellIDs],
                     lat0:,
                     latf:,
                     lon0:,
                     lonf:,
                     t0:,
                     tend:,
                     u:,
                     v: }
                     ...
                     }
                                
                *        *        *        *        *        *
                                
        File Name:   YYYYMMDD_YYYYMMDD.meta
        Attributes:  Start Time
                     End Time
                     File Type
                     Buffer Distance
                     Buffer Time
                     Join Distance
                     Join Time
                     Min Cells
                     Min Iterations
                     Breakup Iterations
                     Number of Cells      Total number of cells processed
                     Completed            Local time the script completed
                                
The file name YYYYMMDD_YYYYMMDD... uses the start time and end time for which
the script was run.  CAUTION: If Best Track is run for multiple time ranges within
the same day, the output files will have the same name and will overwrite the
existing files.

Big Data Method:

If the -ts flag is not used and Best Track is run in "Big Data Mode" 
(see Input), a YYYYMMDD_cells.data, YYYYMMDD_tracks.data, and 
YYYYMMDD.meta file will be produced for each day in the specified time frame for
which data exists.  These files will have the same structure as their "normal
mode" single file method counterparts (see above).

Time Step Method:

If the -ts flag is used at runtime or in the CONFIG file, Best Track will output
one file for every time step in the date range.  For instance, running Best Track
for 1 hour on a dataset with 2 minute time steps will produce 30 files plus a
.meta file.  This is useful if you want to be able to load a subset of processed
data later, but can produce a large quantity of files and will require more 
post-processing by the user later on.  Each file is equivalent to the 
..._cells.data file described for the previous method, but only cells valid at that
time step will be included.  No ..._tracks.data file is made using this method.

        File Name:   YYYYMMDD_HHmmss_cells.data
        Attributes:  See ..._cells.data in Single File Method
        
-----------------------------------------------------------------------------------

READING THE DATA:

All Best Track output files are encoded using the json Python library.  This makes
the files easy to read by both humans and Python scripts.  In order to read in data
from a file, simply do the following (in Python):

        import json
        
        filepath = 'tracks/20130520_20130521_cells.data'
        f = open(filepath)
        cells = json.load(f)
        f.close()
        
This will create a dictionary of storm cells with the same structure specified in 
the file (see OUTPUT).  Note that all keys in the dictionary will be represented
as strings.  Access the individual attributes like:

        speed = cell['0']['speed']
        
This will work for tracks as well.  Note: the .meta file is not json encoded, but
can be easily parsed with a standard file reader if necessary.
