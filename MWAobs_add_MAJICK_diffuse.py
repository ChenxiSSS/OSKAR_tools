#!/usr/bin/python
from subprocess import call
from sys import exit
from optparse import OptionParser
#from os import environ,getcwd,chdir
#from numpy import zeros, pi, sin, cos, real, imag, loadtxt, array, floor, arange

from numpy import pi, arange

#from ephem import Observer
#from cmath import exp
#from jdcal import gcal2jd
try:
	import pyfits as fits
except ImportError:
	from astropy.io import fits

#OSKAR_dir = environ['OSKAR_TOOLS']
R2D = 180.0 / pi
D2R = pi / 180.0
MWA_LAT = -26.7033194444

parser = OptionParser()

parser.add_option('-n','--output_name', help='Enter prefix name for outputs')
parser.add_option('-d','--debug',default=False,action='store_true', help='Enable to debug with print statements')
parser.add_option('-m','--metafits', help='Enter name of metafits file to base obs on')
parser.add_option('-t','--time', help='Enter start,end of sim in seconds from the beginning of the observation (as set by metafits)')
parser.add_option('-x','--twosec', default=False, help='Enable to force a different time cadence - enter the time in seconds')
parser.add_option('-a','--telescope', default='MWA_phase1', help='Enter telescope used for simulation. Default = MWA_phase1')
parser.add_option('-b','--band_num', help='Enter band number to simulate')
parser.add_option('-e', '--base_uvfits', help='Base fits file name and location (e.g. /location/file/uvfits_tag) tag to add diffuse model to')
parser.add_option('-i', '--data_loc', default='./data',	help='Location to output the uvfits to OR location of uvfits if just adding diffuse model. Default = ./data')

options, args = parser.parse_args()
debug = options.debug

def run_command(cmd):
	if debug: print cmd
	call(cmd,shell=True)
	
##Open the metafits file and get the relevant info
try:
	import pyfits
except ImportError:
	import astropy.io.fits as pyfits

try:
	f=pyfits.open(options.metafits)
except Exception,e:
	print 'Unable to open metafits file %s: %s' % (options.metafits,e)
	exit(1)
	
def test_avail(key):
	if not key in f[0].header.keys():
		print 'Cannot find %s in %s' % (key,options.metafits)
		exit(1)

for key in ['DATE-OBS','FREQCENT','FINECHAN','INTTIME','BANDWDTH']:
	test_avail(key)


intial_date = f[0].header['DATE-OBS']
dump_time = float(f[0].header['INTTIME'])

if options.twosec: dump_time = float(options.twosec)

ch_width = float(f[0].header['FINECHAN'])*1e+3
freqcent = float(f[0].header['FREQCENT'])*1e+6
b_width = float(f[0].header['BANDWDTH'])*1e+6
low_freq = freqcent - (b_width/2) - (ch_width/2)

band_num = int(options.band_num)
base_freq = ((band_num - 1)*(b_width/24.0)) + low_freq

start_tstep,end_tstep = map(float,options.time.split(','))
tsteps = arange(start_tstep,end_tstep,dump_time)


good_chans = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29]

for chan in good_chans:
	freq = base_freq + (chan*ch_width)

	sim_command = "time python $MAJICK_DIR/simulate_uvfits.py"
	sim_command += " --freq_start=%.5f" %(freq / 1e+6)
	sim_command += " --num_freqs=1"
	sim_command += " --freq_res=%.5f" %(ch_width / 1e+6)
	sim_command += " --time_start=%.5f " %start_tstep
	sim_command += " --num_times=%d" %len(tsteps)
	sim_command += " --time_res=%.5f" %dump_time
	sim_command += " --date=%s" %intial_date
	sim_command += " --tag_name=%s" %options.output_name
	sim_command += " --base_uvfits=%s" %options.base_uvfits
	sim_command += " --data_loc=%s" %options.data_loc
	sim_command += " --telescope=%s" %options.telescope
	sim_command += " --diffuse"

	run_command(sim_command)


