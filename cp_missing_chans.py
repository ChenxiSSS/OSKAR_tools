from subprocess import call
from sys import exit
from optparse import OptionParser
from os import environ,getcwd,chdir
from numpy import arange
from glob import glob
try:
	import pyfits as fits
except ImportError:
	from astropy.io import fits

parser = OptionParser()

parser.add_option('-n','--output_name', help='Enter prefix name for outputs')
parser.add_option('-d','--debug',default=False,action='store_true', help='Enable to debug with print statements')
parser.add_option('-o','--data_dir', help='Where to output the finished uvfits')
parser.add_option('-m','--metafits', help='Enter name of metafits file to base obs on')
parser.add_option('-t','--time', help='Enter start,end of sim in seconds from the beginning of the observation (as set by metafits)')
parser.add_option('-x','--twosec', default=False, help='Enable to force a different time cadence - enter the time in seconds')
parser.add_option('-a','--telescope', default='MWA_phase1', help='Enter telescope used for simulation. Default = MWA_phase1')
parser.add_option('-b','--band_num', help='Enter band number to fill missing channels for')

options, args = parser.parse_args()
debug = options.debug


def run_command(cmd):
	if debug: print cmd
	call(cmd,shell=True)
	
def add_time(date_time,time_step):
	'''Take the time string format that oskar uses ('23-08-2013 17:54:32.0'), and add a time time_step (seconds).
	Return in the same format - NO SUPPORT FOR CHANGES MONTHS CURRENTLY!!'''
	date,time = date_time.split()
	day,month,year = map(int,date.split('-'))
	hours,mins,secs = map(float,time.split(':'))
	##Add time
	secs += time_step
	if secs >= 60.0:
		secs -= 60.0
		mins += 1.0
		if mins >= 60.0:
			mins -= 60.0
			hours +=1.0
			if hours >= 24.0:
				hours -= 24.0
				day += 1
			else:
				pass
		else:
			pass
	else:
		pass
	return '%02d-%02d-%d %d:%02d:%05.2f' %(day,month,year,int(hours),int(mins),secs)

def edit_freq(base_uvfits, new_uvfits, base_freq,ch_width):
	'''Edits the freq of the file, just incase RTS is doing something mental'''
	oskar_file = fits.open(base_uvfits)
	oskar_file[0].header['CRVAL4'] = base_freq + (ch_width / 2.0)
	oskar_file[1].header['FREQ'] = base_freq + (ch_width / 2.0)
	oskar_file.writeto(new_uvfits,clobber=True)
	oskar_file.close()

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

for key in ['DATE-OBS','RA','DEC','FREQCENT','FINECHAN','INTTIME','BANDWDTH']:
	test_avail(key)


obs_date = f[0].header['DATE-OBS']
##Change in to oskar date format
date,time = obs_date.split('T')
year,month,day = date.split('-')
oskar_date = "%s-%s-%s %s" %(day,month,year,time)

dump_time = float(f[0].header['INTTIME'])

if options.twosec: dump_time = float(options.twosec)

ch_width = float(f[0].header['FINECHAN'])*1e+3
freqcent = float(f[0].header['FREQCENT'])*1e+6
b_width = float(f[0].header['BANDWDTH'])*1e+6
low_freq = freqcent - (b_width/2) - (ch_width/2)
ra_point = float(f[0].header['RA'])
dec_point = float(f[0].header['DEC'])


##Sidereal seconds per solar seconds - ie if 1s passes on
##the clock, sky has moved by 1.00274 secs of angle
SOLAR2SIDEREAL = 1.00274

##Unflagged channel numbers
good_chans = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29]
#good_chans = [2,3]

##Flagged channel numbers
bad_chans = [0,1,16,30,31]
	
band_num = int(options.band_num)
base_freq = ((band_num - 1)*(b_width/24.0)) + low_freq

start_tstep,end_tstep = map(float,options.time.split(','))
tsteps = arange(start_tstep,end_tstep,dump_time)

cwd = getcwd()
tmp_dir = cwd+'/tmp'
#cmd = "mkdir %s" %tmp_dir
#run_command(cmd)
data_dir = options.data_dir

##Even though we flag 5 channels, RTS still needs to read them in,
##so just copy the base channel
chdir(data_dir)
for chan in bad_chans:
		##Take the band base_freq and add on fine channel freq
		freq = base_freq + (chan*ch_width)
		##For each time step
		for tstep in tsteps:
			##If less than second time step, RTS needs a different naming convention
			if dump_time < 1:
				new_uvfits = "%s_%.3f_%05.2f.uvfits" %(options.output_name,freq/1e+6,tstep)
				base_uvfits = "%s_%.3f_%05.2f.uvfits" %(options.output_name,(base_freq + (good_chans[0]*ch_width))/1e+6,tstep)
			else:
				new_uvfits = "%s_%.3f_%02d.uvfits" %(options.output_name,freq/1e+6,int(tstep))
				base_uvfits = "%s_%.3f_%02d.uvfits" %(options.output_name,(base_freq + (good_chans[0]*ch_width))/1e+6,int(tstep))
			#cmd = "cp %s %s" %(first_good_chan_name,uvfits_name)
			#run_command(cmd)
			##After creating file, change frequency in header of both data and antenna files
			edit_freq(base_uvfits, new_uvfits, base_freq, ch_width)
			
