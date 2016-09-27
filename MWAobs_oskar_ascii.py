#!/usr/bin/python
from subprocess import call
from sys import exit
from optparse import OptionParser
from os import environ,getcwd,chdir
from numpy import zeros, pi, sin, cos, real, imag, loadtxt, array, floor, arange
from ephem import Observer
from cmath import exp
from jdcal import gcal2jd
try:
	import pyfits as fits
except ImportError:
	from astropy.io import fits

OSKAR_dir = environ['OSKAR_TOOLS']
R2D = 180.0 / pi
D2R = pi / 180.0
MWA_LAT = -26.7033194444

parser = OptionParser()

parser.add_option('-i', '--ini_file', default=False, help='Enter template oskar .ini - defaults to the template .ini located in $OSKAR_TOOLS/telescopes/--telescope')
parser.add_option('-n','--name', help='Enter prefix name for outputs')
parser.add_option('-s','--srclist', help='Enter location and name of the RTS srclist to use as a sky model')
parser.add_option('-d','--debug',default=False,action='store_true', help='Enable to debug with print statements')
parser.add_option('-o','--data_dir', help='Where to output the finished uvfits')
parser.add_option('-m','--metafits', help='Enter name of metafits file to base obs on')
parser.add_option('-t','--time', help='Enter start,end of sim in seconds from the beginning of the observation (as set by metafits)')
parser.add_option('-x','--twosec', default=False, help='Enable to force a different time cadence - enter the time in seconds')
parser.add_option('-f','--healpix', default=False, help='Enter healpix tag to use base images')
parser.add_option('-a','--telescope', default='MWA_phase1', help='Enter telescope used for simulation. Default = mwa_phase1')
parser.add_option('-b','--band_num', help='Enter band number to simulate')

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
		##Find out full minutes extra and take away
		ext_mins = int(secs / 60.0)
		secs -= 60*ext_mins
		mins += ext_mins
		if mins >= 60.0:
			ext_hours = int(mins / 60.0)
			mins -= 60*ext_hours
			hours += ext_hours
			if hours >= 24.0:
				ext_days = int(hours / 24.0)
				hours -= 24*ext_days
				day += ext_days
			else:
				pass
		else:
			pass
	else:
		pass
	return '%02d-%02d-%d %d:%02d:%05.2f' %(day,month,year,int(hours),int(mins),secs)

def calc_jdcal(date):
	dmy, hms = date.split()
	
	day,month,year = map(int,dmy.split('-'))
	hour,mins,secs = map(float,hms.split(':'))

	##For some reason jdcal gives you the date in two pieces
	##Gives you the time up until midnight of the day
	jd1,jd2 = gcal2jd(year,month,day)
	
	jd3 = (hour + (mins / 60.0) + (secs / 3600.0)) / 24.0
	
	jd = jd1 + jd2 + jd3
	
	##The header of the uvdata file takes the integer, and
	##then the fraction goes into the data array for PTYPE5
	
	return floor(jd), jd - floor(jd)
	

##OSKAR phase tracks, but the MWA correlator doesn't, so we have to undo
##any phase tracking done
def rotate_phase(w,visibility):
	'''Undoes any phase tracking applied to data - to phase track, a phase was applied
	to counter the delay term caused by w term of baseline - so just apply the opposite
	w term. Negative sign was decided through experiment :-S'''

	##I think this has to be negative because OSKAR has negative w when compared
	##to MAPS; works with the extra -w I add in later
	sign = 1
	PhaseConst = 1j * 2 * pi * sign
	
	##theory - so normal phase delay is caused by path difference across
	##a base line, which is u*l + v*m + w*n
	##To phase track, you insert a phase to make sure there is no w contribution at
	##phase centre; this is when n = 1, so you insert a phase thus:
	##a base line, which is u*l + v*m + w*(n - 1)
	##So we just need to remove the effect of the -w term
	PhaseShift = exp( PhaseConst*(w) )
	
	##For XX,XY,YX,YY:
	for k in xrange(4):
		##Take the real and imaginary parts and form a complex number
		complex_data = complex(visibility[k][0],visibility[k][1])
		##Form new complex by multiplying by phase shift?
		new_data = complex_data * PhaseShift
		##Insert new visibility data
		visibility[k][0] = real(new_data)
		visibility[k][1] = imag(new_data)
		
def get_osk_data(oskar_vis_tag=None,polarisation=None):
	OSK = loadtxt('%s_%s.txt' %(oskar_vis_tag,polarisation))

	O_us = OSK[:,1]
	O_vs = OSK[:,2]
	O_ws = OSK[:,3]
	O_res = OSK[:,4]
	O_ims = OSK[:,5]
	
	return O_us,O_vs,O_ws,O_res,O_ims
		
def create_uvfits(template_uvfits=None, freq_cent=None, ra_point=None, dec_point=None, oskar_vis_tag=None, output_uvfits_name=None,date=None):
	
	int_jd, float_jd = calc_jdcal(date)
	
	template_file = fits.open(template_uvfits)
	template_data = template_file[0].data
	antenna_table = template_file[1].data

	# Create uv structure by hand, probably there is a better way of doing this but the uvfits structure is kind of finicky
	n_freq = 1 # only one frequency per uvfits file as read by the RTS

	n_data = len(template_data)

	v_container = zeros((n_data,1,1,1,n_freq,4,3))
	uu = zeros(n_data)
	vv = zeros(n_data)
	ww = zeros(n_data)
	baseline = zeros(n_data)
	date_array = zeros(n_data)
	
	xx_us,xx_vs,xx_ws,xx_res,xx_ims = get_osk_data(oskar_vis_tag=oskar_vis_tag,polarisation='XX')
	yy_us,yy_vs,yy_ws,yy_res,yy_ims = get_osk_data(oskar_vis_tag=oskar_vis_tag,polarisation='YY')
	xy_us,xy_vs,xy_ws,xy_res,xy_ims = get_osk_data(oskar_vis_tag=oskar_vis_tag,polarisation='XY')
	yx_us,yx_vs,yx_ws,yx_res,yx_ims = get_osk_data(oskar_vis_tag=oskar_vis_tag,polarisation='YX')


	for i in range(len(template_data)):
		xx_list = [xx_res[i],xx_ims[i],1.0]
		yy_list = [yy_res[i],yy_ims[i],1.0]
		xy_list = [xy_res[i],xy_ims[i],1.0]
		yx_list = [yx_res[i],yx_ims[i],1.0]
		
		uvdata = [xx_list,yy_list,xy_list,yx_list]
		uvdata = array(uvdata)
		uvdata.shape = (4,3)
		
		v_container[i] = uvdata
		uu[i] = xx_us[i] / freq_cent
		vv[i] = xx_vs[i] / freq_cent
		ww[i] = xx_ws[i] / freq_cent
		baseline[i] = template_data[i][3]
		date_array[i] = float_jd
		rotate_phase(xx_ws[i],v_container[i][0,0,0,0,:,:])

	##UU, VV, WW don't actually get read in by RTS - might be an issue with
	##miriad/wsclean however, as it looks like oskar w = negative maps w
	uvparnames = ['UU','VV','WW','BASELINE','DATE']
	parvals = [uu,vv,ww,baseline,date_array]
		
	uvhdu = fits.GroupData(v_container,parnames=uvparnames,pardata=parvals,bitpix=-32)
	uvhdu = fits.GroupsHDU(uvhdu)

	###Try to copy MAPS as sensibly as possible
	uvhdu.header['CTYPE2'] = 'COMPLEX '
	uvhdu.header['CRVAL2'] = 1.0
	uvhdu.header['CRPIX2'] = 1.0
	uvhdu.header['CDELT2'] = 1.0

	##This means it's linearly polarised
	uvhdu.header['CTYPE3'] = 'STOKES '
	uvhdu.header['CRVAL3'] = -5.0
	uvhdu.header['CRPIX3'] =  1.0
	uvhdu.header['CDELT3'] = -1.0

	uvhdu.header['CTYPE4'] = 'FREQ'
	###Oskar/CASA for some reason adds half of the frequency specified in the 
	###simulation setup. I think this is happens because CASA is unsure
	###what 'channel' the data is - when you run with multiple channels, they
	###are all set to spw = 0, but the output freq is correct. Somethig funky anyway
	###For one channel, set by hand
	uvhdu.header['CRVAL4'] = freq_cent ##(sim freq + half channel width)
	uvhdu.header['CRPIX4'] = template_file[0].header['CRPIX4']
	uvhdu.header['CDELT4'] = template_file[0].header['CDELT4']

	uvhdu.header['CTYPE5'] = template_file[0].header['CTYPE5']
	uvhdu.header['CRVAL5'] = template_file[0].header['CRVAL5']
	uvhdu.header['CRPIX5'] = template_file[0].header['CRPIX5']
	uvhdu.header['CDELT5'] = template_file[0].header['CDELT5']

	uvhdu.header['CTYPE6'] = template_file[0].header['CTYPE6']
	uvhdu.header['CRVAL6'] = template_file[0].header['CRVAL6']
	uvhdu.header['CRPIX6'] = template_file[0].header['CRPIX6']
	uvhdu.header['CDELT6'] = template_file[0].header['CDELT6']

	uvhdu.header['CTYPE7'] = template_file[0].header['CTYPE7']
	uvhdu.header['CRVAL7'] = template_file[0].header['CRVAL7']
	uvhdu.header['CRPIX7'] = template_file[0].header['CRPIX7']
	uvhdu.header['CDELT7'] = template_file[0].header['CDELT7']

	## Write the parameters scaling explictly because they are omitted if default 1/0

	uvhdu.header['PSCAL1'] = 1.0
	uvhdu.header['PZERO1'] = 0.0
	uvhdu.header['PSCAL2'] = 1.0
	uvhdu.header['PZERO2'] = 0.0
	uvhdu.header['PSCAL3'] = 1.0
	uvhdu.header['PZERO3'] = 0.0
	uvhdu.header['PSCAL4'] = 1.0
	uvhdu.header['PZERO4'] = 0.0
	uvhdu.header['PSCAL5'] = 1.0

	uvhdu.header['PZERO5'] = float(int_jd)

	uvhdu.header['OBJECT']  = 'Undefined'                                                           
	uvhdu.header['OBSRA']   = ra_point                                          
	uvhdu.header['OBSDEC']  = dec_point
	
	##ANTENNA TABLE MODS======================================================================

	template_file[1].header['FREQ'] = freq_cent
	
	##MAJICK uses this date to set the LST
	dmy, hms = date.split()
	day,month,year = map(int,dmy.split('-'))
	hour,mins,secs = map(float,hms.split(':'))
	
	rdate = "%d-%02d-%2dT%2d:%2d:%.2f" %(year,month,day,hour,mins,secs)
	
	template_file[1].header['RDATE'] = rdate

	## Create hdulist and write out file
	hdulist = fits.HDUList(hdus=[uvhdu,template_file[1]])
	hdulist.writeto(output_uvfits_name,clobber=True)
	template_file.close()
	hdulist.close()
	
def make_ini(prefix_name,ra,dec,freq,start_time,sky_osm_name,healpix,template_ini):
	out_file = open("%s.ini" %prefix_name,'w+')
	template = open(template_ini).read().split('\n')
	for line in template:
		if "start_frequency_hz" in line:
			line = "start_frequency_hz=%.10f" %freq
		elif "frequency_inc_hz" in line:
			line = "frequency_inc_hz=%.10f" %ch_width
		elif "phase_centre_ra_deg" in line:
			line = "phase_centre_ra_deg=%.10f" %ra
			#line = "phase_centre_ra_deg=356.492904053"
		elif "phase_centre_dec_deg" in line:
			line = "phase_centre_dec_deg=%.10f" %MWA_LAT
		elif "start_time_utc" in line:
			line = "start_time_utc=%s" %start_time
		elif line.split('=')[0]=="length":
			line = "length=%.10f" %dump_time
		elif "oskar_vis_filename" in line:
			line = "oskar_vis_filename=%s.vis" %prefix_name
		elif "channel_bandwidth_hz" in line:
			line = "channel_bandwidth_hz=%.10f" %ch_width
		elif "time_average_sec" in line:
			line = "time_average_sec=%.10f" %dump_time
		elif "ms_filename" in line:
			line = "ms_filename=%s.ms" %prefix_name
		elif "oskar_sky_model" in line:
			line = "oskar_sky_model\\file=%s" %sky_osm_name
		elif "healpix_fits" in line and "file" in line:
			#print 'here1'
			if healpix:
				#print 'here2'
				heal_name = healpix + '_%.3fMHz.fits' %(freq / 1.0e+6)
				line = "healpix_fits\\file=%s" %heal_name
		elif "input_directory" in line:
			line = 'input_directory=%s' %telescope_dir
		else:
			pass
		out_file.write(line+'\n')
	out_file.close()
	

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
##Change in to oskar date format
date,time = intial_date.split('T')
year,month,day = date.split('-')
oskar_date = "%s-%s-%s %s" %(day,month,year,time)

dump_time = float(f[0].header['INTTIME'])

if options.twosec: dump_time = float(options.twosec)

ch_width = float(f[0].header['FINECHAN'])*1e+3
freqcent = float(f[0].header['FREQCENT'])*1e+6
b_width = float(f[0].header['BANDWDTH'])*1e+6
low_freq = freqcent - (b_width/2) - (ch_width/2)

##ephem Observer class, use this to compute LST from the date of the obs 
MRO = Observer()
##Set the observer at Boolardy
MRO.lat, MRO.long, MRO.elevation = '-26:42:11.95', '116:40:14.93', 0
date,time = intial_date.split('T')
MRO.date = '/'.join(date.split('-'))+' '+time
intial_lst = float(MRO.sidereal_time())*R2D
#intial_ra_point = 356.303842095
intial_ra_point = float(MRO.sidereal_time())*R2D
dec_point = MWA_LAT

healpix = options.healpix
telescope_dir = "%s/telescopes/%s" %(OSKAR_dir,options.telescope)
template_uvfits = "%s/telescopes/%s/template_%s.uvfits" %(OSKAR_dir,options.telescope,options.telescope)

if options.ini_file:
	template_ini = options.ini_file
else:
	template_uvfits = "%s/telescopes/%s/template_%s.ini" %(OSKAR_dir,options.telescope,options.telescope)
	
##Unflagged channel numbers
good_chans = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29]
#good_chans = xrange(32)
#good_chans = [0,1,2,3,4]
#good_chans = [0]

##Flagged channel numbers
#bad_chans = [0,1,16,30,31]

band_num = int(options.band_num)
base_freq = ((band_num - 1)*(b_width/24.0)) + low_freq

start_tstep,end_tstep = map(float,options.time.split(','))
tsteps = arange(start_tstep,end_tstep,dump_time)

cwd = getcwd()
tmp_dir = cwd+'/tmp'
data_dir = options.data_dir

outname = options.name

##Sidereal seconds per solar seconds - ie if 1s passes on
##the clock, sky has moved by 1.00274 secs of angle
SOLAR2SIDEREAL = 1.00274
	
###Go to the temporary dir
chdir(tmp_dir)
	
##For frequency channel in band
for chan in good_chans:
	prefixes = []
	ras = []
	freqs = []
	##Take the band base_freq and add on fine channel freq
	freq = base_freq + (chan*ch_width)
	
	##Sky model is the same for all time steps as OSKAR does the horizon clipping itself - 
	##only generate one for all timesteps
	sky_osm_name = "%s_%.3f.osm" %(outname,freq/1e+6)
	##Create the sky model at the obs frequency - this way we can half mimic spectral curvature
	cmd = "python %s/srclist2osm.py -s %s -o %s -f %.10f" %(OSKAR_dir,options.srclist,sky_osm_name,freq)
	run_command(cmd)
	
	##For each time step
	for tstep in tsteps:
		time = add_time(oskar_date,tstep)
		##Precess ra by time since the beginning of the observation 
		##(convert time to angle, change from seconds to degrees)
		##Include half of the time step
		
		ra = intial_ra_point + (((tstep + dump_time/2.0)*SOLAR2SIDEREAL)*(15/3600.0))
		if ra >=360.0: ra -= 360.0
		
		##Make prefix name for the individual time/freq step
		##If less than second time step, RTS needs a different naming convention
		if dump_time < 1:
			prefix_name = "%s_%.3f_%05.2f" %(outname,freq/1e+6,tstep)
		else:
			prefix_name = "%s_%.3f_%02d" %(outname,freq/1e+6,int(tstep))
		
		##Create ini file to run oskar
		make_ini(prefix_name,ra,MWA_LAT,freq,time,sky_osm_name,healpix,template_ini)
		
		##Run the simulation
		cmd = "oskar_sim_interferometer %s.ini" %prefix_name
		run_command(cmd)
		
		##Convert the *.vis into an ascii that we can convert into uvfits
		cmd = "oskar_vis_to_ascii_table -p 4 --baseline_wavelengths -h -v %s.vis %s_XX.txt" %(prefix_name,prefix_name)
		run_command(cmd)
		cmd = "oskar_vis_to_ascii_table -p 5 --baseline_wavelengths -h -v %s.vis %s_XY.txt" %(prefix_name,prefix_name)
		run_command(cmd)
		cmd = "oskar_vis_to_ascii_table -p 6 --baseline_wavelengths -h -v %s.vis %s_YX.txt" %(prefix_name,prefix_name)
		run_command(cmd)
		cmd = "oskar_vis_to_ascii_table -p 7 --baseline_wavelengths -h -v %s.vis %s_YY.txt" %(prefix_name,prefix_name)
		run_command(cmd)
		
		##Clean up the oskar outputs apart from ms set
		cmd = "rm %s.ini %s.vis" %(prefix_name,prefix_name)
		run_command(cmd)

		##Use the centre of the fine channel
		freq_cent = freq + (ch_width / 2.0)
		
		oskar_vis_tag = "%s/%s" %(tmp_dir,prefix_name)
		output_uvfits_name = "%s/%s.uvfits" %(data_dir,prefix_name)
		
		create_uvfits(template_uvfits=template_uvfits, freq_cent=freq_cent, ra_point=ra, dec_point=MWA_LAT, oskar_vis_tag=oskar_vis_tag, output_uvfits_name=output_uvfits_name, date=time)
		
		cmd = "rm %s*txt" %prefix_name
		run_command(cmd)
		
	cmd = "rm %s" %sky_osm_name
	run_command(cmd)
		
chdir(cwd)
