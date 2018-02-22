#!/usr/bin/env python
from subprocess import call
from sys import exit
from optparse import OptionParser
from numpy import zeros, pi, sin, cos, real, imag, loadtxt, array, floor, arange, ones, where, savetxt
from numpy import exp as n_exp
from ephem import Observer
from cmath import exp
from sys import path as sys_path
try:
    from jdcal import gcal2jd
except ImportError:
    sys_path.append('/lustre/projects/p048_astro/jline/software/jdcal-1.3')
    from jdcal import gcal2jd
from os import environ,getcwd,chdir,makedirs,path
from struct import unpack
try:
    OSKAR_dir = environ['OSKAR_TOOLS']
except:
    OSKAR_dir = '/lustre/projects/p048_astro/jline/software/OSKAR_tools'
sys_path.append(OSKAR_dir)
from MWAobs_oskar_ascii_lib import *
from astropy.io import fits

R2D = 180.0 / pi
D2R = pi / 180.0
MWA_LAT = -26.7033194444
VELC = 299792458.0

parser = OptionParser()


parser.add_option('-a','--telescope', default='%s/telescopes/MWA_phase1' %OSKAR_dir, help='Enter telescope used for simulation. Default = $OSKAR_TOOLS/telescopes/MWA_phase1')
parser.add_option('-b','--band_num', help='Enter band number to simulate')
parser.add_option('-c','--osm', default=False, help='Location of OSKAR osm sky model to use')
parser.add_option('-d','--debug',default=False,action='store_true', help='Enable to debug with print statements')
parser.add_option('-e','--antenna_coord_file',default='%s/telescopes/MWA_phase1/MWATools-antenna_locations.txt' %OSKAR_dir,
                  help='If creating a telescope model with dipole flags, use this as the array layout. Defaults to MWA_phase1 ($OSKAR_TOOLS/telescopes/MWA_phase1/MWATools-antenna_locations.txt)')
parser.add_option('-f','--healpix', default=False, help='Enter healpix tag to use base images NOT CURRENTLY USED')
parser.add_option('-g','--fit_osm', default=False, help='Location of sky parameters to create osm from')
parser.add_option('-i', '--ini_file', default=False, help='Enter template oskar .ini - defaults to the template .ini located in $OSKAR_TOOLS/telescopes/--telescope')
parser.add_option('-j','--flag_dipoles',default=False,action='store_true', help='Add to switch on dipole flagging via the metafits file. NOTE needs a metafits that has the correct data within')

parser.add_option('-o','--data_dir', help='Where to output the finished uvfits - default is ./data',default=False)
parser.add_option('-m','--metafits', help='Enter name of metafits file to base obs on')
parser.add_option('-n','--output_name', help='Enter prefix name for outputs')
parser.add_option('-s','--srclist', default=False, help='Enter location and name of the RTS srclist to use as a sky model')
parser.add_option('-t','--time', help='Enter start,end of sim in seconds from the beginning of the observation (as set by metafits)')
parser.add_option('-x','--time_int', default=False, help='Enable to force a different time cadence from that in the metafits - enter the time in seconds')


options, args = parser.parse_args()
debug = options.debug

def run_command(cmd):
    if debug: print cmd
    call(cmd,shell=True)

try:
    f=fits.open(options.metafits)
except Exception,e:
    print 'Unable to open metafits file %s: %s' % (options.metafits,e)
    exit(1)
    
def test_avail(key):
    if not key in f[0].header.keys():
        print 'Cannot find %s in %s' % (key,options.metafits)
        exit(1)

for key in ['DATE-OBS','FREQCENT','FINECHAN','INTTIME','BANDWDTH','AZIMUTH','ALTITUDE']:
    test_avail(key)


intial_date = f[0].header['DATE-OBS']
##Change in to oskar date format
date,time = intial_date.split('T')
year,month,day = date.split('-')
oskar_date = "%s-%s-%s %s" %(day,month,year,time)

time_int = float(f[0].header['INTTIME'])

if options.time_int: time_int = float(options.time_int)

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
initial_ra_point = float(MRO.sidereal_time())*R2D
dec_point = MWA_LAT

healpix = options.healpix
telescope_dir = options.telescope
telescope_name = options.telescope.split('/')[-1]
template_uvfits = fits.open("%s/template_%s.uvfits" %(telescope_dir,telescope_name))
template_data = template_uvfits[0].data
num_baselines = len(template_data)
num_freq_channels = 32

if options.ini_file:
    template_ini = options.ini_file
else:
    template_ini = "%s/template_%s.ini" %(telescope_dir,telescope_name)
template_ini = open(template_ini).read().split('\n')
    
##Unflagged channel numbers
good_chans = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29]
central_freq_chan = 15
#good_chans = xrange(32)
good_chans = [2,3]
central_freq_chan = 2

##Flagged channel numbers
#bad_chans = [0,1,16,30,31]

band_num = int(options.band_num)
base_freq = ((band_num - 1)*(b_width/24.0)) + low_freq

start_tstep,end_tstep = map(float,options.time.split(','))
tsteps = arange(start_tstep,end_tstep,time_int)

##Total length of data axis
n_data = num_baselines * len(tsteps)

cwd = getcwd()
tmp_dir = cwd+'/tmp'
if not path.exists(tmp_dir):
    makedirs(tmp_dir)

if options.data_dir:
    data_dir = options.data_dir
else:
    data_dir = cwd + '/data'
##Check to see if data directory exists; if not create it
if not path.exists(data_dir):
    makedirs(data_dir)

outname = options.output_name

##Sidereal seconds per solar seconds - ie if 1s passes on
##the clock, sky has moved by 1.00274 secs of angle
SOLAR2SIDEREAL = 1.00274
    
int_jd, float_jd = calc_jdcal(oskar_date)
#print int_jd, float_jd

##Need an array the length of number of baselines worth of the fractional jd date
float_jd_array = ones(num_baselines)*float_jd

##Create empty data structures for final uvfits file
v_container = zeros((n_data,1,1,1,num_freq_channels,4,3))
uus = zeros(n_data)
vvs = zeros(n_data)
wws = zeros(n_data)
baselines_array = zeros(n_data)
date_array = zeros(n_data)

###Go to the temporary dir
chdir(tmp_dir)

##Depending on which type of foreground model is required,
##generate or declare the osm
if options.osm:
    sky_osm_name = options.osm
elif options.srclist:
    ##For frequency channel in band
    for chan in good_chans:
        freq = base_freq + (chan*ch_width)
        ##Sky model is the same for all time steps as OSKAR does the horizon clipping itself - 
        ##only generate one for all timesteps
        sky_osm_name = "%s_%.3f.osm" %(outname,freq/1e+6)
        ##Create the sky model at the obs frequency - this way we can half mimic spectral curvature
        cmd = "python %s/srclist2osm.py -s %s -o %s -f %.10f" %(OSKAR_dir,options.srclist,sky_osm_name,freq)
        run_command(cmd)
elif options.fit_osm:
    ##Fill in with new feature when ready
    pass
else:
    print("No valid sky model option declared; need either --osm, --srclist, --fit_osm")
    print("Exiting now")
    exit(0)
    
###Create a copy of the telescope model===============================
##Add in a permitted_beams.txt file to stop the beam phase tracking
##If requested, flag dipoles via metafits file TODO

azimuth = float(f[0].header['AZIMUTH'])
altitude = float(f[0].header['ALTITUDE'])


##Annoying things happen if we use an old temporary telescope model - usually
##hangs around after an incomplete run. Kill it with fire if exists
if path.exists('%s/telescope_%s_band%02d' %(tmp_dir,outname,band_num)):
    cmd = 'rm -rf %s/telescope_%s_band%02d' %(tmp_dir,outname,band_num)
    run_command(cmd)
else:
    pass

if options.flag_dipoles:
    create_flagged_telescope(metafits=options.metafits,antenna_coord_file=options.antenna_coord_file,azimuth=azimuth,altitude=altitude,
                                telescope_dir='%s/telescope_%s_band%02d' %(tmp_dir,outname,band_num))
else:
    ##Just copy the template telescope
    cmd = 'cp -r %s %s/telescope_%s_band%02d' %(telescope_dir,tmp_dir,outname,band_num)
    run_command(cmd)
    #OSKAR naturally beam forms towards phase centre
    #Read the pointing from the metafits and force observation towards single pointing
    
    permitted = open('%s/station/permitted_beams.txt' %(telescope_dir),'w+')
    permitted.write('%.5f %.5f' %(azimuth,altitude))
    permitted.close

    permitted = open('%s/permitted_beams.txt' %(telescope_dir),'w+')
    permitted.write('%.5f %.5f' %(azimuth,altitude))
    permitted.close



##Reset the telescope_dir to copied location
telescope_dir = '%s/telescope_%s_band%02d' %(tmp_dir,outname,band_num)

#@profile
def the_main_loop(tsteps=None):
    if options.osm:
        
        ##Make prefix name for the individual time step
        ##If less than second time step, RTS needs a different naming convention
        #prefix_name = "%s_band%02d_%.3f" %(outname,band_num,tstep)
        prefix_name = "%s_band%02d_t%.1f-%.1f" %(outname,band_num,tsteps[0],tsteps[-1])
        
        ##Start at the first good_chan freq, and do enough chans to cover first good chan to last good chan
        num_channels = len(good_chans)
        oskar_channels = range(good_chans[0],good_chans[-1]+1)
        
        #print oskar_channels
        
        oskar_inds = [oskar_channels.index(good_chan) for good_chan in good_chans]
        num_oskar_channels = len(oskar_channels)
        #print "NUM CHANNELS",num_channels
        osk_low_freq = base_freq + (good_chans[0]*ch_width)
        #osk_low_freq = base_freq
        
        num_time_steps = len(tsteps)
        obs_time_length = num_time_steps * time_int
        
        #oskar_date
        
        make_ini(prefix_name=prefix_name,ra=initial_ra_point,dec=MWA_LAT,freq=osk_low_freq,start_time=oskar_date,
                    sky_osm_name=options.osm,healpix=healpix,num_channels=num_oskar_channels,
                    template_ini=template_ini,ch_width=ch_width,time_int=time_int,
                    telescope_dir=telescope_dir,num_time_steps=num_time_steps,obs_time_length=obs_time_length)
        
        ##Run the simulation
        cmd = "oskar_sim_interferometer --quiet %s.ini" %prefix_name
        run_command(cmd)
        
        ##Read in the data directly from the binary file
        ##This file contains all frequency channels for this time step
        #num_vis = num_baselines * num_oskar_channels * num_time_steps
        num_vis = num_oskar_channels * num_time_steps
        
        uu,vv,ww,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims = read_oskar_binary(filename="%s.vis" %prefix_name,num_time_steps=num_time_steps,num_channels=num_oskar_channels,num_baselines=num_baselines)

        ##Useful for checking what is actually coming out of OSKAR
        #cmd = "oskar_vis_to_ascii_table -p 4 %s.vis" %prefix_name
        #run_command(cmd)
        
        ##Clean up the oskar outputs
        cmd = "rm %s.ini %s.vis" %(prefix_name,prefix_name)
        run_command(cmd)
        
        ###For each time step
        for time_ind,tstep in enumerate(tsteps):
            ##Precess ra by time since the beginning of the observation 
            ##(convert time to angle, change from seconds to degrees)
            ##Include half of the time step
            time = add_time(oskar_date,tstep)
            
            ##Need to know where in the uvfits file structure to put data
            array_time_loc = num_baselines*time_ind

            ##If we are using a sky model with fixed spectral indexes, can run all
            ##frequencies using the same osm model, and so only need to run OSKAR once
            ##per time step
        
            for oskar_ind,chan in zip(oskar_inds,good_chans):
                freq = base_freq + (chan*ch_width)
                
                time_ind_lower = time_ind*num_baselines
                
                osk_low = (num_oskar_channels*time_ind_lower) + oskar_ind * num_baselines
                osk_high = (num_oskar_channels*time_ind_lower) + (oskar_ind + 1) * num_baselines
                ##Stored in metres in binary, convert to wavelengths
                chan_uu = uu[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
                chan_vv = vv[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
                chan_ww = ww[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
                
                ##Select the correct visibilities from the binary data
                chan_xx_res = xx_res[osk_low:osk_high]
                chan_xx_ims = xx_ims[osk_low:osk_high]
                chan_xy_res = xy_res[osk_low:osk_high]
                chan_xy_ims = xy_ims[osk_low:osk_high]
                chan_yx_res = yx_res[osk_low:osk_high]
                chan_yx_ims = yx_ims[osk_low:osk_high]
                chan_yy_res = yy_res[osk_low:osk_high]
                chan_yy_ims = yy_ims[osk_low:osk_high]
                
                ##Make complex numpy arrays
                comp_xx = make_complex(chan_xx_res,chan_xx_ims)
                comp_xy = make_complex(chan_xy_res,chan_xy_ims)
                comp_yx = make_complex(chan_yx_res,chan_yx_ims)
                comp_yy = make_complex(chan_yy_res,chan_yy_ims)
                
                ##Remove the phase tracking added in by OSKAR
                rotated_xx = rotate_phase(wws=chan_ww,visibilities=comp_xx)
                rotated_xy = rotate_phase(wws=chan_ww,visibilities=comp_xy)
                rotated_yx = rotate_phase(wws=chan_ww,visibilities=comp_yx)
                rotated_yy = rotate_phase(wws=chan_ww,visibilities=comp_yy)
                
                #print oskar_ind,chan_ww[112],freq,comp_xx[112],rotated_xx[112]

                ##Use the centre of the fine channel
                freq_cent = freq + (ch_width / 2.0)
                
                ##Populate the v_container at the correct spots
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,0,0] = real(rotated_xx)
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,0,1] = imag(rotated_xx)
                
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,1,0] = real(rotated_yy)
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,1,1] = imag(rotated_yy)
                
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,2,0] = real(rotated_xy)
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,2,1] = imag(rotated_xy)
                
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,3,0] = real(rotated_yx)
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,3,1] = imag(rotated_yx)
                
                ##Set the weights of everything to ones
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,0,2] = ones(len(chan_uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,1,2] = ones(len(chan_uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,2,2] = ones(len(chan_uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,3,2] = ones(len(chan_uu))
                
                ##Only set the u,v,w to the central frequency channel
                ##Header of uvfits has to match central_freq_chan + 1 (uvfits 1 ordered, python zero ordered)
                if chan == central_freq_chan:
                    central_freq_chan_value = freq_cent
                    ##u,v,w stored in seconds by uvfits files
                    chan_uu /= freq_cent
                    chan_vv /= freq_cent
                    chan_ww /= freq_cent
                    
                    uus[array_time_loc:array_time_loc+num_baselines] = chan_uu
                    vvs[array_time_loc:array_time_loc+num_baselines] = chan_vv
                    wws[array_time_loc:array_time_loc+num_baselines] = chan_ww
                    
                    ##Fill in the baselines using the first time and freq uvfits
                    baselines_array[array_time_loc:array_time_loc+num_baselines] = array(template_data['BASELINE'])
                    
                    ##Fill in the fractional julian date, after adding on the appropriate amount of
                    ##time - /(24*60*60) because julian number is a fraction of a whole day
                    adjust_float_jd_array = float_jd_array + (tstep / (24.0*60.0*60.0))
                    date_array[array_time_loc:array_time_loc+num_baselines] = adjust_float_jd_array
            
    ##If we want other sky model behaviours, i.e. curvature to the spectrum,
    ##must generate a sky model for every fine channel. Two methods: RTS style
    ##extrapolation between points, or use a fit of a 2nd order polynomial
    ##TODO method two doesn't exist currently
    else:
        ###For frequency channel in band
        for chan in good_chans:
            ##Take the band base_freq and add on fine channel freq
            freq = base_freq + (chan*ch_width)
            prefix_name = "%s_%.3f" %(outname,freq/1e+6)
            
            ##Create ini file to run oskar
            sky_osm_name = "%s_%.3f.osm" %(outname,freq/1e+6)
           
            num_time_steps = len(tsteps)
            obs_time_length = num_time_steps * time_int

            make_ini(prefix_name=prefix_name,ra=initial_ra_point,dec=MWA_LAT,freq=freq,start_time=oskar_date,
                    sky_osm_name=sky_osm_name,healpix=healpix,num_channels=1,
                    template_ini=template_ini,ch_width=ch_width,time_int=time_int,
                    telescope_dir=telescope_dir,num_time_steps=num_time_steps,obs_time_length=obs_time_length)
            
            ##Run the simulation
            cmd = "oskar_sim_interferometer --quiet %s.ini" %prefix_name
            run_command(cmd)
            
            ##Read in the data directly from the binary file
            ##Only one freq channel so number of vis is number of baselines
            #uu,vv,ww,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims = read_oskar_binary(filename="%s.vis" %prefix_name,num_vis=num_baselines,num_baselines=num_baselines)

            ##Read in the data directly from the binary file
            ##This file contains a single freq channel with all time steps
            #num_vis = num_baselines * num_oskar_channels * num_time_steps

            uu,vv,ww,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims = read_oskar_binary(filename="%s.vis" %prefix_name,num_time_steps=num_time_steps,num_channels=1,num_baselines=num_baselines)

            ##Clean up the oskar outputs
            cmd = "rm %s.ini %s.vis" %(prefix_name,prefix_name)
            run_command(cmd)

            ###For each time step
            for time_ind,tstep in enumerate(tsteps):
                ##Precess ra by time since the beginning of the observation 
                ##(convert time to angle, change from seconds to degrees)
                ##Include half of the time step
                time = add_time(oskar_date,tstep)
                
                ##Need to know where in the uvfits file structure to put data
                array_time_loc = num_baselines*time_ind
                time_ind_lower = time_ind*num_baselines
                
                ##Stored in metres in binary, convert to wavelengths
                chan_uu = uu[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
                chan_vv = vv[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
                chan_ww = ww[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
                
                ##Select the correct visibilities from the binary data
                chan_xx_res = xx_res[time_ind_lower:time_ind_lower+num_baselines]
                chan_xx_ims = xx_ims[time_ind_lower:time_ind_lower+num_baselines]
                chan_xy_res = xy_res[time_ind_lower:time_ind_lower+num_baselines]
                chan_xy_ims = xy_ims[time_ind_lower:time_ind_lower+num_baselines]
                chan_yx_res = yx_res[time_ind_lower:time_ind_lower+num_baselines]
                chan_yx_ims = yx_ims[time_ind_lower:time_ind_lower+num_baselines]
                chan_yy_res = yy_res[time_ind_lower:time_ind_lower+num_baselines]
                chan_yy_ims = yy_ims[time_ind_lower:time_ind_lower+num_baselines]
                
                ##Make complex numpy arrays
                comp_xx = make_complex(chan_xx_res,chan_xx_ims)
                comp_xy = make_complex(chan_xy_res,chan_xy_ims)
                comp_yx = make_complex(chan_yx_res,chan_yx_ims)
                comp_yy = make_complex(chan_yy_res,chan_yy_ims)
                
                ##Remove the phase tracking added in by OSKAR
                rotated_xx = rotate_phase(wws=chan_ww,visibilities=comp_xx)
                rotated_xy = rotate_phase(wws=chan_ww,visibilities=comp_xy)
                rotated_yx = rotate_phase(wws=chan_ww,visibilities=comp_yx)
                rotated_yy = rotate_phase(wws=chan_ww,visibilities=comp_yy)
                
                ##Use the centre of the fine channel
                freq_cent = freq + (ch_width / 2.0)
                
                ##Populate the v_container at the correct spots
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,0,0] = real(rotated_xx)
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,0,1] = imag(rotated_xx)
                
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,1,0] = real(rotated_yy)
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,1,1] = imag(rotated_yy)
                
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,2,0] = real(rotated_xy)
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,2,1] = imag(rotated_xy)
                
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,3,0] = real(rotated_yx)
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,3,1] = imag(rotated_yx)
                
                ##Set the weights of everything to ones
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,0,2] = ones(len(chan_uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,1,2] = ones(len(chan_uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,2,2] = ones(len(chan_uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,3,2] = ones(len(chan_uu))
                
                ##Only set the u,v,w to the central frequency channel
                ##Header of uvfits has to match central_freq_chan + 1 (uvfits 1 ordered, python zero ordered)
                if chan == central_freq_chan:
                    central_freq_chan_value = freq_cent
                    ##u,v,w stored in seconds by uvfits files
                    chan_uu /= freq_cent
                    chan_vv /= freq_cent
                    chan_ww /= freq_cent
                    
                    uus[array_time_loc:array_time_loc+num_baselines] = chan_uu
                    vvs[array_time_loc:array_time_loc+num_baselines] = chan_vv
                    wws[array_time_loc:array_time_loc+num_baselines] = chan_ww
                    
                    ##Fill in the baselines using the first time and freq uvfits
                    baselines_array[array_time_loc:array_time_loc+num_baselines] = array(template_data['BASELINE'])
                    
                    ##Fill in the fractional julian date, after adding on the appropriate amount of
                    ##time - /(24*60*60) because julian number is a fraction of a whole day
                    adjust_float_jd_array = float_jd_array + (tstep / (24.0*60.0*60.0))
                    date_array[array_time_loc:array_time_loc+num_baselines] = adjust_float_jd_array
                
    output_uvfits_name = "%s/%s_t%02d_f%.3f_band%02d.uvfits" %(data_dir,outname,time_int,ch_width/1e+6,band_num)
    print "Now creating uvfits file %s" %output_uvfits_name
    create_uvfits(v_container=v_container,freq_cent=central_freq_chan_value,ra_point=initial_ra_point,
                output_uvfits_name=output_uvfits_name,uu=uus,vv=vvs,ww=wws,baselines_array=baselines_array,
                date_array=date_array,date=oskar_date,central_freq_chan=central_freq_chan,ch_width=ch_width,
                template_uvfits=template_uvfits,int_jd=int_jd)
    print "Written file %s" %output_uvfits_name

###Whack into a function so we can line_profile it
the_main_loop(tsteps=tsteps)

if options.osm:
    pass
else:
    for chan in good_chans:
        ##Take the band base_freq and add on fine channel freq
        freq = base_freq + (chan*ch_width)
        ##Create ini file to run oskar
        cmd = "rm %s_%.3f.osm" %(outname,freq/1e+6)
        run_command(cmd)
        
cmd = 'rm -r %s/telescope_%s_band%02d' %(tmp_dir,outname,band_num)
run_command(cmd)
        
chdir(cwd)
template_uvfits.close()
