#!/usr/bin/env python
from subprocess import call
from sys import exit
from optparse import OptionParser
from numpy import zeros, pi, sin, cos, real, imag, loadtxt, array, floor, arange, ones, where, savetxt, mod, load
from numpy import exp as n_exp
from ephem import Observer
from cmath import exp
from sys import path as sys_path
from scipy.interpolate import interp1d
try:
    from jdcal import gcal2jd
except ImportError:
    ##horrible gstar specific fix
    sys_path.append('/lustre/projects/p048_astro/jline/software/jdcal-1.3')
    from jdcal import gcal2jd
from os import environ,getcwd,chdir,makedirs,path
from struct import unpack
try:
    OSKAR_dir = environ['OSKAR_TOOLS']
except:
    ##horrible gstar specific fix
    OSKAR_dir = '/lustre/projects/p048_astro/jline/software/OSKAR_tools'
sys_path.append(OSKAR_dir)
from MWAobs_oskar_lib import *
from astropy.io import fits

R2D = 180.0 / pi
D2R = pi / 180.0
MWA_LAT = -26.7033194444
VELC = 299792458.0

parser = OptionParser()


parser.add_option('--telescope', default='%s/telescopes/MWA_phase1' %OSKAR_dir, help='Enter telescope used for simulation. Default = $OSKAR_TOOLS/telescopes/MWA_phase1')
parser.add_option('--band_num', help='Enter band number to simulate')
parser.add_option('--osm', default=False, help='Location of OSKAR osm sky model to use')
parser.add_option('--debug',default=False,action='store_true', help='Enable to debug with print statements')
parser.add_option('--healpix', default=False, help='Enter healpix tag to use base images NOT CURRENTLY USED')
parser.add_option('--fit_osm', default=False, help='Location of sky parameters to create osm from NOT CURRENTLY USED')
parser.add_option( '--ini_file', default=False, help='Enter template oskar .ini - defaults to the template .ini located in $OSKAR_TOOLS/telescopes/--telescope')
parser.add_option('--flag_dipoles',default=False,action='store_true', help='Add to switch on dipole flagging via the metafits file. NOTE needs a metafits that has the correct data within')
parser.add_option('--do_phase_track',default=False,action='store_true',
    help='Add to leave the default phase tracking done by OSKAR (phase track the RA,DEC in the metafits). OPTIONAL - can also set the phase centre explicitly using --phase_centre')
parser.add_option('--phase_centre',default=False,
    help='Set phase centre and leave in phase tracking in final uvfits. Usage: --phase_centre=ra,dec with ra,dec in deg')

parser.add_option('--data_dir', help='Where to output the finished uvfits - default is ./data',default=False)
parser.add_option('--metafits', help='Enter name of metafits file to base obs on')
parser.add_option('--output_name', help='Enter prefix name for outputs')
parser.add_option('--srclist', default=False, help='Enter location and name of the RTS srclist to use as a sky model')
parser.add_option('--time', help='Enter start,end of sim in seconds from the beginning of the observation (as set by metafits)')
parser.add_option('--time_int', default=False, help='Enable to force a different time cadence from that in the metafits - enter the time in seconds')
parser.add_option('--freq_int', default=False, help='Enable to force a different fine channel width from that in the metafits - enter the frequency in Hz')
parser.add_option('--chips_settings', default=False, action='store_true',
    help='Swtiches on a default CHIPS resolution and uvfits weightings - 8s, 80kHz integration with the normal 5 40kHz channels missing. OVERRIDES other time/freq int settings')
parser.add_option('--full_chips', default=False, action='store_true',
    help='Instead of missing freq channels, do a complete simulation when doing a CHIPS simulation')
parser.add_option('--all_chans', default=False, action='store_true',
    help='Instead of missing fine freq channels, do a complete 32 chan simulation')

parser.add_option('--oskar_gsm', default=False, action='store_true',help='Add to include the gsm as calculated by OSKAR')
parser.add_option('--oskar_gsm_SI', default=-2.5, help='The spectral index to give the gsm made by OSKAR')
parser.add_option('--oskar_gsm_file', default='%s/gsm/component_maps_408locked.dat' %OSKAR_dir, help='The de Oliveira-Costa model file to use - defaults to $OSKAR_TOOLS/gsm/component_maps_408locked.dat')
parser.add_option('--retain_vis_file',default=False,action='store_true', help='Add to not delete the oskar binary .vis files')
parser.add_option('--retain_ini_file',default=False,action='store_true', help='Add to not delete the oskar binary .ini files')

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

for key in ['DATE-OBS','FREQCENT','FINECHAN','INTTIME','BANDWDTH','AZIMUTH','ALTITUDE','RA','DEC']:
    test_avail(key)

##Get the east, north, height antenna positions from the metafits
##Tile positions are stored for both XX/YY pols so only
##want to select half of them
east = f[1].data['East']
north = f[1].data['North']
height = f[1].data['Height']

##Create and fill a layout array
layout_array = zeros((len(east)/2,3))
selection = arange(0,len(east),2)

##There is some kind of co-ordinate difference between the RTS and OSKAR, so
##need to make all coords negative here
## (trial and error found this gives the expected u,v,w coords)
layout_array[:,0] = -east[selection]
layout_array[:,1] = -north[selection]
layout_array[:,2] = -height[selection]

initial_date = f[0].header['DATE-OBS']
##Change in to oskar date format
date,time = initial_date.split('T')
year,month,day = date.split('-')
oskar_date = "%s-%s-%s %s" %(day,month,year,time)

time_int = float(f[0].header['INTTIME'])
if options.time_int: time_int = float(options.time_int)

ch_width = float(f[0].header['FINECHAN'])*1e+3
freqcent = float(f[0].header['FREQCENT'])*1e+6
b_width = float(f[0].header['BANDWDTH'])*1e+6
base_low_freq = freqcent - (b_width/2) - (ch_width/2)
low_freq = base_low_freq

if options.freq_int:
    ch_width = float(options.freq_int)
    low_freq = base_low_freq + (ch_width / 2.0)

if options.chips_settings:
    # print("HERE MAAAAAN WOT")
    ch_width = 80e+3
    time_int = 8.0
    low_freq = base_low_freq - (ch_width / 2.0)

##ephem Observer class, use this to compute LST from the date of the obs
MRO = Observer()
##Set the observer at Boolardy
MRO.lat, MRO.long, MRO.elevation = '-26:42:11.95', '116:40:14.93', 0
date,time = initial_date.split('T')
MRO.date = '/'.join(date.split('-'))+' '+time
initial_lst = float(MRO.sidereal_time())*R2D
#initial_ra_point = float(MRO.sidereal_time())*R2D
#dec_point = MWA_LAT
initial_ra_point = float(f[0].header['RA'])
dec_point = float(f[0].header['DEC'])

healpix = options.healpix
telescope_dir = options.telescope
telescope_name = options.telescope.split('/')[-1]
template_uvfits = fits.open("%s/template_%s.uvfits" %(telescope_dir,telescope_name))
template_data = template_uvfits[0].data
template_baselines = template_uvfits[0].data['BASELINE'].copy()
num_baselines = len(template_data)
num_freq_channels = int(1.28e+6 / ch_width)

X,Y,Z = enh2xyz(east[selection], north[selection],height[selection])

template_uvfits[1].data['STABXYZ'][:,0] = X
template_uvfits[1].data['STABXYZ'][:,1] = Y
template_uvfits[1].data['STABXYZ'][:,2] = Z

x_lengths = []
y_lengths = []
z_lengths = []

for ant1 in arange(len(X) - 1):
    for ant2 in arange(ant1+1,len(X)):
        x_length = X[ant1] - X[ant2]
        y_length = Y[ant1] - Y[ant2]
        z_length = Z[ant1] - Z[ant2]
        x_lengths.append(x_length)
        y_lengths.append(y_length)
        z_lengths.append(z_length)


x_lengths = array(x_lengths)
y_lengths = array(y_lengths)
z_lengths = array(z_lengths)

oskar_gsm = options.oskar_gsm
oskar_gsm_file = options.oskar_gsm_file
oskar_gsm_SI = float(options.oskar_gsm_SI)

if options.ini_file:
    template_ini = options.ini_file
else:
    template_ini = "%s/template_%s.ini" %(telescope_dir,telescope_name)
template_ini = open(template_ini).read().split('\n')

##Unflagged channel numbers for 40kHz observation
if num_freq_channels == 32:
    if options.all_chans:
        good_chans = arange(32)
    else:
        good_chans = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29]
    central_freq_chan = 15
##Ignores first and last channels for CHIPS settings
elif options.chips_settings:
    if options.full_chips:
        good_chans = range(0,16)
    else:
        good_chans = range(1,15)
    central_freq_chan = 8
##Anything else just simulate them all
else:
    good_chans = range(0,num_freq_channels)
    central_freq_chan = num_freq_channels / 2

#good_chans = xrange(32)
#good_chans = [2,3]
#central_freq_chan = 2

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

if options.do_phase_track:
    undo_phase_track = False
else:
    undo_phase_track = True

if options.phase_centre:
    initial_ra_point,dec_point = map(float,options.phase_centre.split(','))
    undo_phase_track = False


ha_point = initial_lst - initial_ra_point

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
    create_flagged_telescope(metafits=options.metafits,layout_array=layout_array,azimuth=azimuth,altitude=altitude,
                                telescope_dir=telescope_dir,new_telescope_dir='%s/telescope_%s_band%02d' %(tmp_dir,outname,band_num))
    ##Reset the telescope_dir to copied location
    telescope_dir = '%s/telescope_%s_band%02d' %(tmp_dir,outname,band_num)
else:
    ##Just copy the template telescope
    cmd = 'cp -r %s %s/telescope_%s_band%02d' %(telescope_dir,tmp_dir,outname,band_num)
    run_command(cmd)

    savetxt('%s/telescope_%s_band%02d/layout.txt' %(tmp_dir,outname,band_num),layout_array)

    ##Reset the telescope_dir to copied location
    telescope_dir = '%s/telescope_%s_band%02d' %(tmp_dir,outname,band_num)

    #OSKAR naturally beam forms towards phase centre
    #Read the pointing from the metafits and force observation towards single pointing
    permitted = open('%s/station/permitted_beams.txt' %(telescope_dir),'w+')
    permitted.write('%.5f %.5f' %(azimuth,altitude))
    permitted.close

    permitted = open('%s/permitted_beams.txt' %(telescope_dir),'w+')
    permitted.write('%.5f %.5f' %(azimuth,altitude))
    permitted.close


##The OSKAR analytic beam can only be normalised to
##the phase centre - this is bad for no phase tracking,
##as the beam normalisation changes with time for a stationary
##beam. To solve for this, I've made images of the beam created
##from 50 - 300MHz at 0.5MHz intervals, and measured the max gain
##From this we can fit a spline function, and empirically get a gain
##correction, so all beams are normalised to 1 at zenith

beam_info = load('%s/OSKAR_beam_gains.npz' %(telescope_dir))['beam_info']

freqs = beam_info[:,0]
maxs = beam_info[:,1]

##This is now a function of freq - pass it to add_data_to_uvfits to
##do the correction when adding data to uvfits
beam_gain = interp1d(freqs,array(maxs),kind='cubic')


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

        oskar_inds = [oskar_channels.index(good_chan) for good_chan in good_chans]
        num_oskar_channels = len(oskar_channels)
        #print "NUM CHANNELS",num_channels
        osk_low_freq = base_freq + (good_chans[0]*ch_width)
        #osk_low_freq = base_freq

        num_time_steps = len(tsteps)
        obs_time_length = num_time_steps * time_int

        #oskar_date

        make_ini(prefix_name=prefix_name,ra=initial_ra_point,dec=dec_point,freq=osk_low_freq,start_time=oskar_date,
                    sky_osm_name=options.osm,healpix=healpix,num_channels=num_oskar_channels,
                    template_ini=template_ini,ch_width=ch_width,time_int=time_int,
                    telescope_dir=telescope_dir,num_time_steps=num_time_steps,obs_time_length=obs_time_length,
                    oskar_gsm=oskar_gsm,oskar_gsm_file=oskar_gsm_file,oskar_gsm_SI=oskar_gsm_SI)

        # ##Run the simulation
        # if undo_phase_track:
        #     ##Run the simulation
        #     environ["LD_LIBRARY_PATH"] = environ["OSK_NOTRACK_LD"]
        #     cmd = "%s/oskar_sim_interferometer --quiet %s.ini" %(environ["OSK_NOTRACK_BIN"],prefix_name)
        #     run_command(cmd)
        #
        # else:
        ##Run the simulation
        environ["LD_LIBRARY_PATH"] = environ["OSK_TRACK_LD"]
        cmd = "%s/oskar_sim_interferometer --quiet %s.ini" %(environ["OSK_TRACK_BIN"],prefix_name)
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
        if options.retain_vis_file and options.retain_ini_file:
            pass
        else:
            if options.retain_vis_file:
                cmd = "rm %s.ini" %(prefix_name)
            elif options.retain_ini_file:
                cmd = "rm %s.vis" %(prefix_name)
            else:
                cmd = "rm %s.ini %s.vis" %(prefix_name,prefix_name)
            run_command(cmd)

        ###For each time step
        for time_ind,tstep in enumerate(tsteps):
            ##Precess ra by time since the beginning of the observation
            ##(convert time to angle, change from seconds to degrees)
            ##Include half of the time step
            time = add_time(oskar_date,tstep)

            ##If we are using a sky model with fixed spectral indexes, can run all
            ##frequencies using the same osm model, and so only need to run OSKAR once
            ##per time step
            for oskar_ind,chan in zip(oskar_inds,good_chans):
                freq = base_freq + (chan*ch_width)
                freq_cent = freq + (ch_width / 2.0)

                if chan == central_freq_chan:

                    central_freq_chan_value = freq_cent

                add_data_to_uvfits(v_container=v_container,time_ind=time_ind,num_baselines=num_baselines,template_baselines=template_baselines,
                    chan=chan,num_oskar_channels=num_oskar_channels,oskar_ind=oskar_ind,uu=uu,vv=vv,ww=ww,float_jd_array=float_jd_array,
                    baselines_array=baselines_array,date_array=date_array,undo_phase_track=undo_phase_track,xx_res=xx_res,xx_ims=xx_ims,
                    xy_res=xy_res,xy_ims=xy_ims,yx_res=yx_res,yx_ims=yx_ims,yy_res=yy_res,yy_ims=yy_ims,
                    x_lengths=x_lengths,y_lengths=y_lengths,z_lengths=z_lengths,uus=uus,vvs=vvs,wws=wws,
                    tstep=tstep,freq=freq,ch_width=ch_width,central_freq_chan=central_freq_chan,chips_settings=options.chips_settings,
                    full_chips=options.full_chips,gain_correction=beam_gain,
                    ha_point=ha_point*D2R,dec_point=dec_point*D2R)


    ##If we want other sky model behaviours, i.e. curvature to the spectrum,
    ##must generate a sky model for every fine channel. Two methods: RTS style
    ##extrapolation between points, or use a fit of a 2nd order polynomial
    ##TODO method two doesn't exist currently
    else:
        ###For frequency channel in band
        for chan in good_chans:

            if chan == central_freq_chan:
                freq_cent = freq + (ch_width / 2.0)
                central_freq_chan_value = freq_cent

            ##Take the band base_freq and add on fine channel freq
            freq = base_freq + (chan*ch_width)
            prefix_name = "%s_%.3f" %(outname,freq/1e+6)

            ##Create ini file to run oskar
            sky_osm_name = "%s_%.3f.osm" %(outname,freq/1e+6)

            num_time_steps = len(tsteps)
            obs_time_length = num_time_steps * time_int

            make_ini(prefix_name=prefix_name,ra=initial_ra_point,dec=dec_point,freq=freq,start_time=oskar_date,
                    sky_osm_name=sky_osm_name,healpix=healpix,num_channels=1,
                    template_ini=template_ini,ch_width=ch_width,time_int=time_int,
                    telescope_dir=telescope_dir,num_time_steps=num_time_steps,obs_time_length=obs_time_length,
                    oskar_gsm=oskar_gsm,oskar_gsm_file=oskar_gsm_file,oskar_gsm_SI=oskar_gsm_SI)

            # ##Run the simulation
            # if undo_phase_track:
            #     ##Run the simulation
            #     environ["LD_LIBRARY_PATH"] = environ["OSK_NOTRACK_LD"]
            #     cmd = "%s/oskar_sim_interferometer --quiet %s.ini" %(environ["OSK_NOTRACK_BIN"],prefix_name)
            #     run_command(cmd)
            #
            # else:
            #     ##Run the simulation
            environ["LD_LIBRARY_PATH"] = environ["OSK_TRACK_LD"]
            cmd = "%s/oskar_sim_interferometer --quiet %s.ini" %(environ["OSK_TRACK_BIN"],prefix_name)
            run_command(cmd)

            ##Read in the data directly from the binary file
            ##Only one freq channel so number of vis is number of baselines
            #uu,vv,ww,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims = read_oskar_binary(filename="%s.vis" %prefix_name,num_vis=num_baselines,num_baselines=num_baselines)

            ##Read in the data directly from the binary file
            ##This file contains a single freq channel with all time steps
            #num_vis = num_baselines * num_oskar_channels * num_time_steps

            uu,vv,ww,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims = read_oskar_binary(filename="%s.vis" %prefix_name,num_time_steps=num_time_steps,num_channels=1,num_baselines=num_baselines)

            ##Clean up the oskar outputs
            if options.retain_vis_file and options.retain_ini_file:
                pass
            else:
                if options.retain_vis_file:
                    cmd = "rm %s.ini" %(prefix_name)
                elif options.retain_ini_file:
                    cmd = "rm %s.vis" %(prefix_name)
                else:
                    cmd = "rm %s.ini %s.vis" %(prefix_name,prefix_name)
                run_command(cmd)

            ###For each time step
            for time_ind,tstep in enumerate(tsteps):
                ##Precess ra by time since the beginning of the observation
                ##(convert time to angle, change from seconds to degrees)
                ##Include half of the time step
                add_data_to_uvfits(v_container=v_container,time_ind=time_ind,num_baselines=num_baselines,template_baselines=template_baselines,
                    chan=chan,uu=uu,vv=vv,ww=ww,float_jd_array=float_jd_array,baselines_array=baselines_array,date_array=date_array,
                    undo_phase_track=undo_phase_track,xx_res=xx_res,xx_ims=xx_ims,xy_res=xy_res,xy_ims=xy_ims,yx_res=yx_res,
                    yx_ims=yx_ims,yy_res=yy_res,yy_ims=yy_ims,x_lengths=x_lengths,y_lengths=y_lengths,z_lengths=z_lengths,
                    uus=uus,vvs=vvs,wws=wws,tstep=tstep,freq=freq,ch_width=ch_width,central_freq_chan=central_freq_chan,
                    chips_settings=options.chips_settings,full_chips=options.full_chips,gain_correction=beam_gain,
                    ha_point=ha_point*D2R,dec_point=dec_point*D2R)

    if options.chips_settings:
        output_uvfits_name = "%s/%s_chips-t%02d_f%.3f_band%02d.uvfits" %(data_dir,outname,time_int,ch_width/1e+6,band_num)
    else:
        output_uvfits_name = "%s/%s_t%02d_f%.3f_band%02d.uvfits" %(data_dir,outname,time_int,ch_width/1e+6,band_num)
    print "Now creating uvfits file %s" %output_uvfits_name
    create_uvfits(v_container=v_container,freq_cent=central_freq_chan_value,ra_point=initial_ra_point,dec_point=dec_point,
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
