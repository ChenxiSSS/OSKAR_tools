#!/usr/bin/env python
from subprocess import call
from sys import exit
from optparse import OptionParser
from os import environ,getcwd,chdir,makedirs,path
from numpy import zeros, pi, sin, cos, real, imag, loadtxt, array, floor, arange, ones, where, mod, ndarray
from numpy import exp as n_exp
from ephem import Observer
from cmath import exp
from jdcal import gcal2jd
from struct import unpack
try:
    import pyfits as fits
except ImportError:
    from astropy.io import fits
import matplotlib.pyplot as plt

OSKAR_dir = environ['OSKAR_TOOLS']
R2D = 180.0 / pi
D2R = pi / 180.0
MWA_LAT = -26.7033194444
VELC = 299792458.0

parser = OptionParser()

parser.add_option('-a','--telescope', default='%s/telescopes/MWA_phase1' %OSKAR_dir, help='Enter telescope used for simulation. Default = $OSKAR_TOOLS/telescopes/MWA_phase1')
parser.add_option('-b','--band_num', help='Enter band number to simulate')
parser.add_option('-c','--osm', default=False, help='Location of OSKAR osm sky model to use')
parser.add_option('-d','--debug',default=False,action='store_true', help='Enable to debug with print statements')

parser.add_option('-f','--healpix', default=False, help='Enter healpix tag to use base images')
parser.add_option('-g','--fit_osm', default=False, help='Location of sky parameters to create osm from')
parser.add_option('-i', '--ini_file', default=False, help='Enter template oskar .ini - defaults to the template .ini located in $OSKAR_TOOLS/telescopes/--telescope')
parser.add_option('-m','--metafits', help='Enter name of metafits file to base obs on')
parser.add_option('-n','--output_name', help='Enter prefix name for outputs')
parser.add_option('-o','--data_dir', help='Where to output the finished uvfits - default is ./data',default=False)
parser.add_option('-p','--phase_centre', help='Set the phase centre of the output visibilities. Write as RA,DEC in deg',default=False)
parser.add_option('-s','--srclist', default=False, help='Enter location and name of the RTS srclist to use as a sky model')
parser.add_option('-t','--time', help='Enter start,end of sim in seconds from the beginning of the observation (as set by metafits)')
parser.add_option('-x','--time_int', default=False, help='Enable to force a different time cadence - enter the time in seconds')

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

for key in ['DATE-OBS','FREQCENT','FINECHAN','INTTIME','BANDWDTH','DEC','RA']:
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
#intial_ra_point = float(MRO.sidereal_time())*R2D
intial_ra_point = float(f[0].header['RA'])
dec_point = float(f[0].header['DEC'])

healpix = options.healpix
telescope_dir = options.telescope
telescope_name = options.telescope.split('/')[-1]
template_uvfits = fits.open("%s/template_%s.uvfits" %(telescope_dir,telescope_name))
template_data = template_uvfits[0].data
num_baselines = len(template_data)
num_freq_channels = 32

##Need to compute new u,v,w coords if rephasing - use X,Y,Z coords to do so
if options.phase_centre:
    antennas = {}
    xyz = template_uvfits[1].data['STABXYZ'].copy()
    for i in xrange(len(xyz)):
        antennas['ANT%03d' %(i+1)] = xyz[i]
    baselines = template_uvfits[0].data['BASELINE'].copy()

    x_lengths = []
    y_lengths = []
    z_lengths = []

    for baseline in baselines:
        ant2 = mod(baseline, 256)
        ant1 = (baseline - ant2)/256
        x_length,y_length,z_length = antennas['ANT%03d' %ant1] - antennas['ANT%03d' %ant2]
        x_lengths.append(x_length)
        y_lengths.append(y_length)
        z_lengths.append(z_length)
        
    final_phase_ra,final_phase_dec = map(float,options.phase_centre.split(','))


#antenna_table = template_uvfits[1].data

if options.ini_file:
    template_ini = options.ini_file
else:
    template_ini = "%s/template_%s.ini" %(telescope_dir,telescope_name)
template_ini = open(template_ini).read().split('\n')
    
##Unflagged channel numbers
good_chans = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29]
#good_chans = xrange(32)
#good_chans = [2]

##Flagged channel numbers
#bad_chans = [0,1,16,30,31]

central_freq_chan = 15

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

def get_uvw_freq_local(x_length=None,y_length=None,z_length=None,dec=None,ha=None,wavelength=None):
    '''Takes the baseline length in meters and uses the frequency to calc
    u,v,w in units of wavelength'''
    
    X = x_length / wavelength
    Y = y_length / wavelength
    Z = z_length / wavelength
    
    u = sin(ha)*X + cos(ha)*Y
    v = -sin(dec)*cos(ha)*X + sin(dec)*sin(ha)*Y + cos(dec)*Z
    w = cos(dec)*cos(ha)*X - cos(dec)*sin(ha)*Y + sin(dec)*Z
    
    return u,v,w
    
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
    return '%02d-%02d-%d %02d:%02d:%05.2f' %(day,month,year,int(hours),int(mins),secs)

def calc_jdcal(date):
    dmy, hms = date.split()
    
    day,month,year = map(int,dmy.split('-'))
    hour,mins,secs = map(float,hms.split(':'))

    ##For some reason jdcal gives you the date in two pieces
    ##Gives you the time up until midnight of the day
    jd1,jd2 = gcal2jd(year,month,day)
    jd3 = (hour + (mins / 60.0) + (secs / 3600.0)) / 24.0

    jd = jd1 + jd2 + jd3
    
    jd_day = jd1 + floor(jd2)
    jd_fraction = (jd2 - floor(jd2)) + jd3
    
    ##The header of the uvdata file takes the integer, and
    ##then the fraction goes into the data array for PTYPE5
    return jd_day, jd_fraction

##OSKAR phase tracks, but the MWA correlator doesn't, so we have to undo
##any phase tracking done
def rotate_phase(wws=None,visibilities=None,new_wws=None):
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
    if type(new_wws) != ndarray:
        phase_rotate = n_exp( PhaseConst * wws)
        rotated_visis = visibilities * phase_rotate
    ##however, if we want to phase to a new phase centre, must apply w terms in
    ##that direction
    else:
        phase_rotate = n_exp(PhaseConst * wws)
        rotated_visis = visibilities * phase_rotate
        
        phase_rotate = n_exp(PhaseConst * -new_wws)
        rotated_visis *= phase_rotate
        #print phase_rotate
    
    
    #rotated_visis = visibilities
    
    return rotated_visis

def get_uvw_data(d_single=None,d_double=None,num_baselines=None,block=None):
    '''Takes the format of the OSKAR binary data, and pulls out the u, v, or w data
    Returns the data in metres'''
    uvw_data = zeros(num_baselines)
    
    if d_single:
        ##single precision, number is 4 bytes long
        uvw_inds = arange(0,4*num_baselines,4)
    elif d_double:
        ##double precision, number is 8 bytes long
        uvw_inds = arange(0,8*num_baselines,8)
    else:
        print "J Line hasn't coded to deal with the data type in the OSKAR binary. Have coded in single and double precision (float and double) but not int or char. Exiting"
        exit(0)
    
    for ind in xrange(num_baselines):
        if d_single:
            uvw_data[ind] = unpack('<f',block[uvw_inds[ind]:uvw_inds[ind]+4])[0]
        elif d_double:
            uvw_data[ind] = unpack('<d',block[uvw_inds[ind]:uvw_inds[ind]+8])[0]
    
    return uvw_data

def read_oskar_binary(filename=None,num_vis=None,num_baselines=None):
    '''Open the native OSKAR binary format and extracts the visibilities.
    Returns the u,v,w coords and complex XX,XY,YX,YY arrays
    Will ONLY work with the default OSKAR settings found in the template .ini
    files. TODOs are needs to find the data element type (e.g. float or double),
    the dimensions of the data etc. Currently hard-coded to defaults'''
    
    ##The binary format is detailed here: http://oskar.oerc.ox.ac.uk/sites/default/files/2.6/doc/OSKAR-2.6-Binary-File-Format.pdf
    
    ##This opens the file up into a massive string I think,
    ##with each character representing a byte of data
    binary_file = open(filename, mode='rb').read()

    ##OSKAR splits the header as first 64 bytes
    ##The rest in the main body
    #header = binary_file[:64]
    #binary_body = binary_file[64:]

    ###Should tell us OSKAR binary format
    #binary_format = header[9]
    #print binary_format

    ##Gotta find the start of each chunk using the block size from
    ##the tag of each chunk - found from last 8 bytes of the tag
    
    ##First tag happens at byte 65
    tag_indexes = [64]
    tag_index = 64
    block_size = 0
    next_block_size = 0
    
    ##Loops over each chunk (tag + block), reading how long
    ##each block is using the last 4 bytes of the tag
    ##Essentially each tag spilts the data up into different data types
    ##e.g. metadata, observation settings, actually visibilities etc
    ##Each tag is 20 bytes in size
    while tag_index+next_block_size+20 < len(binary_file):
        this_tag = binary_file[tag_index:tag_index+20]
        block_size = unpack('<Q',this_tag[12:])[0]
        tag_index += (block_size + 20)
        tag_indexes.append(tag_index)
        
        next_tag = binary_file[tag_index:tag_index+20]
        next_block_size = unpack('<Q',next_tag[12:])[0]

    ##emtpy arrays a needed values
    #wavelength = VELC / freq
    uu = zeros(num_vis)
    vv = zeros(num_vis)
    ww = zeros(num_vis)
    
    ##Each tag_index marks where one chunk starts - one chunk = one tag + one block
    for tag_index in tag_indexes:
        ##OSKAR tag is 20 bytes long, including 'TGB'. Ignore first 3 characters
        tag =  binary_file[tag_index+3:tag_index+20]
        
        ##Use unpack to turn binary into numbers
        ##The first arg identifies the type of number
        ## '<Q' = little-endian 8-byte integer
        ## '=b' = native signed char
        ## '<I' = little-endian 4-byte integer

        ##Get data info from the tag
        size_of_one_element_of_payload_data = unpack('=b',tag[0])[0]
        chunk_flag = unpack('=b',tag[1])[0]
        
        ##Each bit of the chunk_flag (1 byte) encodes info.
        ##Convert to a string, and map to ints. I think the order
        ##of the string is opposite of order written in OSKAR docs
        ##I think.
        extended_tag,crc_code,payload_endian,meh1,meh2,meh3,meh4,meh5 = map(int,format(chunk_flag,'08b'))
        
        ##Don't actually know what the extended tag does - hopefully never shows up!
        if extended_tag:
            group_name_size_in_bytes = unpack('=b',tag[3])[0]
            tag_name_size_in_bytes = unpack('=b',tag[4])[0]
        else:
            group_ID = unpack('=b',tag[3])[0]
            tag_ID = unpack('=b',tag[4])[0]

        data_type_code = unpack('=b',tag[2])[0]
        meh1,d_matrix,d_complex,meh2,d_double,d_single,d_int,d_char = map(int,format(data_type_code,'08b'))
        ##The bits above specify the data type, e.g. single or double precision
        ##Changes the number of bytes a number is stored by and is read in
        ##by unpack
        
        user_specified_index = unpack('<I',tag[5:9])[0]
        block_size = unpack('<Q',tag[9:])[0]
        ##CRC code is some integer that helps you check for corrupted data
        ##No idea how to read it, but want to cut it off the end if it's there
        if crc_code:
            block =  binary_file[tag_index+20:tag_index+20+block_size-4]
        else:
            block =  binary_file[tag_index+20:tag_index+20+block_size]
        
        ##In the OSKAR binary format, u,v,w stored under group_ID = 12, as tag_ID 4,5,6
        ##respectively. Stored in metres so divide by wavelength for natural units
        ##Each number is 4 bytes long
        
        if group_ID == 12 and tag_ID == 4:
            uu = get_uvw_data(d_single=d_single,d_double=d_double,num_baselines=num_baselines,block=block)
            
        elif group_ID == 12 and tag_ID == 5:
            vv = get_uvw_data(d_single=d_single,d_double=d_double,num_baselines=num_baselines,block=block)
            
        elif group_ID == 12 and tag_ID == 6:
            ww = get_uvw_data(d_single=d_single,d_double=d_double,num_baselines=num_baselines,block=block)
        
        ##In the OSKAR binary format, XX,YY,XY,YX are group_ID = 12, tag_ID = 3
        elif group_ID == 12 and tag_ID == 3:
            
            ##*8 because we 4 complex numbers so 4 real 4 imag
            all_nums = zeros(8*num_vis)
            
            if d_single:
                ##*4 because floats are 4 bytes long
                visi_inds = arange(0,8*4*num_vis,4)
            elif d_double:
                ##*4 because floats are 8 bytes long
                visi_inds = arange(0,8*8*num_vis,8)
            
            for ind in xrange(8*num_vis):
                if d_single:
                    all_nums[ind] = unpack('<f',block[visi_inds[ind]:visi_inds[ind]+4])[0]
                elif d_double:
                    all_nums[ind] = unpack('<f',block[visi_inds[ind]:visi_inds[ind]+8])[0]
            
            ##Split the data up by polarisation
            xx_res = all_nums[arange(0,8*num_vis,8)]
            xx_ims = all_nums[arange(1,8*num_vis,8)]
            xy_res = all_nums[arange(2,8*num_vis,8)]
            xy_ims = all_nums[arange(3,8*num_vis,8)]
            yx_res = all_nums[arange(4,8*num_vis,8)]
            yx_ims = all_nums[arange(5,8*num_vis,8)]
            yy_res = all_nums[arange(6,8*num_vis,8)]
            yy_ims = all_nums[arange(7,8*num_vis,8)]
            
        else:
            pass
            
    return uu,vv,ww,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims
        
def make_complex(re=None,im=None):
    '''Takes two arrays, and returns a complex array with re real values and im imarginary values'''
    comp = array(re,dtype=complex)
    comp += 1j * im
    
    return comp

def create_uvfits(freq_cent=None, ra_point=intial_ra_point, dec_point=dec_point, output_uvfits_name=None,uu=None,vv=None,ww=None,baselines_array=None,date_array=None,date=oskar_date):
    
    ##UU, VV, WW don't actually get read in by RTS - might be an issue with
    ##miriad/wsclean however, as it looks like oskar w = negative maps w
    uvparnames = ['UU','VV','WW','BASELINE','DATE']
    parvals = [uu,vv,ww,baselines_array,date_array]
        
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
    uvhdu.header['CRVAL4'] = freq_cent  ##Middle pixel value in Hz
    uvhdu.header['CRPIX4'] = int(central_freq_chan) + 1 ##Middle pixel number
    uvhdu.header['CDELT4'] = ch_width

    uvhdu.header['CTYPE5'] = template_uvfits[0].header['CTYPE5']
    uvhdu.header['CRVAL5'] = template_uvfits[0].header['CRVAL5']
    uvhdu.header['CRPIX5'] = template_uvfits[0].header['CRPIX5']
    uvhdu.header['CDELT5'] = template_uvfits[0].header['CDELT5']

    uvhdu.header['CTYPE6'] = template_uvfits[0].header['CTYPE6']
    uvhdu.header['CRVAL6'] = ra_point
    uvhdu.header['CRPIX6'] = template_uvfits[0].header['CRPIX6']
    uvhdu.header['CDELT6'] = template_uvfits[0].header['CDELT6']

    uvhdu.header['CTYPE7'] = template_uvfits[0].header['CTYPE7']
    uvhdu.header['CRVAL7'] = dec_point
    uvhdu.header['CRPIX7'] = template_uvfits[0].header['CRPIX7']
    uvhdu.header['CDELT7'] = template_uvfits[0].header['CDELT7']

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

    template_uvfits[1].header['FREQ'] = freq_cent
    
    ##MAJICK uses this date to set the LST
    dmy, hms = date.split()
    day,month,year = map(int,dmy.split('-'))
    hour,mins,secs = map(float,hms.split(':'))
    
    rdate = "%d-%02d-%02dT%02d:%02d:%.2f" %(year,month,day,hour,mins,secs)
    
    template_uvfits[1].header['RDATE'] = rdate

    ## Create hdulist and write out file
    hdulist = fits.HDUList(hdus=[uvhdu,template_uvfits[1]])
    hdulist.writeto(output_uvfits_name,clobber=True)
    hdulist.close()
    
def make_ini(prefix_name=None,ra=None,dec=None,freq=None,start_time=None,sky_osm_name=None,healpix=None,num_channels=None):
    out_file = open("%s.ini" %prefix_name,'w+')
    for line in template_ini:
        if "num_channels" in line:
            line = "num_channels=%d" %num_channels
        if "start_frequency_hz" in line:
            line = "start_frequency_hz=%.10f" %freq
        elif "frequency_inc_hz" in line:
            line = "frequency_inc_hz=%.10f" %ch_width
        elif "phase_centre_ra_deg" in line:
            line = "phase_centre_ra_deg=%.10f" %ra
            #line = "phase_centre_ra_deg=0.0"
        elif "phase_centre_dec_deg" in line:
            line = "phase_centre_dec_deg=%.10f" %dec
        elif "start_time_utc" in line:
            line = "start_time_utc=%s" %start_time
        elif line.split('=')[0]=="length":
            line = "length=%.10f" %time_int
        elif "oskar_vis_filename" in line:
            line = "oskar_vis_filename=%s.vis" %prefix_name
        elif "channel_bandwidth_hz" in line:
            line = "channel_bandwidth_hz=%.10f" %ch_width
        elif "time_average_sec" in line:
           line = "time_average_sec=%.10f" %time_int
        #elif "ms_filename" in line:
            #line = "ms_filename=%s.ms" %prefix_name
        elif "oskar_sky_model" in line:
            line = "oskar_sky_model\\file=%s" %sky_osm_name
        elif "healpix_fits" in line and "file" in line:
            if healpix:
                heal_name = healpix + '_%.3fMHz.fits' %(freq / 1.0e+6)
                line = "healpix_fits\\file=%s" %heal_name
        elif "input_directory" in line:
            line = 'input_directory=%s' %telescope_dir
        else:
            pass
        out_file.write(line+'\n')
#    out_file.write('pointing_file=/home/jline/Documents/shintaro_foregrounds/quick_OSKAR/pointing_file.txt\n')
    out_file.close()
    
int_jd, float_jd = calc_jdcal(oskar_date)
print int_jd, float_jd

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
    
#def get_osk_data(oskar_vis_tag=None,polarisation=None):
    #OSK = loadtxt('%s_%s.txt' %(oskar_vis_tag,polarisation))

    #O_us = OSK[:,1]
    #O_vs = OSK[:,2]
    #O_ws = OSK[:,3]
    #O_res = OSK[:,4]
    #O_ims = OSK[:,5]
    
    #return O_us,O_vs,O_ws,O_res,O_ims
    
#@profile
def the_main_loop(tsteps=None):
    ##For each time step
    for time_ind,tstep in enumerate(tsteps):
        time = add_time(oskar_date,tstep)
        ##Precess ra by time since the beginning of the observation 
        ##(convert time to angle, change from seconds to degrees)
        ##Include half of the time step
        
        ##DO NOT ADD ON HALF A TIME STEP - I *think* OSKAR does this internally
        #ra = intial_ra_point + (((tstep + time_int/2.0)*SOLAR2SIDEREAL)*(15/3600.0))
        ra = intial_ra_point + (((tstep)*SOLAR2SIDEREAL)*(15/3600.0))
        if ra >=360.0: ra -= 360.0
        
        ##Need the lst to work out the hour angle if rephasing final data
        if options.phase_centre:
            lst_central = intial_lst + (((tstep + time_int/2.0)*SOLAR2SIDEREAL)*(15/3600.0))
            #lst_central = intial_lst + (((tstep)*SOLAR2SIDEREAL)*(15/3600.0))
            if lst_central >=360.0: lst_central -= 360.0
            ha_phase = lst_central - final_phase_ra
        
        ##Need to know where in the uvfits file structure to put data
        array_time_loc = num_baselines*time_ind

        ##If we are using a sky model with fixed spectral indexes, can run all
        ##frequencies using the same osm model, and so only need to run OSKAR once
        ##per time step
        if options.osm:
            
            ##Make prefix name for the individual time step
            ##If less than second time step, RTS needs a different naming convention
            prefix_name = "%s_band%02d_%.3f" %(outname,band_num,tstep)
            
            ##Start at the first good_chan freq, and do enough chans to cover first good chan to last good chan
            num_channels = good_chans[-1] - good_chans[0]
            oskar_channels = range(good_chans[0],good_chans[0]+num_channels+1)
            oskar_inds = [oskar_channels.index(good_chan) for good_chan in good_chans]
            num_oskar_channels = len(oskar_channels)
            #print "NUM CHANNELS",num_channels
            
            low_channel_freq = base_freq + (good_chans[0]*ch_width)
            
            make_ini(prefix_name=prefix_name,ra=ra,dec=dec_point,freq=low_channel_freq,start_time=time,sky_osm_name=options.osm,healpix=healpix,num_channels=num_oskar_channels)
            
            ##Run the simulation
            cmd = "oskar_sim_interferometer %s.ini" %prefix_name
            run_command(cmd)
            
            ##Read in the data directly from the binary file
            ##This file contains all frequency channels for this time step
            num_vis = num_baselines * num_oskar_channels
            uu,vv,ww,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims = read_oskar_binary(filename="%s.vis" %prefix_name,num_vis=num_vis,num_baselines=num_baselines)
            
            ##Clean up the oskar outputs
            cmd = "rm %s.ini %s.vis" %(prefix_name,prefix_name)
            run_command(cmd)
            
            for oskar_ind,chan in zip(oskar_inds,good_chans):
                
                freq = base_freq + (chan*ch_width)
                ##Use the centre of the fine channel
                freq_cent = freq + (ch_width / 2.0)
                
                osk_low = oskar_ind * num_baselines
                osk_high = (oskar_ind + 1) * num_baselines
                ##Stored in metres in binary, convert to wavelengths
                chan_uu = uu / (VELC / freq)
                chan_vv = vv / (VELC / freq)
                chan_ww = ww / (VELC / freq)
                
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
                
                ##If phasing to new phase centre, undo OSKAR phase tracking, and add in
                ##desired phase tracking by calculating u,v,ws towards phase centre
                if options.phase_centre:
                    wavelength = VELC / freq
                    
                    u_phases, v_phases, w_phases = get_uvw_freq_local(x_length=array(x_lengths),y_length=array(y_lengths),z_length=array(z_lengths),dec=final_phase_dec*D2R,ha=ha_phase*D2R,wavelength=wavelength)
                
                    rotated_xx = rotate_phase(wws=chan_ww,visibilities=comp_xx,new_wws=w_phases)
                    rotated_xy = rotate_phase(wws=chan_ww,visibilities=comp_xy,new_wws=w_phases)
                    rotated_yx = rotate_phase(wws=chan_ww,visibilities=comp_yx,new_wws=w_phases)
                    rotated_yy = rotate_phase(wws=chan_ww,visibilities=comp_yy,new_wws=w_phases)
                    
                    ##Final u,v,w should be toward phase centre
                    chan_uu = u_phases
                    chan_vv = v_phases
                    chan_ww = w_phases
                else:
                    ##Just remove the phase tracking added in by OSKAR
                    rotated_xx = rotate_phase(wws=chan_ww,visibilities=comp_xx)
                    rotated_xy = rotate_phase(wws=chan_ww,visibilities=comp_xy)
                    rotated_yx = rotate_phase(wws=chan_ww,visibilities=comp_yx)
                    rotated_yy = rotate_phase(wws=chan_ww,visibilities=comp_yy)
    
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
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,0,2] = ones(len(uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,1,2] = ones(len(uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,2,2] = ones(len(uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,3,2] = ones(len(uu))
                
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
        ##extrepolation between points, or use a fit of a 2nd order polynomial
        else:
            ###For frequency channel in band
            for chan in good_chans:
                ##Take the band base_freq and add on fine channel freq
                freq = base_freq + (chan*ch_width)
                if time_int < 1:
                    prefix_name = "%s_%.3f_%05.2f" %(outname,freq/1e+6,tstep)
                else:
                    prefix_name = "%s_%.3f_%02d" %(outname,freq/1e+6,int(tstep))
                
                ##Create ini file to run oskar
                sky_osm_name = "%s_%.3f.osm" %(outname,freq/1e+6)
                make_ini(prefix_name=prefix_name,ra=ra,dec=MWA_LAT,freq=freq,start_time=time,sky_osm_name=sky_osm_name,healpix=healpix,num_channels=1)
                ##Run the simulation
                cmd = "oskar_sim_interferometer %s.ini" %prefix_name
                run_command(cmd)
                
                ##Read in the data directly from the binary file
                ##Only one freq channel so number of vis is number of baselines
                uu,vv,ww,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims = read_oskar_binary(filename="%s.vis" %prefix_name,num_vis=num_baselines,num_baselines=num_baselines)
                
                ##Stored in metres in binary, convert to wavelengths
                uu /= (VELC / freq)
                vv /= (VELC / freq)
                ww /= (VELC / freq)
                
                ##Make complex numpy arrays
                comp_xx = make_complex(xx_res,xx_ims)
                comp_xy = make_complex(xy_res,xy_ims)
                comp_yx = make_complex(yx_res,yx_ims)
                comp_yy = make_complex(yy_res,yy_ims)
                
                ##Remove the phase tracking added in by OSKAR
                rotated_xx = rotate_phase(wws=ww,visibilities=comp_xx)
                rotated_xy = rotate_phase(wws=ww,visibilities=comp_xy)
                rotated_yx = rotate_phase(wws=ww,visibilities=comp_yx)
                rotated_yy = rotate_phase(wws=ww,visibilities=comp_yy)
    
                ##Clean up the oskar outputs
                cmd = "rm %s.ini %s.vis" %(prefix_name,prefix_name)
                run_command(cmd)

                ##Use the centre of the fine channel
                freq_cent = freq + (ch_width / 2.0)
                
                oskar_vis_tag = "%s/%s" %(tmp_dir,prefix_name)
                output_uvfits_name = "%s/%s.uvfits" %(data_dir,prefix_name)
                
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
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,0,2] = ones(len(uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,1,2] = ones(len(uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,2,2] = ones(len(uu))
                v_container[array_time_loc:array_time_loc+num_baselines,0,0,0,chan,3,2] = ones(len(uu))
                
                ##Only set the u,v,w to the central frequency channel
                ##Header of uvfits has to match central_freq_chan + 1 (uvfits 1 ordered, python zero ordered)
                if chan == central_freq_chan:
                    central_freq_chan_value = freq_cent
                    ##u,v,w stored in seconds by uvfits files
                    uu /= freq_cent
                    vv /= freq_cent
                    ww /= freq_cent
                    
                    uus[array_time_loc:array_time_loc+num_baselines] = uu
                    vvs[array_time_loc:array_time_loc+num_baselines] = vv
                    wws[array_time_loc:array_time_loc+num_baselines] = ww
                    
                    ##Fill in the baselines using the first time and freq uvfits
                    baselines_array[array_time_loc:array_time_loc+num_baselines] = array(template_data['BASELINE'])
                    
                    ##Fill in the fractional julian date, after adding on the appropriate amount of
                    ##time - /(24*60*60) because julian number is a fraction of a whole day
                    adjust_float_jd_array = float_jd_array + (tstep / (24.0*60.0*60.0))
                    date_array[array_time_loc:array_time_loc+num_baselines] = adjust_float_jd_array
                
    output_uvfits_name = "%s/%s_t%02d_f%.3f_band%02d.uvfits" %(data_dir,outname,time_int,ch_width/1e+6,band_num)
    print "Now creating uvfits file %s" %output_uvfits_name
    create_uvfits(freq_cent=central_freq_chan_value,output_uvfits_name=output_uvfits_name,uu=uus,vv=vvs,ww=wws,baselines_array=baselines_array,date_array=date_array)
    print "Written file %s" %output_uvfits_name

##Whack into a function so we can line_profile it
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
        
chdir(cwd)
template_uvfits.close()
