#!/usr/bin/env python
from subprocess import call
from sys import exit
from optparse import OptionParser
from numpy import zeros, pi, sin, cos, real, imag, loadtxt, array, floor, arange, ones, where
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
    
from astropy.io import fits

R2D = 180.0 / pi
D2R = pi / 180.0
MWA_LAT = -26.7033194444
VELC = 299792458.0

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
    '''Takes a string format date and returns julian date in
    the two chunks a uvfits file likes'''
    
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
def rotate_phase(wws=None,visibilities=None):
    '''Undoes any phase tracking applied to data - to phase track, a phase was applied
    to counter the delay term caused by w term of baseline - so just apply the opposite
    w term.'''

    sign = 1
    PhaseConst = 1j * 2 * pi * sign
    
    ##theory - so normal phase delay is caused by path difference across
    ##a base line, which is u*l + v*m + w*n
    ##To phase track, you insert a phase to make sure there is no w contribution at
    ##phase centre; this is when n = 1, so you insert a phase thus:
    ##a base line, which is u*l + v*m + w*(n - 1)
    ##So we just need to remove the effect of the -w term
    phase_rotate = n_exp( PhaseConst * wws)
    rotated_visis = visibilities * phase_rotate
    
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

def create_uvfits(v_container=None,freq_cent=None, ra_point=None, dec_point=MWA_LAT, 
                  output_uvfits_name=None,uu=None,vv=None,ww=None,
                  baselines_array=None,date_array=None,date=None,
                  central_freq_chan=None,ch_width=None,template_uvfits=None,
                  int_jd=None):
    '''Takes visibility date and writes out a uvfits format file'''
    
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
    hdulist.writeto(output_uvfits_name,overwrite=True)
    hdulist.close()
    
def make_ini(prefix_name=None,ra=None,dec=None,freq=None,start_time=None,sky_osm_name=None,healpix=None,num_channels=None,template_ini=None,ch_width=None,time_int=None,telescope_dir=None):
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
            line = "phase_centre_dec_deg=%.10f" %MWA_LAT
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
    out_file.close()