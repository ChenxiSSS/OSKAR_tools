#!/usr/bin/env python
from subprocess import call
from sys import exit
from optparse import OptionParser
from numpy import zeros, pi, sin, cos, real, imag, loadtxt, array, floor, arange, ones, where, ceil
from numpy import exp as n_exp
from ephem import Observer
from cmath import exp
from sys import path as sys_path
try:
    from jdcal import gcal2jd
except ImportError:
    sys_path.append('/lustre/projects/p048_astro/jline/software/jdcal-1.3')
    from jdcal import gcal2jd
    
from os import environ,getcwd,chdir,makedirs,path,walk
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


def get_uvw_zenith(x_length=None,y_length=None,z_length=None):
    '''Takes the baseline length in meters and returns the u,v,w at zenith'''
    ha = 0.0
    dec = MWA_LAT*D2R
    
    u = sin(ha)*x_length + cos(ha)*y_length
    v = sin(dec)*cos(ha)*x_length + sin(dec)*sin(ha)*y_length + cos(dec)*z_length
    w = cos(dec)*cos(ha)*x_length - cos(dec)*sin(ha)*y_length + sin(dec)*z_length
    
    return u,v,w

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
    
    #rotated_visis = visibilities
    
    return rotated_visis

def get_uvw_data(d_single=None,d_double=None,num_baselines=None,block=None,num_time_steps=None):
    '''Takes the format of the OSKAR binary data, and pulls out the u, v, or w data
    Returns the data in metres'''
    uvw_data = zeros(num_baselines*num_time_steps)
    
    if d_single:
        ##single precision, number is 4 bytes long
        uvw_inds = arange(0,4*num_baselines*num_time_steps,4)
    elif d_double:
        ##double precision, number is 8 bytes long
        uvw_inds = arange(0,8*num_baselines*num_time_steps,8)
    else:
        print "J Line hasn't coded to deal with the data type in the OSKAR binary. Have coded in single and double precision (float and double) but not int or char. Exiting"
        exit(0)
    
    for ind in xrange(num_baselines*num_time_steps):
        if d_single:
            uvw_data[ind] = unpack('<f',block[uvw_inds[ind]:uvw_inds[ind]+4])[0]
        elif d_double:
            uvw_data[ind] = unpack('<d',block[uvw_inds[ind]:uvw_inds[ind]+8])[0]
    
    return uvw_data


def recover_block(tag_index=None,binary_file=None,return_data_type=False):
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
    
    if return_data_type:
        return block,d_single,d_double
    else:
        return block

def read_oskar_binary(filename=None,num_time_steps=None,num_baselines=None,num_channels=None):
    '''Open the native OSKAR binary format and extracts the visibilities.
    Returns the u,v,w coords and complex XX,XY,YX,YY arrays
    Will ONLY work with the default OSKAR settings found in the template .ini
    files. TODOs are needs to find the data element type (e.g. float or double),
    the dimensions of the data etc. Currently hard-coded to defaults'''
    
    ##The binary format is detailed here: http://oskar.oerc.ox.ac.uk/sites/default/files/2.7/doc/OSKAR-2.6-Binary-File-Format.pdf
    
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
        
    tag_indexes = array(tag_indexes).astype(int)

    ##emtpy arrays a needed values
    #wavelength = VELC / freq
    uus = zeros(num_time_steps*num_baselines)
    vvs = zeros(num_time_steps*num_baselines)
    wws = zeros(num_time_steps*num_baselines)
    
    xx_res = zeros(num_time_steps*num_baselines*num_channels)
    xx_ims = zeros(num_time_steps*num_baselines*num_channels)
    xy_res = zeros(num_time_steps*num_baselines*num_channels)
    xy_ims = zeros(num_time_steps*num_baselines*num_channels)
    yx_res = zeros(num_time_steps*num_baselines*num_channels)
    yx_ims = zeros(num_time_steps*num_baselines*num_channels)
    yy_res = zeros(num_time_steps*num_baselines*num_channels)
    yy_ims = zeros(num_time_steps*num_baselines*num_channels)
    
    group_IDs = []
    tag_IDs = []
    crc_codes = []
    
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
            
        group_IDs.append(group_ID)
        tag_IDs.append(tag_ID)

    group_IDs = array(group_IDs)
    tag_IDs = array(tag_IDs)
    crc_codes = array(crc_codes)
    
    
    ##These two number help calculate how many visibility blocks there are
    max_num_time_samps_ind = where((group_IDs == 11) & (tag_IDs == 7))
    block_max_num_time_samps = recover_block(tag_index=tag_indexes[max_num_time_samps_ind][0],binary_file=binary_file)
    max_num_time_samps = unpack('<I',block_max_num_time_samps)[0]
    
    
    total_num_usable_time_samps_ind = where((group_IDs == 11) & (tag_IDs == 8))
    block_total_num_usable_time_samps = recover_block(tag_index=tag_indexes[total_num_usable_time_samps_ind][0],binary_file=binary_file)
    total_num_usable_time_samps = unpack('<I',block_total_num_usable_time_samps)[0]
    
    ##The actual number of visi blocks goddam
    num_visi_blocks = ceil(float(total_num_usable_time_samps) / float(max_num_time_samps))

    ##Where the visibility block data heads live
    dimension_start_and_size_inds = where((group_IDs == 12) & (tag_IDs == 1))[0]
    
    ##Do a sanity check and see whether the number of visibility headers we can find
    ##matches the predicted number of visibility data block
    if len(dimension_start_and_size_inds) != num_visi_blocks:
        print 'WARNING - len(dimension_start_and_size_inds) != num_visi_blocks\nUnclear how many blocks of visibility data are in the binary file'
    
    time_steps_per_block = []
    
    for dim_ind in dimension_start_and_size_inds:
        block = recover_block(tag_index=tag_indexes[dim_ind],binary_file=binary_file)
        
        ##From OSKAR manual, data on each of the visibility blocks
        #Global start time index for the first time in the block, relative to observation start.
        #Global start channel index for the first channel in the block (usually 0; reserved for future use).
        #Number of usable time samples in the block.
        #Number of usable frequency channels in the block (usually the total number; reserved for future #use).
        #Number of cross-correlated baselines.
        #Number of stations.
        
        glob_start_time = unpack('<I',block[:4])[0]
        start_channel_index = unpack('<I',block[4:8])[0]
        usable_time_samps = unpack('<I',block[8:12])[0]
        usable_freq_chans = unpack('<I',block[12:16])[0]
        num_cross_baselines = unpack('<I',block[16:20])[0]
        num_stations  = unpack('<I',block[20:24])[0]
        
        time_steps_per_block.append(usable_time_samps)
        
        
    ##Where the visibility block data lives (XX,XY,YX,YY)
    visi_data_block_inds = where((group_IDs == 12) & (tag_IDs == 3))[0]
    
    ##Use this number to work out where to insert the visibilities
    ##into the blank visi containers
    filled_visis = 0

    ##For each visibility data block, use the expected number of time steps
    ##to pull out the correct number of visibilities
    for visi_ind,block_time_steps in zip(visi_data_block_inds,time_steps_per_block):
        block,d_single,d_double = recover_block(tag_index=tag_indexes[visi_ind],binary_file=binary_file,return_data_type=True)
        
        num_vis = block_time_steps*num_baselines*num_channels
    
        ###*8 because we 4 complex numbers so 4 real 4 imag
        all_nums = zeros(8*num_vis)
        if d_single:
            ##*4 because floats are 4 bytes long
            visi_inds = arange(0,8*4*num_vis,4)
        elif d_double:
            ##*8 because doubles are 8 bytes long
            visi_inds = arange(0,8*8*num_vis,8)
        
        for ind in xrange(8*num_vis):
            if d_single:
                all_nums[ind] = unpack('<f',block[visi_inds[ind]:visi_inds[ind]+4])[0]
            elif d_double:
                all_nums[ind] = unpack('<f',block[visi_inds[ind]:visi_inds[ind]+8])[0]
        
        ###Split the data up by polarisation
        
        scale = 1.0
        
        xx_res[filled_visis:filled_visis+num_vis] = all_nums[arange(0,8*num_vis,8)]*scale
        xx_ims[filled_visis:filled_visis+num_vis] = all_nums[arange(1,8*num_vis,8)]*scale
        xy_res[filled_visis:filled_visis+num_vis] = all_nums[arange(2,8*num_vis,8)]*scale
        xy_ims[filled_visis:filled_visis+num_vis] = all_nums[arange(3,8*num_vis,8)]*scale
        yx_res[filled_visis:filled_visis+num_vis] = all_nums[arange(4,8*num_vis,8)]*scale
        yx_ims[filled_visis:filled_visis+num_vis] = all_nums[arange(5,8*num_vis,8)]*scale
        yy_res[filled_visis:filled_visis+num_vis] = all_nums[arange(6,8*num_vis,8)]*scale
        yy_ims[filled_visis:filled_visis+num_vis] = all_nums[arange(7,8*num_vis,8)]*scale
        
        filled_visis += num_vis
            
        
    ##In the OSKAR binary format, u,v,w stored under group_ID = 12, as tag_ID 4,5,6
    ##respectively. Stored in metres so divide by wavelength for natural units
    ##Each number is 4 bytes long
    uu_data_block_inds = where((group_IDs == 12) & (tag_IDs == 4))[0]
    vv_data_block_inds = where((group_IDs == 12) & (tag_IDs == 5))[0]
    ww_data_block_inds = where((group_IDs == 12) & (tag_IDs == 6))[0]
    
    
    
    def fill_uvw_data(tag_index=None,num_time_steps=None,filled_uvws=None,uvw_container=None):
        block,d_single,d_double = recover_block(tag_index=tag_indexes[tag_index],binary_file=binary_file,return_data_type=True)
        uvw = get_uvw_data(d_single=d_single,d_double=d_double,num_baselines=num_baselines,block=block,num_time_steps=num_time_steps)
        uvw_container[filled_uvws:filled_uvws+(num_baselines*num_time_steps)] = uvw
        
    filled_uvws = 0
    for block_ind in xrange(int(num_visi_blocks)):
        fill_uvw_data(tag_index=uu_data_block_inds[block_ind],num_time_steps=time_steps_per_block[block_ind],filled_uvws=filled_uvws,uvw_container=uus)
        fill_uvw_data(tag_index=vv_data_block_inds[block_ind],num_time_steps=time_steps_per_block[block_ind],filled_uvws=filled_uvws,uvw_container=vvs)
        fill_uvw_data(tag_index=ww_data_block_inds[block_ind],num_time_steps=time_steps_per_block[block_ind],filled_uvws=filled_uvws,uvw_container=wws)
    
        filled_uvws += (num_baselines*time_steps_per_block[block_ind])
        
    return uus,vvs,wws,xx_res,xx_ims,xy_res,xy_ims,yx_res,yx_ims,yy_res,yy_ims
        
def make_complex(re=None,im=None):
    '''Takes two arrays, and returns a complex array with re real values and im imarginary values'''
    comp = array(re,dtype=complex)
    comp += 1j * im
    
    return comp

def create_uvfits(v_container=None,freq_cent=None, ra_point=None, dec_point=None, 
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
    
def make_ini(prefix_name=None,ra=None,dec=None,freq=None,start_time=None,sky_osm_name=None,healpix=None,
            num_channels=None,template_ini=None,ch_width=None,time_int=None,telescope_dir=None,
            num_time_steps=None,obs_time_length=None,oskar_gsm=False,oskar_gsm_file=None,
            oskar_gsm_SI=None,save_ms=False):
    out_file = open("%s.ini" %prefix_name,'w+')
    for line in template_ini:
        if "num_channels" in line:
            line = "num_channels=%d" %num_channels
        elif "start_frequency_hz" in line:
            line = "start_frequency_hz=%.10f" %freq
        elif "frequency_inc_hz" in line:
            line = "frequency_inc_hz=%.10f" %ch_width
        elif "phase_centre_ra_deg" in line:
            line = "phase_centre_ra_deg=%.10f" %ra
        elif "phase_centre_dec_deg" in line:
            line = "phase_centre_dec_deg=%.10f" %dec
        elif "start_time_utc" in line:
            line = "start_time_utc=%s" %start_time
        elif line.split('=')[0]=="length":
            line = "length=%.10f" %obs_time_length
        elif "oskar_vis_filename" in line:
            line = "oskar_vis_filename=%s.vis" %prefix_name
        elif "channel_bandwidth_hz" in line:
            line = "channel_bandwidth_hz=%.10f" %ch_width
        elif "time_average_sec" in line:
           line = "time_average_sec=%.10f" %time_int
        elif "oskar_sky_model" in line:
            line = "oskar_sky_model\\file=%s" %sky_osm_name
        elif "input_directory" in line:
            line = 'input_directory=%s' %telescope_dir
        elif "num_time_steps" in line:
            line = "num_time_steps=%d" %num_time_steps
            
        ##Inspect the optional lines, and include if necessary
        elif '#' in line:
            
            if 'gsm' in line:
                if oskar_gsm:
                    if 'gsm/file' in line:
                        line = 'gsm/file=%s' %oskar_gsm_file
                    elif 'gsm/spectral_index' in line:
                        line = 'gsm/spectral_index=%.5f' %oskar_gsm_SI
                else:
                    line = 'pass'
                    
            elif 'ms_filename' in line:
                if save_ms:
                    line = 'ms_filename=%s.ms' %prefix_name
                else:
                    line = 'pass'
                    
            #elif "healpix_fits" in line and "file" in line:
                #if healpix:
                    #heal_name = healpix + '_%.3fMHz.fits' %(freq / 1.0e+6)
                    #line = "healpix_fits\\file=%s" %heal_name
            else:
                line = 'pass'
            
        else:
            pass
            
        if line == 'pass':
            pass
        else:
            out_file.write(line+'\n')
    
    out_file.close()
    
def create_flagged_telescope(metafits=None,antenna_coord_file=None,azimuth=None,altitude=None,telescope_dir=None):
    '''Takes an MWA layout text file with columns "label east north height" and uses them
    to create an OSKAR telescope model. Uses the metafits file to flag out the dipoles
    based on the XX flags - OSKAR can only flag whole elements, rather than individual pols'''
    
    hdu = fits.open(metafits)

    tiles =  hdu[1].data['Tile']
    dipole_delays = hdu[1].data['Delays']

    hdu.close()

    ##Order of delays is XX1,YY1,XX2,YY2 etc so select just XX
    XX = arange(0,256,2).astype(int)
    
    ##If a delay is set to 32, the dipole is flagged
    tile_inds,dipole_flags = where(dipole_delays[XX,:] == 32)
    
    tile_flags = tiles[tile_inds]
    
    ##Make flags direction
    flags_dict = {}
    ##Set up lists to contain flags for tiles
    for tile in set(tile_flags):
        flags_dict['Tile%d' %int(tile)] = []
    
    ##Populate flags. Allows for more than one flag per tile
    for tile,dipole in zip(tile_flags,dipole_flags):
        flags_dict['Tile%d' %int(tile)].append(int(dipole))
    
    
    
    ##Open the text file containing the MWA tile locations, and assign column names and type to the data
    tile_locs = loadtxt(antenna_coord_file, dtype={'names': ('label', 'east', 'north', 'height'), 'formats': ('S7', float, float, float)})

    ##Spacing between dipoles is 1.1m on the ground mesh.
    ##Station layout coords are relative to station centre
    ##Farthest dipole centre from station centre is then 1.65

    ##Set up a dictionary containg labels for the dipoles and their position within the tile
    station_dict = {}

    ##Doing it this way, if below is in east along x axis, north along y axis
    ##  12  13  14  15
    ##  8   9   10  11
    ##  4   5   6   7
    ##  0   1   2   3
    ##From the beam maps I have made, I *think* this is the ordering of the
    ##dipoles in the FEE. I'm going to assume it's the same in the metafits flagging

    east_start = -1.65
    north_start = -1.65
    for i in xrange(4):
        for j in xrange(4):
            ##coo-ords written out as east, north
            station_dict['%d' %(i*4 + j)] = [east_start,north_start]
            east_start += 1.1
        east_start = -1.65
        north_start += 1.1
        
    ##Put telescope model into named telescope_dir - create
    ##telescope_dir if necessary
    if not path.exists(telescope_dir):
        makedirs(telescope_dir)
    else:
        print('Telescope %s already exists: overwriting contents!!' %telescope_dir)
        
    ##Create overall tile layout file - specifies locations of tiles

    out_file = open('%s/layout.txt' %telescope_dir,'w+')

    for station in xrange(len(tile_locs['label'])):
        ##Write out the telescope layout.txt file in
        ##east,north,height coords.

        ##There is some kind of co-ordinate difference between the RTS and OSKAR, so
        ##need to make all coords negative here
        ## (trial and error found this gives the expected u,v,w coords)
        north = -tile_locs['north'][station]
        east = -tile_locs['east'][station]
        height = -tile_locs['height'][station]
        
        out_file.write('%s %.5f %.5f\n' %(east,north,height))
            
    out_file.close()
    
    ##Add the permitted_beams file so beam only points where we care
    permitted = open('%s/permitted_beams.txt' %(telescope_dir),'w+')
    permitted.write('%.5f %.5f' %(azimuth,altitude))
    permitted.close
    
    ##Tell OSKAR where the telescope lives.
    ##TODO HARDCODED TO MRO AT THE MOMENT
    position = open('%s/position.txt' %(telescope_dir),'w+')
    position.write('116.670813889 -26.703319')
    position.close
        
    ##Build indidual station layouts
    for station in range(1,len(tile_locs['label'])+1):
        if not path.exists('%s/station%03d' %(telescope_dir,station)):
            makedirs('%s/station%03d' %(telescope_dir,station))
        
        ##Just in case add the permitted_beams file so beam only points where we care as well
        permitted = open('%s/station%03d/permitted_beams.txt' %(telescope_dir,station),'w+')
        permitted.write('%.5f %.5f' %(azimuth,altitude))
        permitted.close
        
        out_file = open('%s/station%03d/layout.txt' %(telescope_dir,station),'w+')
        
        if tile_locs['label'][station-1] in flags_dict:
            print('Flagging %s (station%03d)' %(tile_locs['label'][station-1],station))
            flags = flags_dict[tile_locs['label'][station-1]]
            for dipole in xrange(16):
                if dipole not in flags:
                    east,north = station_dict['%d' %dipole]
                    out_file.write('%.3f %.3f\n' %(east,north))
        else:
            for dipole in xrange(16):
                east,north = station_dict['%d' %dipole]
                out_file.write('%.3f %.3f\n' %(east,north))
        
        out_file.close()
        
def add_data_to_uvfits(v_container=None,time_ind=None,num_baselines=None,template_baselines=None,
    chan=None,num_oskar_channels='False',oskar_ind=False,uu=None,vv=None,ww=None,float_jd_array=None,
    baselines_array=None,date_array=None,undo_phase_track=True,xx_res=None,xx_ims=None,
    xy_res=None,xy_ims=None,yx_res=None,yx_ims=None,yy_res=None,yy_ims=None,x_lengths=None,y_lengths=None,
    z_lengths=None,uus=None,vvs=None,wws=None,tstep=None,freq=None,ch_width=None,central_freq_chan=None,
    chips_settings=False):
    
    time_ind_lower = time_ind*num_baselines
    
    ##Stored in metres in OSKAR binary, convert to wavelengths
    chan_uu = uu[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
    chan_vv = vv[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
    chan_ww = ww[time_ind_lower:time_ind_lower+num_baselines] / (VELC / freq)
    
    ##Depending on if OSKAR was run with chans separately or together, select channels appropriately
    if oskar_ind == 'False':
        osk_low = time_ind_lower
        osk_high = time_ind_lower+num_baselines
        
    else:
        
        osk_low = (num_oskar_channels*time_ind*num_baselines) + (oskar_ind*num_baselines)
        osk_high = osk_low + num_baselines
        
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
    
    ##Remove the phase tracking added in by OSKAR if requested
    if undo_phase_track:
        final_xx = rotate_phase(wws=chan_ww,visibilities=comp_xx)
        final_xy = rotate_phase(wws=chan_ww,visibilities=comp_xy)
        final_yx = rotate_phase(wws=chan_ww,visibilities=comp_yx)
        final_yy = rotate_phase(wws=chan_ww,visibilities=comp_yy)
    else:
        final_xx = comp_xx
        final_xy = comp_xy
        final_yx = comp_yx
        final_yy = comp_yy
    
    #print oskar_ind,chan_ww[112],freq,comp_xx[112],rotated_xx[112]

    ##Use the centre of the fine channel
    freq_cent = freq + (ch_width / 2.0)
    
    ##If doing a mock CHIPS obs, the central channel is an average
    ##of a full and an empty channel, so need to set weights to 0.5
    if chips_settings:
        if chan == central_freq_chan:
            scale = 0.5
        else:
            scale = 1.0
    else:
        scale = 1.0
    
    ##Populate the v_container at the correct spots
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,0,0] = real(final_xx)*scale
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,0,1] = imag(final_xx)*scale
    
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,1,0] = real(final_yy)*scale
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,1,1] = imag(final_yy)*scale
    
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,2,0] = real(final_xy)*scale
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,2,1] = imag(final_xy)*scale
    
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,3,0] = real(final_yx)*scale
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,3,1] = imag(final_yx)*scale
    
    ##Set the weights of everything to ones
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,0,2] = ones(len(chan_uu))*scale
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,1,2] = ones(len(chan_uu))*scale
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,2,2] = ones(len(chan_uu))*scale
    v_container[time_ind_lower:time_ind_lower+num_baselines,0,0,0,chan,3,2] = ones(len(chan_uu))*scale
    
    ##Only set the u,v,w to the central frequency channel
    ##Header of uvfits has to match central_freq_chan + 1 (uvfits 1 ordered, python zero ordered)
    if chan == central_freq_chan:
        central_freq_chan_value = freq_cent
        
        ##u,v,w stored in seconds by uvfits files
        ##if removing phase tracking, set 0 coords of u,v,w to zenith
        if undo_phase_track:
            final_uu,final_vv,final_ww = get_uvw_zenith(x_length=x_lengths,y_length=y_lengths,z_length=z_lengths)
        else:
            final_uu = chan_uu / freq_cent
            final_vv = chan_vv / freq_cent
            final_ww = chan_ww / freq_cent
        
        uus[time_ind_lower:time_ind_lower+num_baselines] = final_uu
        vvs[time_ind_lower:time_ind_lower+num_baselines] = final_vv
        wws[time_ind_lower:time_ind_lower+num_baselines] = final_ww
        
        ##Fill in the baselines using the first time and freq uvfits
        baselines_array[time_ind_lower:time_ind_lower+num_baselines] = template_baselines
        
        ##Fill in the fractional julian date, after adding on the appropriate amount of
        ##time - /(24*60*60) because julian number is a fraction of a whole day
        adjust_float_jd_array = float_jd_array + (float(tstep) / (24.0*60.0*60.0))
        
        date_array[time_ind_lower:time_ind_lower+num_baselines] = adjust_float_jd_array
        