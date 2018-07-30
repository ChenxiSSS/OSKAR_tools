#!/usr/bin/env python
import argparse
import healpy as hp
#from numpy import arange,pi,ones,cos,sin,arccos,ndarray,where
from numpy import *
from sys import exit

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

k_B = 1.38064852e-23
VELC = 299792458.

def convert2Jy(hp_array,freq,unit):
    '''Converts healix in mK or K into Jy'''
    nside = hp.get_nside(hp_array)
    solid_angle = hp.nside2pixarea(nside)
    wavelen = VELC / freq
    ##10e-26 because Jy
    S = (2*float(k_B)*hp_array*solid_angle) / (wavelen**2*10e-26) 
    
    if unit == 'mK':
        return S*1e-3
    elif unit == 'K':
        return S
    else:
        exit('Unknown unit type - can only have mK or K for --unit')
        
##Equation 27 in Condon et al 2012
def condon_fixbase(bmax=None,freq=None):
    '''Caclulates confusion noise for particular max baseline bmax (m)
    and frequency freq (Hz)'''
    eightarcsec = (8. / 3600.0) * (pi / 180.)
    reso = (VELC / freq) / bmax
    noise = 1.2 * (freq / 3.02e+9)**-0.7 * (reso / eightarcsec) ** (10.0/3.0)
    return noise * 1e-9

def arcdist(RA1,RA2,Dec1,Dec2):  
    '''Calculates the arcdist between 2 points in deg'''
    dr = pi/180.0
    in1 = (90.0 - Dec1)*dr
    in2 = (90.0 - Dec2)*dr
    RA_d = (RA1 - RA2)*dr
    cosalpha = cos(in1)*cos(in2) + sin(in1)*sin(in2)*cos(RA_d)
    ##Sometimes get floating point errors if two sources at exactly
    ##the same position, so account for this:
    if type(cosalpha) == ndarray:
        cosalpha[where(cosalpha > 1.0)] == 1.0
        cosalpha[where(cosalpha < -1.0 )] == -1.0
    else:
        if cosalpha>1.0: cosalpha = 1.0
        elif cosalpha<-1.0: cosalpha = -1.0
    
    alpha = arccos(cosalpha)
    return alpha/dr
    
def gauss_roll_taper(ra_cent=None,dec_cent=None,ras=None,decs=None,high_cut=None,low_cut=None):
    
    taper = ones(ras.shape)
    ang_dist = arcdist(ra_cent,ras,dec_cent,decs)
    
    sigma = (high_cut+low_cut) / 10.0
    
    taper *= exp( -((ang_dist-low_cut)*(ang_dist-low_cut)) / ( 2.0*sigma*sigma) )
    
    taper[where(ang_dist > high_cut)] = 0.0
    taper[where(ang_dist < low_cut)] = 1.0
    
    ##For low baseline sigma if we want in future
    #taper *= 1.0 - exp( - uv_sq / ( 2.0 * short_baseline_sigma*short_baseline_sigma ) );
    #if ( uv_sq < baseline_min*baseline_min ) taper = 0.0;
    
    return taper


parser = argparse.ArgumentParser(description='Convert a celestial healpix map into an OSKAR .osm file')
parser.add_argument('--healpix', help='Healpix to convert to osm - currently ASSUMES celestial coords')
parser.add_argument('--out', help='Enter output name for osm file')
parser.add_argument('--unit',default='Jy', help='Enter units of healpix map - defaults to Jy. Other accepted values are K, mK')
parser.add_argument('--frequency',default=False, help='If converting from mK or K, or applying a taper, need the frequency of the healpix map')
parser.add_argument('--taper_coords',default=False, help='Include to add gaussian taper to field - enter RA,Dec centre of field here (deg) e.g. --taper_coords=0.0,-27.0')
parser.add_argument('--taper_outer_low',default=7.5, help='Angular distance (deg) from ra0,dec0 where the taper should start - default 2.5')
parser.add_argument('--taper_outer_high',default=10.0, help='Angular distance (deg) from ra0,dec0 where the taper should end - default 10')
parser.add_argument('--use_minimum_flux_as_floor',default=False,action='store_true',
                    help='Add to use the minimum pixel value for taper instead of confusion noise (i.e. if converting a 21cmFAST healpix into an osm)')
parser.add_argument('--max_baseline',default=3e+3, help='Max baseline length with which to calculate confusion noise (m) for tapering')
parser.add_argument('--plot',default=False,action='store_true',help='Make a simple plot of ra/dec to see sky coverage')
parser.add_argument('--auto_zero',default=False,action='store_true',help='Automatically makes the average flux zero - this is important to avoid sharp edges in flux when tapering (especially with MWA sims)')

args = parser.parse_args()

unit = args.unit
low_cut = float(args.taper_outer_low)
high_cut = float(args.taper_outer_high)
frequency = args.frequency

##Check arguments are all good
if unit != 'mK' and unit != 'K' and unit != 'Jy':
    print 'HERE'
    exit('Unknown unit type - can only have Jy, mK or K for --unit')
else:
    if unit == 'mK' or unit == 'K':
        if not frequency:
            exit('Must enter a value for --frequency when converting from mK or K, or when tapering')
            
if args.taper_coords:
    try:
        ra0,dec0 = map(float,args.taper_coords.split(','))
    except:
        exit('Cannot convert --taper_coords into ra,dec - must be two numbers separated by a comma')

##Read in healpix map and get nside
##Convert in Jy if necessary
flux = hp.read_map(args.healpix)

if unit == 'mK' or unit == 'K':
    flux = convert2Jy(flux,float(frequency),unit)

nside = hp.get_nside(flux)
##Find out ra,dec and convert to degrees
dec,ra = hp.pix2ang(nside,arange(hp.nside2npix(nside)))
dec = dec*(180./pi) - 90.0
ra *= (180/pi)

if args.taper_coords:
    ra0,dec0 = map(float,args.taper_coords.split(','))

    #if args.use_minimum_flux_as_floor:
        #noise_floor = flux.min()
        #print 'Noise floor is',noise_floor,flux.max()
    #else:
        #noise_floor = condon_fixbase(bmax=float(args.max_baseline),freq=float(frequency))
        
    #if noise_floor <= 0:
        #adjust_noise = abs(noise_floor)
        ##adjust_noise = 1.0
        
    new_flux = flux * gauss_roll_taper(ra_cent=ra0,dec_cent=dec0,ras=ra,decs=dec,high_cut=high_cut,low_cut=low_cut)
    taper = where(new_flux != 0)
        
    
    #fig = plt.figure(figsize=(10,10))
    #ax = fig.add_subplot(111)
    
    #ras = linspace(45,75,1000)
    #decs = ones(1000)*-26.7
    
    #ax.plot(ras,gauss_roll_taper(ra_cent=ra0,dec_cent=dec0,ras=ras,decs=decs,high_cut=high_cut,low_cut=low_cut))
    
    #fig.savefig('taper_used_RAslice.png',bbox_inches='tight')
    
        
    ##else:
        ##cut = noise_floor * gauss_roll_taper(ra_cent=ra0,dec_cent=dec0,ras=ra,decs=dec,high_cut=high_cut,low_cut=low_cut)
        ##taper = where(flux - cut > 0)
        
    ra = ra[taper]
    dec = dec[taper]
    
    if args.auto_zero:
        flux = new_flux[taper] - mean(new_flux[taper])
    else:
        flux = new_flux[taper]
    
print('There are %d sources in the tapered osm' %len(ra))

if '.osm' in args.out:
    outfile = open(args.out,'w+')
else:
    outfile = open(args.out+'.osm','w+')

for ind,S in enumerate(flux):
    outfile.write('%.10f %.10f %.10f\n' %(ra[ind],dec[ind],S))
    
    
if args.plot:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.scatter(ra,dec,c=flux,s=1)
    
    fig.savefig('%s.png' %args.out.strip('.osm'),bbox_inches='tight')


