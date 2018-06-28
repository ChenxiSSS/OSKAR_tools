#!/usr/bin/env python
import argparse
import healpy as hp
from numpy import arange,pi

parser = argparse.ArgumentParser(description='Convert a celestial healpix map into an OSKAR .osm file')
parser.add_argument('--healpix', help='Healpix to convert to osm - currently ASSUMES celestial coords')
parser.add_argument('--out', help='Enter output name for osm file')

args = parser.parse_args()

data = hp.read_map(args.healpix)
nside = hp.get_nside(data)

dec,ra = hp.pix2ang(nside,arange(hp.nside2npix(nside)))

dec = dec*(180./pi) - 90.0
ra *= (180/pi)


if '.osm' in args.out:
    outfile = open(args.out,'w+')
else:
    outfile = open(args.out+'.osm','w+')


for ind,flux in enumerate(data):
    outfile.write('%.10f %.10f %.10f\n' %(ra[ind],dec[ind],flux))