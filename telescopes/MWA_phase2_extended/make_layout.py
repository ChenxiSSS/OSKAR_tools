'''File that makes an OSKAR telescope. In the top level, need a layout text file that 
describes the positions of each station (tile). Use the antenna_locations.txt file from MWA_Tools
Assuming that all stations are the same, you can just specify one station layout - this is contained
in the dir station. Within station, there must be another layout.txt file describing the layout
of the antennas in the station (a 4 x 4 tile)'''

#import matplotlib.pyplot as plt
import os
from subprocess import call

##Text file containing e,n,h coords
array_layout = 'Phase2_LB_Surveyed_Coordinates.csv'
#anntenna_locs = loadtxt(array_layout,usecols=(2,3,4))

out_file = open('layout.txt','w+')

East = []
North = []
Height = []

lines = [line for line in open(array_layout,'r').read().split('\n') if line != '']
for line in lines:
    a,b,c,h,n,e = line.split(',')
    East.append(float(e))
    North.append(float(n))
    Height.append(float(h))

        
        #make all coords negative (trial and error found this gives the
        #expected u,v,w coords)
    north = -float(n)
    east = -float(e)
    height = -float(h)

    out_file.write('%s %.5f %.5f\n' %(east,north,height))
        
out_file.close()
        
##Spacing between dipoles is 1.1m on the ground mesh.
##Station layout coords are relative to station centre
##Farthest dipole centre from station centre is then 1.65
##Assume flat layout (height = 0)
east_start = -1.65
north_start = -1.65

call('mkdir -p station',shell=True)
out_file = open('./station/layout.txt','w+')

for i in xrange(4):
    for j in xrange(4):
        out_file.write('%.3f %.3f\n' %(east_start,north_start))
        east_start += 1.1
    east_start = -1.65
    north_start += 1.1
    
out_file.close()




