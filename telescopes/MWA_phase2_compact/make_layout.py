'''File that makes an OSKAR telescope. In the top level, need a layout text file that 
describes the positions of each station (tile). Use the antenna_locations.txt file from MWA_Tools
Assuming that all stations are the same, you can just specify one station layout - this is contained
in the dir station. Within station, there must be another layout.txt file describing the layout
of the antennas in the station (a 4 x 4 tile)'''

#import matplotlib.pyplot as plt
import os

##Use the antenna locations given in MWA_Tools
template = open('HexMeasured_phaseIIEOR.csv').read().split('\n')
#template = open('Phase2_EOR_Proposed_Coords.csv').read().split('\n')
out_file = open('layout.txt','w+')

East = []
North = []

for line in template:
    if '#' in line or line=='':
        pass
    else:
        ##Write out the telescope layout.txt file in
        ##east,north,height coords.
        name,east,north,height = line.split()
#        name,lon,lat,east,north,height = line.split(',')
        East.append(east)
        North.append(north)
        
        #make all coords negative (trial and error found this gives the
        #expected u,v,w coords)
        north = -float(north)
        east = -float(east)
        height = -float(height)
        
        out_file.write('%s %.5f %.5f\n' %(east,north,height))
        
out_file.close()
        
##Spacing between dipoles is 1.1m on the ground mesh.
##Station layout coords are relative to station centre
##Farthest dipole centre from station centre is then 1.65
##Assume flat layout (height = 0)
east_start = -1.65
north_start = -1.65

#os.mkdir('station')
out_file = open('./station/layout.txt','w+')

for i in xrange(4):
    for j in xrange(4):
        out_file.write('%.3f %.3f\n' %(east_start,north_start))
        east_start += 1.1
    east_start = -1.65
    north_start += 1.1
    
out_file.close()




