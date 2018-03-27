from numpy import *
import matplotlib.pyplot as plt
from glob import glob
from subprocess import call


dirs = glob('/home/jline/Dropbox/work/SKA_visit/SKA_all_stations.tm/*')

for dir in dirs:
    if 'txt' in dir:
        pass
    else:
        station = dir.split('/')[-1]
        if len(station.split('_')) > 1:
            front,back = station.split('_')
            letter = front[0]
            number = int(front[1:])
            back = int(back)
            
            new_station = '%s%03d_%03d' %(letter,number,back)
            
        else:
            new_station = station
        
        
    cmd = 'cp -r %s /home/jline/software/OSKAR_tools/telescopes/SKA/%s' %(dir,new_station)
    call(cmd,shell=True)