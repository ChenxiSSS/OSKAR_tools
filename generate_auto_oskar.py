#!/usr/bin/env python
from optparse import OptionParser
from sys import exit
import os
from numpy import ceil,floor

OSKAR_dir = os.environ['OSKAR_TOOLS']

def write_oskar(wd=None, metafits=None, srclist=None, point_source_tag=None, time=None, job_bands=None, data_dir=None, telescope=None, time_int=None, ini_file=None, jobs_per_GPU=None):
    '''Writes a bash script for each course band to run OSKAR'''
    
    start, finish = map(float,time.split(','))
    num_time_steps = int((finish - start) / float(time_int))
    
    ##Takes around 7.5 mins to do one time step with 27 channels
    hours = ceil((num_time_steps * 7.5) / 60.0)
    
    ##Set up controlling qsub script to run current batch of OSKAR processes
    file_name = 'qsub_%s_bands%02d-%02d_t%d-%d.sh' %(point_source_tag,job_bands[0],job_bands[-1],start,finish)
    out_file = open(file_name,'w+')
    out_file.write('#!/bin/bash\n')
    out_file.write('#PBS -l nodes=1:gpus=2:ppn=4\n')
    out_file.write('#PBS -l walltime=%02d:00:00\n' %hours)
    out_file.write('#PBS -m e\n')
    out_file.write('#PBS -q gstar\n')
    out_file.write('#PBS -A p048_astro\n')

    #out_file.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
    out_file.write('source /home/jline/.bash_profile\n')
    out_file.write('cd %s\n' %wd)
    
    
    half_of_jobs = len(job_bands) / 2

    ##Setup half of jobs to run on one GPU, half on the other
    run1_name = 'run_%s_bands%02d-%02d_t%d-%d.sh' %(point_source_tag,job_bands[0],job_bands[half_of_jobs-1],start,finish)
    run1 = open(run1_name,'w+')
    #run1.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
    
    run2_name = 'run_%s_bands%02d-%02d_t%d-%d.sh' %(point_source_tag,job_bands[half_of_jobs],job_bands[-1],start,finish)
    run2 = open(run2_name,'w+')
    
    
    def add_modules(filename=None):
        #filename.write('/home/jline/.bash_profile\n')
        filename.write('module load python\n')
        filename.write('module load gcc\n')
        filename.write('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-5.5/lib64\n')
        filename.write('export PATH=/lustre/projects/p048_astro/jline/software/OSKAR-2.7.0-Source/build/gstar_install/bin/:$PATH\n')
        filename.write('export LD_LIBRARY_PATH=/lustre/projects/p048_astro/jline/software/OSKAR-2.7.0-Source/build/gstar_install/lib:$LD_LIBRARY_PATH\n')
        
    add_modules(filename=run1)
    add_modules(filename=run2)
    
    def write_oskar_command(band_num=None,runfile=None):
        oskar_options = "--metafits=%s --srclist=%s --output_name=%s --time=%s --band_num=%s --debug --data_dir=%s" %(metafits, srclist, point_source_tag, time, band_num, data_dir)
        if telescope:
            oskar_options += ' --telescope=%s' %telescope
        if time_int:
            oskar_options += ' --time_int=%s' %time_int
        if ini_file:
            oskar_options += ' --ini_file=%s' %ini_file
        
        runfile.write('time %s/MWAobs_oskar_ascii.py %s &\n' %(OSKAR_dir, oskar_options))
        
    for band_num in job_bands[:half_of_jobs]:
        write_oskar_command(band_num=band_num,runfile=run1)
        
    for band_num in job_bands[half_of_jobs:]:
        write_oskar_command(band_num=band_num,runfile=run2)
        
    run1.close()
    run2.close()
    
    ##Do the fancy CUDA mpi control thing and actually run all the jobs
    out_file.write('time source %s/quick_oskar.sh %s %s\n' %(OSKAR_dir,run1_name,run2_name))
    out_file.write('wait\n')
    out_file.write('rm %s %s %s' %(file_name,run1_name,run2_name))
    out_file.close()
    
    return file_name

#def write_add_diffuse(wd=None, metafits=None, diffuse_tag=None, time=None, band_num=None, data_dir=None, telescope=None, time_int=None, base_uvfits=None):
    #''' '''
    
    #start, finish = map(float,time.split(','))
    #num_time_steps = int((finish - start) / float(time_int))
    
    ###Takes around 10 mins to do one time step with 27 channels
    #hours = ceil((num_time_steps * 15) / 60.0)
    
    #file_name = 'qsub_%s_band%02d_t%d-%d.sh' %(diffuse_tag,band_num,start,finish)
    #out_file = open(file_name,'w+')
    #out_file.write('#!/bin/bash\n')
    #out_file.write('#PBS -l nodes=1\n')
    #out_file.write('#PBS -l walltime=%02d:00:00\n' %hours)
    #out_file.write('#PBS -m e\n')
    #out_file.write('#PBS -q sstar\n')
    #out_file.write('#PBS -A p048_astro\n')

    #out_file.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
    #out_file.write('source /home/jline/.bash_profile\n')
    #out_file.write('cd %s\n' %wd)
    
    #majick_options = "--metafits=%s --output_name=%s --time=%s --band_num=%s --debug --data_loc=%s --base_uvfits=%s" %(metafits, diffuse_tag, time, band_num, data_dir, base_uvfits)
    
    #if telescope:
        #majick_options += ' --telescope=%s' %telescope
    #if time_int:
        #majick_options += ' --twosec=%s' %time_int
    
    #out_file.write('time python $OSKAR_TOOLS/MWAobs_add_MAJICK_diffuse.py ' + majick_options + '\n')
    
    #out_file.close()
    #return file_name
    

parser = OptionParser()

parser.add_option('-o', '--output_dir', default=False, help='Enter output data directory')
parser.add_option('-m', '--metafits', default=False, help='Enter metafits file to base observation on')
parser.add_option('-t', '--time', default=False, help='Enter start,finish times relative to metafits date (seconds - i.e. 0 to start at beginnning')
parser.add_option('-c', '--time_int',default=False, help='Enter time_int of correlator to simulate (s) to override what is in metafits (i.e --time_int=2). Defaults to what is in the metafits')
parser.add_option('-s','--srclist', default=False, help='Enter RTS srclist to use as sky model')
parser.add_option('-p','--point_source_tag', default=False, help='Enter uvfits tag for the point source only model')
parser.add_option('-d','--diffuse_tag', default=False, help='Enter uvfits tag for the point source + diffuse models')
parser.add_option('-l','--diffuse_output_dir', default=False, help='Enter output data directory for diffuse - defaults to same location as point source model')
parser.add_option('-b', '--band_nums', default='all', help='Defaults to running all 24 course bands. Alternatively, enter required numbers delineated by commas, e.g. --band_nums=1,7,9')
parser.add_option('-a','--telescope', default=False, help='Enter telescope tag to use - default is MWA_phase1')
parser.add_option('-i', '--ini_file', default=False, help='Enter template oskar .ini - defaults to the template .ini located in $OSKAR_TOOLS/telescopes/--telescope')
parser.add_option('-j', '--jobs_per_GPU', default=4, help='How many OSKAR jobs to run per GPU - default is 4')

options, args = parser.parse_args()

def false_test(option,name):
    if option == False:
        print '-----------------------------------------------------'
        print '%s not entered but is required. Exiting now!!' %name
        print '-----------------------------------------------------'
        exit()
    else:
        return option

##Get inputs
time = false_test(options.time,'"time"')
metafits = false_test(options.metafits,'"metafits"')
output_dir = false_test(options.output_dir,'"output_dir"')
srclist = false_test(options.srclist,'"srclist"')
point_source_tag = false_test(options.point_source_tag,'"point_source_tag"')

if options.band_nums == 'all':
    band_nums = range(1,25)
else:
    try:
        band_nums = map(int,options.band_nums.split(','))
    except:
        print '-----------------------------------------------------'
        print 'Failed to convert --band_nums into something sensible. Exiting now!!'
        print '-----------------------------------------------------'
        exit()

telescope = options.telescope
time_int = options.time_int
ini_file = options.ini_file
jobs_per_GPU = int(options.jobs_per_GPU)

diffuse_tag = options.diffuse_tag

if options.diffuse_output_dir:
    diffuse_output_dir = options.diffuse_output_dir
else:
    diffuse_output_dir = output_dir

cwd = os.getcwd()

wd = cwd+'/qsub_oskar'

if not os.path.exists(wd):
    os.makedirs(wd)
    
os.chdir(wd)

if not os.path.exists(wd+'/tmp'):
    os.makedirs(wd+'/tmp')
    
point_qsubs = []
diffuse_qsubs = []

##So run 8 OSKAR jobs per nvidia-cuda-mps-control
##4 jobs per GPU
##Each requires a separate PBS job

num_jobs = int(floor(len(band_nums) / (2*jobs_per_GPU)))

print num_jobs

if num_jobs == 0:
    num_jobs = 1

for job_num in xrange(num_jobs):
    job_bands = band_nums[job_num*(2*jobs_per_GPU):(job_num+1)*(2*jobs_per_GPU)]
    print job_bands
    
    point_qsub = write_oskar(wd=wd, metafits=metafits, srclist=srclist, point_source_tag=point_source_tag, time=time, job_bands=job_bands, data_dir=output_dir, telescope=telescope, time_int=time_int, ini_file=ini_file)
    point_qsubs.append(point_qsub)
    
    #if diffuse_tag:
        #base_uvfits = output_dir + '/' + point_source_tag
        #diffuse_qsub1 = write_add_diffuse(wd=wd, metafits=metafits, diffuse_tag=diffuse_tag, time=time, band_num=band_num1, data_dir=diffuse_output_dir, telescope=telescope, time_int=time_int, base_uvfits=base_uvfits)
        #diffuse_qsubs.append(diffuse_qsub1)
        
        #diffuse_qsub2 = write_add_diffuse(wd=wd, metafits=metafits, diffuse_tag=diffuse_tag, time=time, band_num=band_num2, data_dir=diffuse_output_dir, telescope=telescope, time_int=time_int, base_uvfits=base_uvfits)
        #diffuse_qsubs.append(diffuse_qsub2)
        
os.chdir(cwd)

###Write out a controlling bash script to launch all the jobs
out_file = open('run_all_oskarsim.sh','w+')
out_file.write('#!/bin/bash\n')


#if diffuse_tag:
    #for point, diffuse in zip(point_qsubs,diffuse_qsubs):
        #out_file.write('POINT_RUN=$(qsub %s/%s | cut -d "." -f 1)\n' %(wd,point))
        #out_file.write('echo "%s is job "$POINT_RUN\n' %point)
        #out_file.write('DIFFUSE_RUN=$(qsub -Wdepend=afterok:$POINT_RUN %s/%s | cut -d "." -f 1)\n' %(wd,diffuse))
        #out_file.write('echo "%s is job "$DIFFUSE_RUN\n' %diffuse)
        
#else:
for point in point_qsubs:
    out_file.write('POINT_RUN=$(qsub %s/%s | cut -d "." -f 1)\n' %(wd,point))
    out_file.write('echo "%s is job "$POINT_RUN\n' %point)
        
#out_file.close()
