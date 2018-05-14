#!/usr/bin/env python
from optparse import OptionParser
from sys import exit
import os
from numpy import ceil,floor

OSKAR_dir = os.environ['OSKAR_TOOLS']

def write_oskar(wd=None, metafits=None, srclist=None, oskar_uvfits_tag=None, time=None, band_num=None, 
                data_dir=None, telescope=None, time_int=None, ini_file=None, jobs_per_GPU=None,
                flag_dipoles=None, cluster=None,retain_vis_file=None,retain_ini_file=None,
                do_phase_track=False):):
    '''Writes a bash script for each course band to run OSKAR'''
    
    start, finish = map(float,time.split(','))
    num_time_steps = int((finish - start) / float(time_int))
    
    if srclist:
        ##Takes around 7.5 mins to do one time step with 27 channels
        hours = ceil((num_time_steps * 7.5) / 60.0)
    elif osm:
        ##Quicker, run all channels with same sky model
        hours = ceil((num_time_steps * 3.0) / 60.0)
    
    ##Set up controlling qsub script to run current batch of OSKAR processes
    file_name = 'sbatch_%s_band%02d_t%d-%d.sh' %(oskar_uvfits_tag,band_num,start,finish)


    out_file = open(file_name,'w+')
    if cluster == 'spartan':

        out_file.write('#!/bin/sh\n')
        out_file.write('#SBATCH -p gpgpu\n')
        out_file.write('#SBATCH --time=%02d:00:00\n' %hours)
        out_file.write('#SBATCH -A punim0411\n')
        out_file.write('#SBATCH --gres=gpu:1\n')
        #out_file.write('#SBATCH --mincpus=16\n')
        out_file.write('#SBATCH --mem=16384\n')
        out_file.write('source /data/cephfs/punim0411/software/OSKAR_tools/cluster_modles/load_spartan.sh\n')
        out_file.write('source /data/cephfs/punim0411/software/OSKAR_tools/init_OSKAR_tools.sh\n')

    elif cluster == 'ozstar':

        out_file.write('#!/bin/bash -l\n')
        out_file.write('#SBATCH --partition=skylake-gpu\n')
        out_file.write('#SBATCH --time=%02d:00:00\n' %hours)
        out_file.write('#SBATCH --account=oz048\n')
        out_file.write('#SBATCH --gres=gpu:1\n')
        out_file.write('#SBATCH --mem=16384\n')
        out_file.write('source /home/jline/software/OSKAR_tools/cluster_modles/load_ozstar.sh\n')
        out_file.write('source /home/jline/software/OSKAR_tools/init_OSKAR_tools.sh\n')

    out_file.write('cd %s\n' %wd)

    
    #half_of_jobs = len(job_bands) / 2

    ###Setup half of jobs to run on one GPU, half on the other
    #run1_name = 'run_%s_bands%02d-%02d_t%d-%d.sh' %(oskar_uvfits_tag,job_bands[0],job_bands[half_of_jobs-1],start,finish)
    #run1 = open(run1_name,'w+')
    ##run1.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
    
    #run2_name = 'run_%s_bands%02d-%02d_t%d-%d.sh' %(oskar_uvfits_tag,job_bands[half_of_jobs],job_bands[-1],start,finish)
    #run2 = open(run2_name,'w+')
    
    
    #def write_oskar_command(band_num=None,runfile=None):
    oskar_options = "--metafits=%s --output_name=%s --time=%s --band_num=%s --debug --data_dir=%s" %(metafits, oskar_uvfits_tag, time, band_num, data_dir)
    if srclist:
        oskar_options += ' --srclist=%s' %srclist
    elif osm:
        oskar_options += ' --osm=%s' %osm
    else:
        print "Neither a srclist nor an osm model provided - you have nothing to simulate! Exiting"
        exit()
    
    if telescope:
        oskar_options += ' --telescope=%s' %telescope
    if time_int:
        oskar_options += ' --time_int=%s' %time_int
    if ini_file:
        oskar_options += ' --ini_file=%s' %ini_file
    if flag_dipoles:
        oskar_options += ' --flag_dipoles'
    if retain_vis_file:
        oskar_options += ' --retain_vis_file'
    if retain_ini_file:
        oskar_options += ' --retain_ini_file'
    if do_phase_track:
        oskar_options += ' --do_phase_track'
    
    out_file.write('time %s/MWAobs_oskar.py %s &\n' %(OSKAR_dir, oskar_options))
        
    #for band_num in job_bands[:half_of_jobs]:
        #write_oskar_command(band_num=band_num,runfile=run1)
        
    #for band_num in job_bands[half_of_jobs:]:
        #write_oskar_command(band_num=band_num,runfile=run2)
        
    #run1.close()
    #run2.close()
    
    ###Do the fancy CUDA mpi control thing and actually run all the jobs
    #out_file.write('time source %s/run_mps_oskar_slurm.sh %s/%s %s/%s\n' %(OSKAR_dir,wd,run1_name,wd,run2_name))
    out_file.write('wait\n')
    out_file.write('rm %s/%s' %(wd,file_name))
    out_file.close()
    
    return file_name

#def write_add_diffuse(wd=None, metafits=None, majick_tag=None, time=None, band_num=None, data_dir=None, telescope=None, time_int=None, base_uvfits=None):
    #''' '''
    
    #start, finish = map(float,time.split(','))
    #num_time_steps = int((finish - start) / float(time_int))
    
    ###Takes around 10 mins to do one time step with 27 channels
    #hours = ceil((num_time_steps * 15) / 60.0)
    
    #file_name = '_%s_band%02d_t%d-%d.sh' %(majick_tag,band_num,start,finish)
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
    
    #majick_options = "--metafits=%s --output_name=%s --time=%s --band_num=%s --debug --data_loc=%s --base_uvfits=%s" %(metafits, majick_tag, time, band_num, data_dir, base_uvfits)
    
    #if telescope:
        #majick_options += ' --telescope=%s' %telescope
    #if time_int:
        #majick_options += ' --twosec=%s' %time_int
    
    #out_file.write('time python $OSKAR_TOOLS/MWAobs_add_MAJICK_diffuse.py ' + majick_options + '\n')
    
    #out_file.close()
    #return file_name
    

parser = OptionParser()

parser.add_option('-a','--telescope', default=False, help='Enter telescope tag to use - default is MWA_phase1')
parser.add_option('-b', '--band_nums', default='all', help='Defaults to running all 24 course bands. Alternatively, enter required numbers delineated by commas, e.g. --band_nums=1,7,9')
parser.add_option('-c', '--time_int',default=False, help='Enter time_int of correlator to simulate (s) to override what is in metafits (i.e --time_int=2). Defaults to what is in the metafits')
parser.add_option('-d','--majick_tag', default=False, help='Enter uvfits tag for the point source + diffuse models')
parser.add_option('-e','--osm', default=False, help='Alternatively just use an OSKAR .osm model')
parser.add_option('-f','--flag_dipoles',default=False,action='store_true', help='Add to switch on dipole flagging via the metafits file. NOTE needs a metafits that has the correct data within')
parser.add_option('-g', '--cluster', default='spartan', help='Enter the super cluster name you are on - default is spartan, options are spartan or ozstar')


parser.add_option('-i', '--ini_file', default=False, help='Enter template oskar .ini - defaults to the template .ini located in $OSKAR_TOOLS/telescopes/--telescope')
parser.add_option('-l','--diffuse_output_dir', default=False, help='Enter output data directory for diffuse - defaults to same location as point source model')
parser.add_option('-m', '--metafits', default=False, help='Enter metafits file to base observation on')
parser.add_option('-o', '--output_dir', default=False, help='Enter output data directory')
parser.add_option('-p','--oskar_uvfits_tag', default=False, help='Enter uvfits tag for the point source only model')
parser.add_option('-s','--srclist', default=False, help='Enter RTS srclist to use as sky model')
parser.add_option('-t', '--time', default=False, help='Enter start,finish times relative to metafits date (seconds - i.e. 0 to start at beginnning')

parser.add_option('--retain_vis_file',default=False,action='store_true', help='Add to not delete the oskar binary .vis files')
parser.add_option('--retain_ini_file',default=False,action='store_true', help='Add to not delete the oskar binary .ini files')
parser.add_option('--do_phase_track',default=False,action='store_true', help='Add to leave on the phase tracking done by OSKAR')




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
#srclist = false_test(options.srclist,'"srclist"')
time_int = false_test(options.time_int,'"time_int"')
oskar_uvfits_tag = false_test(options.oskar_uvfits_tag,'"oskar_uvfits_tag"')

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
ini_file = options.ini_file
srclist = options.srclist
osm = options.osm
flag_dipoles = options.flag_dipoles
retain_vis_file = options.retain_vis_file
retain_ini_file = options.retain_ini_file
do_phase_track = options.do_phase_track

if srclist:
    pass
elif osm:
    pass
else:
    print "Neither a srclist nor an osm model provided - you have nothing to simulate! Exiting"
    exit()

majick_tag = options.majick_tag

if options.diffuse_output_dir:
    diffuse_output_dir = options.diffuse_output_dir
else:
    diffuse_output_dir = output_dir

cwd = os.getcwd()

wd = cwd+'/slurm_oskar'

if not os.path.exists(wd):
    os.makedirs(wd)
    
os.chdir(wd)

if not os.path.exists(wd+'/tmp'):
    os.makedirs(wd+'/tmp')
    
oskar_slurms = []
majick_slurms = []

    
    
for band_num in band_nums:
    oskar_slurm = write_oskar(wd=wd, metafits=metafits, srclist=srclist, oskar_uvfits_tag=oskar_uvfits_tag, time=time, band_num=band_num, data_dir=output_dir,
        telescope=telescope, time_int=time_int, ini_file=ini_file, flag_dipoles=flag_dipoles, cluster=options.cluster, retain_vis_file=retain_vis_file,
        retain_ini_file=retain_ini_file, do_phase_track=do_phase_track)
    
    oskar_slurms.append(oskar_slurm)
    
    #if majick_tag:
        #base_uvfits = output_dir + '/' + oskar_uvfits_tag
        #majick_slurm1 = write_add_diffuse(wd=wd, metafits=metafits, majick_tag=majick_tag, time=time, band_num=band_num1, data_dir=diffuse_output_dir, telescope=telescope, time_int=time_int, base_uvfits=base_uvfits)
        #majick_slurms.append(majick_slurm1)
        
        #majick_slurm2 = write_add_diffuse(wd=wd, metafits=metafits, majick_tag=majick_tag, time=time, band_num=band_num2, data_dir=diffuse_output_dir, telescope=telescope, time_int=time_int, base_uvfits=base_uvfits)
        #majick_slurms.append(majick_slurm2)
        
os.chdir(cwd)

###Write out a controlling bash script to launch all the jobs
out_file = open('run_all_oskarsim.sh','w+')
out_file.write('#!/bin/bash\n')


#if majick_tag:
    #for point, diffuse in zip(oskar_slurms,majick_slurms):
        #out_file.write('OSKAR_RUN=$(qsub %s/%s | cut -d "." -f 1)\n' %(wd,point))
        #out_file.write('echo "%s is job "$OSKAR_RUN\n' %point)
        #out_file.write('MAJICK_RUN=$(qsub -Wdepend=afterok:$OSKAR_RUN %s/%s | cut -d "." -f 1)\n' %(wd,diffuse))
        #out_file.write('echo "%s is job "$MAJICK_RUN\n' %diffuse)
        
#else:
for point in oskar_slurms:
    out_file.write('OSKAR_RUN=$(sbatch %s/%s | cut -d "." -f 1)\n' %(wd,point))
    out_file.write('echo "%s is job "$OSKAR_RUN\n' %point)
        
#out_file.close()
