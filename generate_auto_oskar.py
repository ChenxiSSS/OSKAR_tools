from optparse import OptionParser
from sys import exit
import os
from numpy import ceil

def write_oskar(wd=None, metafits=None, srclist=None, point_source_tag=None, time=None, band_num=None, data_dir=None, telescope=None, cadence=None, ini_file=None):
	'''Writes a bash script for each course band to run OSKAR'''
	
	start, finish = map(float,time.split(','))
	num_time_steps = int((finish - start) / float(cadence))
	
	##Takes around 5 mins to do one time step with 27 channels
	hours = ceil((num_time_steps * 5) / 60.0)
	
	file_name = 'qsub_%s_band%02d_t%d-%d.sh' %(point_source_tag,band_num,start,finish)
	out_file = open(file_name,'w+')
	out_file.write('#!/bin/bash\n')
	out_file.write('#PBS -l nodes=1:gpus=1\n')
	out_file.write('#PBS -l walltime=%02d:00:00\n' %hours)
	out_file.write('#PBS -m e\n')
	out_file.write('#PBS -q sstar\n')
	out_file.write('#PBS -A p048_astro\n')

	##Something one of the gstar people told me a long time ago should be included
	##if you only use one GPU.
	out_file.write('#cat $PBS_GPUFILE\n')
	out_file.write('b="`cat $PBS_GPUFILE | cut -c13`"\n')
	out_file.write('export CUDA_VISIBLE_DEVICES=$b\n')

	out_file.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
	out_file.write('source /home/jline/.bash_profile\n')
	out_file.write('cd %s\n' %wd)
	
	oskar_options = "--metafits=%s --srclist=%s --output_name=%s --time=%s --band_num=%s --debug --data_dir=%s" %(metafits, srclist, point_source_tag, time, band_num, data_dir)
	
	if telescope:
		oskar_options += ' --telescope=%s' %telescope
	if cadence:
		oskar_options += ' --twosec=%s' %cadence
	if ini_file:
		oskar_options += ' --ini_file=%s' %ini_file
	
	out_file.write('time python $OSKAR_TOOLS/MWAobs_oskar_ascii.py ' + oskar_options + '\n')
	
	miss_chan_options = "--metafits=%s --output_name=%s --time=%s --band_num=%d --debug --data_dir=%s" %(metafits, point_source_tag, time, band_num, data_dir)
	if telescope:
		miss_chan_options += ' --telescope=%s' %telescope
	if cadence:
		miss_chan_options += ' --twosec=%s' %cadence
	
	out_file.write('time python $OSKAR_TOOLS/cp_missing_chans.py ' + miss_chan_options + '\n')
		
	out_file.close()
	return file_name

def write_add_diffuse(wd=None, metafits=None, diffuse_tag=None, time=None, band_num=None, data_dir=None, telescope=None, cadence=None, base_uvfits=None):
	''' '''
	
	start, finish = map(float,time.split(','))
	num_time_steps = int((finish - start) / float(cadence))
	
	##Takes around 10 mins to do one time step with 27 channels
	hours = ceil((num_time_steps * 10) / 60.0)
	
	file_name = 'qsub_%s_band%02d_t%d-%d.sh' %(diffuse_tag,band_num,start,finish)
	out_file = open(file_name,'w+')
	out_file.write('#!/bin/bash\n')
	out_file.write('#PBS -l nodes=1\n')
	out_file.write('#PBS -l walltime=%02d:00:00\n' %hours)
	out_file.write('#PBS -m e\n')
	out_file.write('#PBS -q sstar\n')
	out_file.write('#PBS -A p048_astro\n')

	out_file.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
	out_file.write('source /home/jline/.bash_profile\n')
	out_file.write('cd %s\n' %wd)
	
	majick_options = "--metafits=%s --output_name=%s --time=%s --band_num=%s --debug --data_loc=%s --base_uvfits=%s" %(metafits, diffuse_tag, time, band_num, data_dir, base_uvfits)
	
	if telescope:
		majick_options += ' --telescope=%s' %telescope
	if cadence:
		majick_options += ' --twosec=%s' %cadence
	
	out_file.write('time python $OSKAR_TOOLS/MWAobs_add_MAJICK_diffuse.py ' + majick_options + '\n')
	
	miss_chan_options = "--metafits=%s --output_name=%s --time=%s --band_num=%d --debug --data_dir=%s" %(metafits, diffuse_tag, time, band_num, data_dir)
	if telescope:
		miss_chan_options += ' --telescope=%s' %telescope
	if cadence:
		miss_chan_options += ' --twosec=%s' %cadence
	
	out_file.write('time python $OSKAR_TOOLS/cp_missing_chans.py ' + miss_chan_options + '\n')
		
	out_file.close()
	return file_name
	

parser = OptionParser()

parser.add_option('-o', '--output_dir', default=False, help='Enter output data directory')
parser.add_option('-m', '--metafits', default=False, help='Enter metafits file to base observation on')
parser.add_option('-t', '--time', default=False, help='Enter start,finish times relative to metafits date (seconds - i.e. 0 to start at beginnning')
parser.add_option('-c', '--cadence',default=False, help='Enter cadence of correlator to simulate (s) to override what is in metafits (i.e --cadence=2). Defaults to what is in the metafits')
parser.add_option('-s','--srclist', default=False, help='Enter RTS srclist to use as sky model')
parser.add_option('-p','--point_source_tag', default=False, help='Enter uvfits tag for the point source only model')
parser.add_option('-d','--diffuse_tag', default=False, help='Enter uvfits tag for the point source + diffuse models')
parser.add_option('-l','--diffuse_output_dir', default=False, help='Enter output data directory for diffuse - defaults to same location as point source model')
parser.add_option('-b', '--num_bands', default=False, help='Enter number of course channels to simulate')
parser.add_option('-a','--telescope', default=False, help='Enter telescope tag to use - default is MWA_phase1')
parser.add_option('-i', '--ini_file', default=False, help='Enter template oskar .ini - defaults to the template .ini located in $OSKAR_TOOLS/telescopes/--telescope')

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
num_bands = int(false_test(options.num_bands,'"num_bands"'))

telescope = options.telescope
cadence = options.cadence
ini_file = options.ini_file

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

for band_num in xrange(num_bands):
	band_num += 1
	point_qsub = write_oskar(wd=wd, metafits=metafits, srclist=srclist, point_source_tag=point_source_tag, time=time, band_num=band_num, data_dir=output_dir, telescope=telescope, cadence=cadence, ini_file=ini_file)
	point_qsubs.append(point_qsub)
	
	if diffuse_tag:
		base_uvfits = output_dir + '/' + point_source_tag
		diffuse_qsub =  write_add_diffuse(wd=wd, metafits=metafits, diffuse_tag=diffuse_tag, time=time, band_num=band_num, data_dir=diffuse_output_dir, telescope=telescope, cadence=cadence, base_uvfits=base_uvfits)
		diffuse_qsubs.append(diffuse_qsub)
		
os.chdir(cwd)

###Write out a controlling bash script to launch all the jobs
out_file = open('run_all_oskarsim.sh','w+')
out_file.write('#!/bin/bash\n')


if diffuse_tag:
	for point, diffuse in zip(point_qsubs,diffuse_qsubs):
		out_file.write('POINT_RUN=$(qsub %s | cut -d "." -f 1)\n' %point)
		out_file.write('echo "%s is job "$MAIN_RUN\n' %point)
		out_file.write('DIFFUSE_RUN=$(qsub -W --depend:afterok:$POINT_RUN %s | cut -d "." -f 1)\n' %diffuse)
		out_file.write('echo "%s is job "$DIFFUSE_RUN\n' %diffuse)
		
else:
	for point in point_qsubs:
		out_file.write('POINT_RUN=$(qsub %s | cut -d "." -f 1)\n' %point)
		out_file.write('echo "%s is job "$MAIN_RUN\n' %point)
		
out_file.close()