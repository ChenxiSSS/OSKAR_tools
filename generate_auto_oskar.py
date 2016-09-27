from optparse import OptionParser
from sys import exit

parser = OptionParser()

parser.add_option('-g', '--obsID', default=False, help='Enter GPS of obsID')
parser.add_option('-s', '--start', default=False, help='Enter start time relative to metafits date (seconds - i.e. 0 to start at beginnning')
parser.add_option('-f', '--finish',	default=False, help='Enter end time relative to metafits data (seconds - i.e. 112 will finish 112s after metafits date/time)')
parser.add_option('-c', '--cadence',default=False, help='Enter cadence of correlator to simulate (s)')
parser.add_option('-o', '--oskar', default=False, help='Enter template oskar .ini file')
parser.add_option('-m','--srclist', default=False, help='Enter RTS srclist to use as sky model')
parser.add_option('-t','--tag', default=False, help='Enter name tag for files/scripts')
parser.add_option('-d','--MWA_DIR', default='/lustre/projects/p048_astro/jline/OSKAR_sims/PUMA_sims/', help='Destination data directory - must contain obsID dir, containing metafits file')

options, args = parser.parse_args()

##Possibly have this as a environment variable if we ever need to 
##transfer across to another cluster/other people need to use this.
MWA_DIR = options.MWA_DIR
if MWA_DIR[-1] == '/': MWA_DIR =  MWA_DIR[:-1]

def false_test(option,name):
	if option == False:
		print '-----------------------------------------------------'
		print '%s not entered but is required. Exiting now!!' %name
		print '-----------------------------------------------------'
		exit()
	else:
		return option

##Get inputs
gps = false_test(options.obsID,'"obsID"')
start = float(false_test(options.start,'"start"'))
finish = float(false_test(options.finish,'"finish"'))
cadence = float(false_test(options.cadence,'"cadence"'))
oskar = false_test(options.oskar,'"oskar"')
srclist = false_test(options.srclist,'"srclist"')
tag = false_test(options.tag,'"tag"')

def write_main(gps,start,finish,cadence,tag):
	'''Writes the bash script to run the majority of the simulations and
	the conversion of OSKAR outputs in to uvfits files'''
	
	##Safe to have 1 hour to generate 10 time steps for 24 bands, 27 channels 
	##(don't sim flagged channels)
	num_time_steps = int((finish - start) / cadence)
	hours = num_time_steps / 10
	if num_time_steps % 10 > 0: hours +=1
	
	file_name = '%s_%s_oskarsim_%d-%d.sh' %(tag,gps,start,finish)
	out_file = open(file_name,'w+')
	out_file.write('#!/bin/bash\n')
	out_file.write('#PBS -l nodes=1:gpus=1\n')
	out_file.write('#PBS -l walltime=%02d:00:00\n' %hours)
	out_file.write('#PBS -t 0-23\n')
	out_file.write('#PBS -m e\n')
	out_file.write('#PBS -q gstar\n')
	out_file.write('#PBS -A p048_astro\n')

	out_file.write('#cat $PBS_GPUFILE\n')
	out_file.write('b="`cat $PBS_GPUFILE | cut -c13`"\n')
	out_file.write('export CUDA_VISIBLE_DEVICES=$b\n')

	out_file.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
	out_file.write('source /home/jline/.bash_profile\n')
	out_file.write('cd %s\n' %MWA_DIR)
	out_file.write('time python /lustre/projects/p048_astro/jline/OSKAR_sims/MWAobs_oskar.py -m %s/%s/%s.metafits -s %s -i %s -o %s/data/%s -t %.1f,%.1f -x %.1f -d -n %s' %(MWA_DIR,gps,gps,srclist,oskar,MWA_DIR,gps,start,finish,cadence,tag))
	out_file.close()
	return file_name
	
def fill_missing(gps,start,finish,cadence,tag):
	'''Writes the bash script to find any failed channels, run them and
	the conversion of OSKAR outputs in to uvfits files'''
	
	file_name = '%s_%s_fillmissing_%d-%d.sh' %(tag,gps,start,finish)
	out_file = open(file_name,'w+')
	out_file.write('#!/bin/bash\n')
	out_file.write('#PBS -l nodes=1:gpus=1\n')
	out_file.write('#PBS -l walltime=04:00:00\n')
	out_file.write('#PBS -t 0-23\n')
	out_file.write('#PBS -m e\n')
	out_file.write('#PBS -q gstar\n')
	out_file.write('#PBS -A p048_astro\n')

        out_file.write('#cat $PBS_GPUFILE\n')
        out_file.write('b="`cat $PBS_GPUFILE | cut -c13`"\n')
        out_file.write('export CUDA_VISIBLE_DEVICES=$b\n')

	out_file.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
	out_file.write('source /home/jline/.bash_profile\n')
	out_file.write('cd %s\n' %MWA_DIR)
	out_file.write('time python /lustre/projects/p048_astro/jline/OSKAR_sims/fill_missing_chans.py -m %s/%s/%s.metafits -s %s -i %s -o %s/data/%s -t %.1f,%.1f -x %.1f -d -n %s' %(MWA_DIR,gps,gps,srclist,oskar,MWA_DIR,gps,start,finish,cadence,tag))
	out_file.close()
	return file_name

def cp_flagged(gps,start,finish,cadence,tag):
	'''Writes the bash script to make copies of uvfits
	for the flagged channels inputs - RTS will crash if they dont exist'''
	
	file_name = '%s_%s_cpflagged.sh' %(tag,gps)
	out_file = open(file_name,'w+')
	out_file.write('#!/bin/bash\n')
	out_file.write('#PBS -l nodes=1\n')
	out_file.write('#PBS -l walltime=00:30:00\n')
	out_file.write('#PBS -t 0-23\n')
	out_file.write('#PBS -m e\n')
	out_file.write('#PBS -q gstar\n')
	out_file.write('#PBS -A p048_astro\n')
	out_file.write('source /lustre/projects/p048_astro/MWA/bin/activate\n')
	out_file.write('source /home/jline/.bash_profile\n')
	out_file.write('cd %s\n' %MWA_DIR)
	out_file.write('time python /lustre/projects/p048_astro/jline/OSKAR_sims/cp_missing_chans.py -m %s/%s/%s.metafits -s %s -i %s -o %s/data/%s -t %.1f,%.1f -x %.1f -d -n %s' %(MWA_DIR,gps,gps,srclist,oskar,MWA_DIR,gps,start,finish,cadence,tag))
	out_file.close()
	return file_name

##MAIN LOOOP=========================================================================
##Gstar fails if you try to make more than 80 time steps simultaneously (in the way
##that the scritps called above run), so split the processing up accordingly,
##depending on how much data the user is trying to create

num_time_steps = int((finish - start) / cadence)

##In case the time of simulation doesn't match the cadence,
##ensure all requested time is simulated
if (finish - start) % cadence > 0: num_time_steps +=1
num_extra = num_time_steps / 80

##This list will split the required time steps up into
##groups of 80 as applicable
time_boundaries = [start]

##If extra scripts are neccessary (num_extra > 0), insert into time_boundaries
for i in xrange(num_extra): time_boundaries.append(start+((i+1)*(80*cadence)))

##Finally, add the finish time
time_boundaries.append(finish)

main_files = []
missing_files = []

for i in xrange(len(time_boundaries) - 1):
	##Write the bash scripts to run each 80 timestep chunk
	scr_start = time_boundaries[i]
	scr_finish = time_boundaries[i+1]
	main_files.append(write_main(gps,scr_start,scr_finish,cadence,tag))
	missing_files.append(fill_missing(gps,scr_start,scr_finish,cadence,tag))
	
cpflagged = cp_flagged(gps,start,finish,cadence,tag)


##Write out a controlling bash script to launch all the jobs
out_file = open('run_all_oskarsim_%s.sh' %gps,'w+')
out_file.write('#!/bin/bash\n')

for main_file,missing_file in zip(main_files[:1],missing_files[:1]):
	out_file.write('MAIN_RUN=$(qsub %s | cut -d "." -f 1)\n' %main_file)
	out_file.write('echo $MAIN_RUN\n')
	out_file.write('MISSING_RUN=$(qsub -W --depend:afterok:$MAIN_RUN %s | cut -d "." -f 1)\n' %missing_file)
	
if len(main_files) > 1:
	for main_file,missing_file in zip(main_files[1:],missing_files[1:]):
		out_file.write('MAIN_RUN=$(qsub -W --depend:afterok:$MISSING_RUN %s | cut -d "." -f 1)\n' %main_file)
		out_file.write('echo $MAIN_RUN\n')
		out_file.write('MISSING_RUN=$(qsub -W --depend:afterok:$MAIN_RUN %s | cut -d "." -f 1)\n' %missing_file)

out_file.write('qsub -W --depend:afterok:$MISSING_RUN %s\n' %cpflagged)
out_file.close()
