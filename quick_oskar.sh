#!/bin/bash

DEVICE_1="$(echo $CUDA_VISIBLE_DEVICES | cut -f1 -d,)"
DEVICE_2="$(echo $CUDA_VISIBLE_DEVICES | cut -f2 -d,)"

SCRIPT_1=$1
SCRIPT_2=$2

######## BP: First part sets up the scheduler. Get GPU IDs from environment variables
for i in $DEVICE_1 $DEVICE_2; do
    mkdir /lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-$i /lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-log-$i
    export CUDA_VISIBLE_DEVICES=$i # SELECT GPU
    export CUDA_MPS_PIPE_DIRECTORY=/lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-$i # NAMED PIPES
    export CUDA_MPS_LOG_DIRECTORY=/lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-log-$i # LOGFILES
    nvidia-cuda-mps-control -d # START MPS DAEMON
done

## IF A CUDA PROGRAM CONNECTS TO THE MPS DAEMON, THE DAEMON CREATES AN
## MPS PROXY SERVER FOR THE CONNECTING CLIENT IF NOT ALREADY PRESENT.
## THE PROXY SERVER USES THE CORRECT GPU DEVICE.

unset CUDA_VISIBLE_DEVICES

#####  BP: Here the actual job is run. The scheduler is chosen by passing the CUDA_MPS_PIPE_DIRECTORY through the mpi -x flag. Probably, you can just set it on the command line before launching oskar

export CUDA_MPS_PIPE_DIRECTORY=/lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-$DEVICE_1
source $SCRIPT_1

export CUDA_MPS_PIPE_DIRECTORY=/lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-$DEVICE_2
source $SCRIPT_2

##Wait means bash will not launch next job until all the above has finished
wait

##Now jobs have finished, kill the daemon
##########BP: Shutdown
for i in $DEVICE_1 $DEVICE_2; do
export CUDA_MPS_PIPE_DIRECTORY=/lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-$i # SELECT MPS DAEMON
echo quit | nvidia-cuda-mps-control # STOP MPS DAEMON
rm -rf /lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-$i /lustre/projects/p048_astro/jline/OSKAR_sims/quick_OSKAR/tmp/mps-log-$i
done
