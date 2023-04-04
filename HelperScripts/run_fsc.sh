#!/bin/bash 
#SBATCH --job-name=fastsimcoal2
#SBATCH --partition=main
#SBATCH -N1
#SBATCH -n1                        ## 28 cores on amarel
#SBATCH --cpus-per-task=28         ## 28 cores on amarel 
#SBATCH --exclusive
############################
#SBATCH --mem=120000               ## 128GB /node of amarel
#SBATCH --time=24:00:00            ## max time is 3 days:  3-00:00:00  or  72:00:00
#SBATCH --export=ALL
#SBATCH --no-requeue       ## optional: this prevent the job from restarting if the node fails or the jobs is preempted.
#SBATCH -o slurm-%j.out   ## the default file name is "slurm-%j.out"
#SBATCH -e slurm-%j.err
#SBATCH --mail-type=END                # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jennifer.hoey@rutgers.edu # Email to which notifications will be sent

#############################################################
# standard output is saved in a file:  slurm-$SLURM_JOBID.out 
############################################################# 

mkdir -p /home/$USER/slurm-out/$SLURM_JOBID
cp /home/$USER/fsc25221 /home/$USER/migrationABC_DNA1millionby100.tpl /home/$USER/migrationABC_5729_2_0.5dispersal.est /home/$USER/slurm-out/$SLURM_JOBID
cd       /home/$USER/slurm-out/$SLURM_JOBID

#Create a personal dirctory on the node scratch disk if needed.
mkdir -p /mnt/scratch/$USER/$SLURM_JOBID


# this will tell you when the job started and the host. 
date=`date "+%Y.%m.%d-%H.%M.%S"`
hostname=`hostname`

# to print the variable -- echo
echo $date
echo $hostname

# obtain the environment variables 
# this is useful for reference and troubleshooting issues.

env >               /home/$USER/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-all.out
env | grep SLURM >  /home/$USER/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-SLURM.out


# start time
date

# start the simulation

srun -N1 -n1 --cpus-per-task=28 --exclusive --mem=120000 ./fsc25221 -t migrationABC_DNA1millionby100.tpl -e migrationABC_5729_2_0.5dispersal.est -n 1 -E 500 -g -q -s 0 -m -x --cores 28 --numBatches 28  > /home/$USER/slurm-out/$SLURM_JOBID/$SLURM_JOBID-fastsimcoal2.$SLURM_NNODES.$SLURM_NODELIST.$date.a.txt

#cp -R /mnt/scratch/$USER/$SLURM_JOBID/* /home/$USER/slurm-out/$SLURM_JOBID/

#end time
date
