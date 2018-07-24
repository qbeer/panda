#!/bin/bash
#SBATCH -J pndana
#SBATCH --time=8:00:00
#SBATCH --get-user-env
#SBATCH -e data/slurmlog/slurm_%j_errout.log
#SBATCH -o data/slurmlog/slurm_%j_errout.log

if [ $# -lt 1 ]; then
  echo -e "\nJob script for submission of PandaRoot FastSim with analysis jobs based on macro 'quickfsimana.C' on KRONOS. *The macro needs to configured beforehand!*\n"
#  echo -e "   ********************************************************"  
#  echo -e "   *** RECOMMENDED: Use anasub.pl for easier submission ***"  
#  echo -e "   *********************************************************\n"  
  echo -e "USAGE: sbatch -a<min>-<max> jobquickfa_kronos.sh <prefix> '<macro>'\n"
  echo -e " <min>     : Minimum job number of files data/<prefix>_<min>_pid.root"
  echo -e " <max>     : Maximum job number of files data/<prefix>_<max>_pid.root"
  echo -e " <prefix>  : Prefix of output files"
  echo -e " <macro>   : Complete call of quickfsimana.C() macro with parameters. Prefix should be PREFIX, run number should be RUN.\n"
  echo -e "Example 1 : sbatch -a1-10 jobquickfa_kronos.sh signal 'quickfsimana.C(\"PREFIX\",\"jpsi2pi.dec\",6.23,\"J/psi->mu+ mu-;pbp->J/psi pi+ pi-\",1000,\"fit4c:mwin=0.6:pidmu=Loose\",0,RUN,1)'"
  echo -e "Example 2 : sbatch -a1-10 jobquickfa_kronos.sh dpmbkg 'quickfsimana.C(\"PREFIX\",\"DPM\",6.23,\"J/psi->mu+ mu-;pbp->J/psi pi+ pi-\",1000,\"fit4c:mwin=0.6:pidmu=Loose\",0,RUN,0)'\n"
  
  exit 1
fi

nyx=$VMCWORKDIR"/tutorials/rho"
_target=$nyx"/data/"

prefix=""
macro=""
run=$SLURM_ARRAY_TASK_ID

# create tmp directory to write stuff in
tmpdir="/tmp/"$USER"_"$SLURM_JOB_ID"_"$run"/"
mkdir $tmpdir
cd $tmpdir

# check input parameters
if test "$1" != ""; then
  prefix=$1
fi

if test "$2" != ""; then
  macro=$2
fi

# variable to hold prefix with run number
prefrun=$prefix"_"$run

# replace place-holders PREFIX and RUN with actual values, in particular add absolute paths
macro=${macro/PREFIX/$prefix}
macro=${macro/RUN/$run}

# if there is a dec file given with relative path, add the absolute path
decname=`expr "$macro" : '.*\"\(.*.dec\)\".*'`

if test "$decname" != ""; then
  if [[ $decname != \/* ]] ; then
    macro=${macro/$decname/$nyx"/"$decname}
	decname=$nyx"/"$decname
  fi
fi

echo "decay file : \"$decname\""

# this is the command to be run
echo "$nyx/$macro"

# output and log files
outname=$prefrun"_ana.root"
logname=$prefrun"_quickfa.log"

# the actual command
root -l -b -q "$nyx/$macro" &> $logname

# check whether everything is there (in slurmlog visible)
ls -lh $tmpdir

# move output and remove tmp dir
rm gphysi.dat
mv $tmpdir/* $_target
rm -rf $tmpdir


