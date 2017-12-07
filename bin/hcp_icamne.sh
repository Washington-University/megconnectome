#!/bin/bash
#
# This bash script allows you to run the analysis pipeline using the
# compiled megconnectome application.
#
# The distribution of the jobs is not performed in this script, but is
# handled in the MATLAB script, which can use the PBS_ARRAYID environment
# variable to select the subject and/or scan. You can run it like this
#
# qsub -l walltime=4:00:00,mem=8gb,vmem=8gb -t 1-100 <scriptname.sh>
#
# This executes the present script in parallel, each instance having a different
# value of PBS_ARRAYID.

#------------------------------
# CODE THAT FINDS AVAILABLE Xvfb PORT AND STARTS THE VIRTUAL BUFFER
xPort=0;
tmpPort=0;
Xfiles=/tmp/.X*-lock
doneFlag=0;
while [ $doneFlag -ne 1 ]
do
        for f in $Xfiles
        do
                tmpPort=$(echo "$f" | awk -F"X" '{print $2}' | awk -F"-" '{print $1}');
                if [[ "$tmpPort" == "$xPort" ]]
                then
                    xPort=$((xPort+1));
                    doneFlag=0;
                    break;
                else
                    doneFlag=1;
                fi
        done
done
Xvfb :$xPort -screen 0  2760x1960x16 &
export curPID=$!
export DISPLAY=:$xPort
#------------------------------

# process the command line arguments, keep spaces intact
shift 0
args=
while [ $# -gt 0 ]; do
  token=`echo "$1" | sed 's/ /\\\\ /g'`   # Add blackslash before each blank
  args="${args} ${token}"
  shift
done

if [ -d /export/matlab/MCR/R2012b/v80 ]; then
# this applies to the CHPC cluster in St Louis
MCRROOT=/export/matlab/MCR/R2012b/v80
fi

if [ -d /opt/MCR/R2012b/v80 ]; then
# this applies to the Nijmegen cluster
MCRROOT=/opt/MCR/R2012b/v80
fi

# determine where the shell and MATLAB scripts are located
MEGCONNECTOME_ROOT=$HOME/matlab/megconnectome
MEGCONNECTOME_BIN=$MEGCONNECTOME_ROOT/bin
MEGCONNECTOME_SCRIPT=$MEGCONNECTOME_ROOT/pipeline_scripts

# determine the name of the pipeline script that should run
PIPELINE=hcp_icamne

$MEGCONNECTOME_BIN/megconnectome.sh $MCRROOT $MEGCONNECTOME_SCRIPT/$PIPELINE.m $args

#------------------------------
#CLOSE VIRTUAL BUFFER 
kill $curPID

