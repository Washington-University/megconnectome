#!/bin/sh
#
# This script rearranges the MEG data and pipeline results into 
# an organization that is suited for release.
#
# Use as
#   hcp_megrelease.sh [options] <inputdir> <outputdir> 
# 
# The inputdir should point to the location of a single subjects
# data, INCLUDING the experiment name.
# 
# The outputdir should point to a directory where the release
# should be compiled, also INCLUDING the experiment name.
# 
# The selection of specific types of data is done with
#   --raw           include the unprocessed raw data in the release
#   --anatomy       include the results from the anatomy pipeline
#   --eprime        include the EPrime log files in the release
#   --datacheck     include the specific pipeline results
#   --baddata       include the specific pipeline results
#   --icaclass      include the specific pipeline results
#   --rmegpreproc   include the specific pipeline results
#   --tmegpreproc   include the specific pipeline results
#   --eravg         include the specific pipeline results
#   --tfavg         include the specific pipeline results
#   --powavg        include the specific pipeline results
#   --icamne        include the specific pipeline results
#   --icablpenv     include the specific pipeline results
#   --icablpcorr    include the specific pipeline results
#   --icaimagcoh    include the specific pipeline results
#   --srcavglcmv    include the specific pipeline results
#   --srcavgdics    include the specific pipeline results
# The default behaviour when none of these options is specified
# is to include all types of data.
# 
# The selection of specific types of scans is done with
#  --noise          include the noise measurements
#  --restin         include the resting state data (note the "g" missing)
#  --wrkmem         include the working memory task data
#  --storym         include the story-math task data
#  --motort         include the motor task data (note the "t" at the end)
# The default behaviour when none of these options is specified
# is to include all types of scans.
#
# For example
#   hcp_megrelease.sh /HCP/BlueArc/meg_scratch/intradb/archive1/HCP_Phase2/arc001/177746_MEG /HCP/BlueArc/meg_scratch/release/177746_MEG 
# to release everything, or
#   hcp_megrelease.sh --raw --anatomy /HCP/BlueArc/meg_scratch/intradb/archive1/HCP_Phase2/arc001/177746_MEG /HCP/BlueArc/meg_scratch/release/177746_MEG
# to release the raw data and the anatomical models.
#
# A more elaborate example to release the raw and anatomical data for multiple subjects is
#
#  for SUBJECT in 104820 122317 177746 500222 
#  do
#    INPUTDIR=/HCP/BlueArc/meg_scratch/intradb/archive1/HCP_Phase2/arc001/"$SUBJECT"_MEG 
#    OUTPUTDIR=/HCP/BlueArc/meg_scratch/release/"$SUBJECT"
#    hcp_megrelease.sh --anatomy --raw $INPUTDIR $OUTPUTDIR
#  done

# the option parsing is based on http://mywiki.wooledge.org/BashFAQ/035

function warning {
  echo WARNING: $* >&2
}

function error {
  echo ERROR: $* >&2
  exit 1
}

function copy {
  if [ ! -d $2 ]; then
    error cannot write to output directory $2
    return
  fi

  if [ ! -f "$1" ] ; then 
   warning $1 is not a regular file
   return
  fi

  if [ ! -r "$1" ] ; then 
    warning cannot read from file $1
   return
  fi

  echo copying $1 to $2
  rsync -a $1 $2
}

# Reset all variables that might be set
raw="true"
anatomy="true"
eprime="true"
datacheck="true"
baddata="true"
icaclass="true"
rmegpreproc="true"
tmegpreproc="true"
eravg="true"
tfavg="true"
powavg="true"
icamne="true"
icablpenv="true"
icablpcorr="true"
icaimagcoh="true"
srcavglcmv="true"
srcavgdics="true"

noise="true"
restin="true"
wrkmem="true"
storym="true"
motort="true"

# first loop, decide between default/all or explicit selection
for ARG in $@
do
    case $ARG in 
      --raw | --anatomy | --eprime | --datacheck | --baddata | --icaclass | --rmegpreproc | --tmegpreproc | --eravg | --tfavg | --powavg | --icamne | --icablpenv | --icablpcorr | --icaimagcoh | --srcavglcmv | --srcavgdics)
          # an explicit subselection has been made for the data type
          raw="false"
          anatomy="false"
          eprime="false"
          datacheck="false"
          baddata="false"
          icaclass="false"
          rmegpreproc="false"
          tmegpreproc="false"
          eravg="false"
          tfavg="false"
          powavg="false"
          icamne="false"
          icablpenv="false"
          icablpcorr="false"
          icaimagcoh="false"
          srcavglcmv="false"
          srcavgdics="false"
          ;;

      --noise | --restin | --wrkmem | --storym | --motort)
          # an explicit subselection has been made for the scan type
          noise="false"
          restin="false"
          wrkmem="false"
          storym="false"
          motort="false"
          ;;

      *)
         break
          ;;
  esac
done

# second loop, parse all arguments
while :
do
    case $1 in
      -h | --help | -\?)
          echo
          cat $0 | head -n 51 | tail -n 49 | tr -d \#
          echo
          exit 0      # This is not an error, User asked help. Don't do "exit 1"
          ;;

      --raw)
          raw="true"
          shift
          ;;

      --anatomy)
          anatomy="true"
          shift
          ;;

      --eprime)
          eprime="true"
          shift
          ;;

      --datacheck)
          datacheck="true"
          shift
          ;;

      --baddata)
          baddata="true"
          shift
          ;;

      --icaclass)
          icaclass="true"
          shift
          ;;

      --rmegpreproc)
          rmegpreproc="true"
          shift
          ;;

      --tmegpreproc)
          tmegpreproc="true"
          shift
          ;;

      --eravg)
          eravg="true"
          shift
          ;;

      --tfavg)
          tfavg="true"
          shift
          ;;

      --powavg)
          powavg="true"
          shift
          ;;

      --icamne)
          icamne="true"
          shift
          ;;

      --icablpenv)
          icablpenv="true"
          shift
          ;;

      --icablpcorr)
          icablpcorr="true"
          shift
          ;;

      --icaimagcoh)
          icaimagcoh="true"
          shift
          ;;

      --srcavglcmv)
          srcavglcmv="true"
          shift
          ;;

      --srcavgdics)
          srcavgdics="true"
          shift
          ;;

       --noise)
          noise="true"
          shift
          ;;

       --restin)
          restin="true"
          shift
          ;;

       --wrkmem)
          wrkmem="true"
          shift
          ;;

       --storym)
          storym="true"
          shift
          ;;

       --motort)
          motort="true"
          shift
          ;;

        --) # End of all options
            shift
            break
            ;;

        -*)
            warning "Unknown option (ignored): $1" >&2
            shift
            ;;

        *)  # no more options. Stop while loop
            break
            ;;
    esac
done

inputdir=$1
outputdir=$2

if [ ! -n "$inputdir" ] ; then 
echo ERROR: the input directory should be specified
exit 1
fi

if [ ! -d "$inputdir" ] ; then 
echo ERROR: the input directory does not exist
exit 1
fi

if [ ! -r "$inputdir" ] ; then 
echo ERROR: cannot read from the input directory
exit 1
fi

if [ ! -n "$outputdir" ] ; then 
echo ERROR: the output directory should be specified
exit 1
fi

if [ ! -d "$outputdir" ] ; then 
mkdir -p $outputdir || ( echo ERROR: the output directory does not exist ; exit 1 )
fi

if [ ! -w "$outputdir" ] ; then 
echo ERROR: cannot write to the output directory
exit 1
fi

subject=`basename $2`
experiment=`basename $1`

echo "inputdir      =" $inputdir
echo "outputdir     =" $outputdir
echo "subject       =" $subject
echo "experiment    =" $experiment

echo "raw           =" $raw
echo "anatomy       =" $anatomy
echo "eprime        =" $eprime
echo "datacheck     =" $datacheck
echo "baddata       =" $baddata
echo "icaclass      =" $icaclass
echo "rmegpreproc   =" $rmegpreproc
echo "tmegpreproc   =" $tmegpreproc
echo "eravg         =" $eravg
echo "tfavg         =" $tfavg
echo "powavg        =" $powavg
echo "icamne        =" $icamne
echo "icablpenv     =" $icablpenv
echo "icablpcorr    =" $icablpcorr
echo "icaimagcoh    =" $icaimagcoh
echo "srcavglcmv    =" $srcavglcmv
echo "srcavgdics    =" $srcavgdics    

echo "noise         =" $noise
echo "restin        =" $restin
echo "wrkmem        =" $wrkmem
echo "storym        =" $storym
echo "motort        =" $motort

###################################################################################################
# create all possible output directories
# once all data files have been copied, all directories that are still empty will be removed
echo creating directories ...

mkdir -p $outputdir/release-notes

# this relies on the convention of specific experiments mapping onto scan numbers
name[1]='Rnoise'
name[2]='Pnoise'
name[3]='Restin'
name[4]='Restin'
name[5]='Restin'
name[6]='Wrkmem'
name[7]='Wrkmem'
name[8]='StoryM'
name[9]='StoryM'
name[10]='Motort'
name[11]='Motort'

for number in 1 2 3 4 5 6 7 8 9 10 11 ; do
mkdir -p $outputdir/unprocessed/MEG/${number}-${name[$number]}/4D/
mkdir -p $outputdir/unprocessed/MEG/${number}-${name[$number]}/EPRIME/
done

mkdir -p $outputdir/MEG/anatomy/provenance
mkdir -p $outputdir/MEG/anatomy/figures/provenance

for scan in Rnoise Pnoise Restin Wrkmem StoryM Motort; do
for pipeline in datacheck baddata icaclass rmegpreproc tmegpreproc eravg tfavg powavg icamne icablpenv icablpcorr icaimagcoh srcavglcmv srcavgdics; do
mkdir -p $outputdir/MEG/$scan/$pipeline/provenance
mkdir -p $outputdir/MEG/$scan/$pipeline/figures/provenance
done
done

###################################################################################################
echo writing release-notes ...

rm -f $outputdir/release-notes/MEG.txt

date                                              >> $outputdir/release-notes/MEG.txt
cat << EOF                                        >> $outputdir/release-notes/MEG.txt

The following data processing pipelines were executed with megconnectome version 2.0 or 2.1
-datacheck
-baddata
-icaclass
-rmegpreproc
-powavg
-tmegpreproc
-tfavg
-eravg

The following data processing pipelines were executed with megconnectome version 2.2
-anatomy
-icamne
-icablpenv
-icablpcorr
-icaimagcoh
-srcavglcmv
-srcavgdics

EOF

cat /HCP/scratch/meg/release/notes/${experiment}  >> $outputdir/release-notes/MEG.txt
cat << EOF                                        >> $outputdir/release-notes/MEG.txt

These data were generated and made available by the Human Connectome Project, WU-Minn Consortium (Principal Investigators: David Van Essen and Kamil Ugurbil; 1U54MH091657), which is funded by the 16 NIH Institutes and Centers that support the NIH Blueprint for Neuroscience Research and by the McDonnell Center for Systems Neuroscience at Washington University.

For additional information on how to acknowledge HCP and cite HCP publications if you have used data provided by the WU-Minn HCP consortium, see http://www.humanconnectome.org/citations.

As a reminder, users of these datasets must comply with the Data Use Terms that were agreed upon when receiving these data.

EOF

###################################################################################################
echo packaging raw data ...

if [ $raw = true ] ; then 

if [ $noise = true ] ; then 
for number in 1 2 ; do
for file in config  c,rfDC  e,rfhp1.0Hz,COH  e,rfhp1.0Hz,COH1; do 
copy $inputdir/SCANS/${number}/4D/"$file" $outputdir/unprocessed/MEG/${number}-${name[$number]}/4D/
done
done
fi

if [ $restin = true ] ; then 
for number in 3 4 5 ; do
for file in config  c,rfDC  e,rfhp1.0Hz,COH  e,rfhp1.0Hz,COH1; do 
copy $inputdir/SCANS/${number}/4D/"$file" $outputdir/unprocessed/MEG/${number}-${name[$number]}/4D/
done
done
fi

if [ $wrkmem = true ] ; then 
for number in 6 7 ; do
for file in config  c,rfDC  e,rfhp1.0Hz,COH  e,rfhp1.0Hz,COH1; do 
copy $inputdir/SCANS/${number}/4D/"$file" $outputdir/unprocessed/MEG/${number}-${name[$number]}/4D/
done
if [ $eprime = true ] ; then 
copy $inputdir/SCANS/${number}/LINKED_DATA/EPRIME/*.xlsx $outputdir/unprocessed/MEG/${number}-${name[$number]}/EPRIME/
copy $inputdir/SCANS/${number}/LINKED_DATA/EPRIME/*.tab  $outputdir/unprocessed/MEG/${number}-${name[$number]}/EPRIME/
fi
done
fi

if [ $storym = true ] ; then 
for number in 8 9 ; do
for file in config  c,rfDC  e,rfhp1.0Hz,COH  e,rfhp1.0Hz,COH1; do 
copy $inputdir/SCANS/${number}/4D/"$file" $outputdir/unprocessed/MEG/${number}-${name[$number]}/4D/
done
if [ $eprime = true ] ; then 
copy $inputdir/SCANS/${number}/LINKED_DATA/EPRIME/*.xlsx $outputdir/unprocessed/MEG/${number}-${name[$number]}/EPRIME/
copy $inputdir/SCANS/${number}/LINKED_DATA/EPRIME/*.tab  $outputdir/unprocessed/MEG/${number}-${name[$number]}/EPRIME/
fi
done
fi

if [ $motort = true ] ; then 
for number in 10 11 ; do
for file in config  c,rfDC  e,rfhp1.0Hz,COH  e,rfhp1.0Hz,COH1; do 
copy $inputdir/SCANS/${number}/4D/"$file" $outputdir/unprocessed/MEG/${number}-${name[$number]}/4D/
done
if [ $eprime = true ] ; then 
copy $inputdir/SCANS/${number}/LINKED_DATA/EPRIME/*.xlsx $outputdir/unprocessed/MEG/${number}-${name[$number]}/EPRIME/
copy $inputdir/SCANS/${number}/LINKED_DATA/EPRIME/*.tab  $outputdir/unprocessed/MEG/${number}-${name[$number]}/EPRIME/
fi
done
fi

fi

###################################################################################################
echo packaging anatomy ...

if [ $anatomy = true ] ; then 

copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_transform.txt        $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_fiducials.txt        $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_landmarks.txt        $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_headmodel.mat        $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_2d.mat $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d4mm.mat $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d6mm.mat $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d8mm.mat $outputdir/MEG/anatomy

copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_transform.txt.xml        $outputdir/MEG/anatomy/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_fiducials.txt.xml        $outputdir/MEG/anatomy/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_landmarks.txt.xml        $outputdir/MEG/anatomy/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_headmodel.mat.xml        $outputdir/MEG/anatomy/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_2d.mat.xml $outputdir/MEG/anatomy/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d4mm.mat.xml $outputdir/MEG/anatomy/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d6mm.mat.xml $outputdir/MEG/anatomy/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d8mm.mat.xml $outputdir/MEG/anatomy/provenance

copy $inputdir/RESOURCES/anatomy/${subject}.L.inflated.4k_fs_LR.surf.gii     $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}.R.inflated.4k_fs_LR.surf.gii     $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}.L.midthickness.4k_fs_LR.surf.gii $outputdir/MEG/anatomy
copy $inputdir/RESOURCES/anatomy/${subject}.R.midthickness.4k_fs_LR.surf.gii $outputdir/MEG/anatomy

for ext in png ; do

copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_headmodel.$ext         $outputdir/MEG/anatomy/figures
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_2d.$ext $outputdir/MEG/anatomy/figures
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d4mm.$ext $outputdir/MEG/anatomy/figures
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d6mm.$ext $outputdir/MEG/anatomy/figures
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d8mm.$ext $outputdir/MEG/anatomy/figures

copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_headmodel.$ext.xml         $outputdir/MEG/anatomy/figures/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_2d.$ext.xml $outputdir/MEG/anatomy/figures/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d4mm.$ext.xml $outputdir/MEG/anatomy/figures/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d6mm.$ext.xml $outputdir/MEG/anatomy/figures/provenance
copy $inputdir/RESOURCES/anatomy/${subject}_MEG_anatomy_sourcemodel_3d8mm.$ext.xml $outputdir/MEG/anatomy/figures/provenance

done

fi

###################################################################################################
# copy txt/mat/png files

pipeline=()
if [ $datacheck   = true ]; then pipeline=("${pipeline[@]}" datacheck  ) ; fi
if [ $baddata     = true ]; then pipeline=("${pipeline[@]}" baddata    ) ; fi
if [ $icaclass    = true ]; then pipeline=("${pipeline[@]}" icaclass   ) ; fi
if [ $rmegpreproc = true ]; then pipeline=("${pipeline[@]}" rmegpreproc) ; fi
if [ $tmegpreproc = true ]; then pipeline=("${pipeline[@]}" tmegpreproc) ; fi
if [ $eravg       = true ]; then pipeline=("${pipeline[@]}" eravg      ) ; fi
if [ $tfavg       = true ]; then pipeline=("${pipeline[@]}" tfavg      ) ; fi
if [ $powavg      = true ]; then pipeline=("${pipeline[@]}" powavg     ) ; fi

scan=()
if [ $noise     = true ]; then scan=("${scan[@]}" Rnoise Pnoise) ; fi
if [ $restin    = true ]; then scan=("${scan[@]}" Restin       ) ; fi
if [ $wrkmem    = true ]; then scan=("${scan[@]}" Wrkmem       ) ; fi
if [ $storym    = true ]; then scan=("${scan[@]}" StoryM       ) ; fi
if [ $motort    = true ]; then scan=("${scan[@]}" Motort       ) ; fi

echo packaging processed data for ${pipeline[*]} ${scan[*]} ...

for s in "${scan[@]}" ; do
for p in "${pipeline[@]}" ; do

for ext in txt mat ; do
for file in $inputdir/RESOURCES/analysis/*"$s"_"$p"*."$ext" ; do 
copy "$file" $outputdir/MEG/$s/$p/
done
for file in $inputdir/RESOURCES/analysis/*"$s"_"$p"*."$ext".xml ; do 
copy "$file" $outputdir/MEG/$s/$p/provenance
done
done

for ext in png ; do
for file in $inputdir/RESOURCES/analysis/*"$s"_"$p"*."$ext" ; do 
copy "$file" $outputdir/MEG/$s/$p/figures/
done
for file in $inputdir/RESOURCES/analysis/*"$s"_"$p"*."$ext".xml ; do 
copy "$file" $outputdir/MEG/$s/$p/figures/provenance/
done
done

done
done

###################################################################################################
# copy nii files

pipeline=()
if [ $icamne     = true ]; then pipeline=("${pipeline[@]}" icamne    ) ; fi
if [ $icablpenv  = true ]; then pipeline=("${pipeline[@]}" icablpenv ) ; fi
if [ $icablpcorr = true ]; then pipeline=("${pipeline[@]}" icablpcorr) ; fi
if [ $icaimagcoh = true ]; then pipeline=("${pipeline[@]}" icaimagcoh) ; fi
if [ $srcavglcmv = true ]; then pipeline=("${pipeline[@]}" srcavglcmv) ; fi
if [ $srcavgdics = true ]; then pipeline=("${pipeline[@]}" srcavgdics) ; fi

scan=()
if [ $restin    = true ]; then scan=("${scan[@]}" Restin       ) ; fi
if [ $wrkmem    = true ]; then scan=("${scan[@]}" Wrkmem       ) ; fi
if [ $storym    = true ]; then scan=("${scan[@]}" StoryM       ) ; fi
if [ $motort    = true ]; then scan=("${scan[@]}" Motort       ) ; fi

echo packaging source reconstructed data for ${pipeline[*]} ${scan[*]} ...

for s in "${scan[@]}" ; do
for p in "${pipeline[@]}" ; do

for ext in nii ; do
for file in $inputdir/RESOURCES/?meg/*"$s"_"$p"*."$ext" ; do 
copy "$file" $outputdir/MEG/$s/$p/
done
for file in $inputdir/RESOURCES/?meg/*"$s"_"$p"*."$ext".xml ; do 
copy "$file" $outputdir/MEG/$s/$p/provenance
done
done

for ext in png ; do
for file in $inputdir/RESOURCES/?meg/*"$s"_"$p"*."$ext" ; do 
copy "$file" $outputdir/MEG/$s/$p/figures/
done
for file in $inputdir/RESOURCES/?meg/*"$s"_"$p"*."$ext".xml ; do 
copy "$file" $outputdir/MEG/$s/$p/figures/provenance/
done
done

done
done


###################################################################################################
echo cleaning up empty directories ...

find $outputdir -depth -type d -exec rmdir --ignore-fail-on-non-empty {} \;

