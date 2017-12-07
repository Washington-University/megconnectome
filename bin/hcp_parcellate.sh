#!/bin/sh

# This script parcellates all dense cifti files in a given directory,
# according to a given parcellation scheme. Optionally the output 
# directory is different from the input directory.
#
# Usage:
#   hcp_parcellate <labelfile> <in-dir> <out-dir> <out-suffix>
#
# <labelfile> points to a dense-label file
# <in-dir>    is the input directory
# <out-dir>   is the output directory
# <out-suffix> is the suffix used in the output filename, to distinguish 
# different parcellation schemes
#
# This script should be run post-packaging, i.e. on a *.zip
# directory.

LABELFILE=$1
DATADIR_IN=$2
DATADIR_OUT=$3
OUTSUFFIX=$4

# deal with the dconn 
LIST=`find $DATADIR_IN -name '*dconn.nii'`
if [ -n "$LIST" ]
then
for f in $LIST
  do
    # the files are located a few directories lower than in $DATADIR_IN
    datadir_in=`dirname "$f"`
    file_in=`basename "$f"`
    
    datadir_out="${datadir_in//$DATADIR_IN/$DATADIR_OUT}"
    file_out="$datadir_out"/"${file_in//dconn/$OUTSUFFIX.pconn}"

    if [ ! -d $datadir_out ]
    then
      echo creating $datadir_out
      mkdir -p $datadir_out
    fi

    echo creating $file_out

    FILE_OUT1=$datadir_out/tempfile.nii
    FILE_OUT2=$file_out
    
    wb_command -cifti-parcellate $f $LABELFILE ROW $FILE_OUT1
    wb_command -cifti-parcellate $FILE_OUT1 $LABELFILE COLUMN $FILE_OUT2
    rm $FILE_OUT1

    # copy the dense label file that has been used into the output directory
    labelfile_out="$datadir_out"/`basename "$LABELFILE"`
    cp $LABELFILE $labelfile_out

  done
fi

# dieal with the dtseries, which may also have dscalar files in the dtseries dir
LIST=`find $DATADIR_IN -name '*dtseries.nii'`
if [ -n "$LIST" ]
then
for f in $LIST 
  do
    # the files are located a few directories lower than in $DATADIR_IN
    datadir_in=`dirname "$f"`
    file_in=`basename "$f"`
    
    datadir_out="${datadir_in//$DATADIR_IN/$DATADIR_OUT}"
    file_out="$datadir_out"/"${file_in//dtseries/$OUTSUFFIX.ptseries}"
    file_out="${file_out//dscalar/$OUTSUFFIX.pscalar}"

    if [ ! -d $datadir_out ]
    then
      echo creating $datadir_out
      mkdir -p $datadir_out
    fi

    echo creating $file_out

    FILE_OUT=$file_out
    
    wb_command -cifti-parcellate $f $LABELFILE COLUMN $FILE_OUT
    
    # copy the dense label file that has been used into the output directory
    labelfile_out="$datadir_out"/`basename "$LABELFILE"`
    cp $LABELFILE $labelfile_out

  done
fi

# deal with the dscalar
LIST=`find $DATADIR_IN -name '*dscalar.nii'`
if [ -n "$LIST" ]
then
for f in $LIST
  do
    # the files are located a few directories lower than in $DATADIR_IN
    datadir_in=`dirname "$f"`
    file_in=`basename "$f"`
    
    datadir_out="${datadir_in//$DATADIR_IN/$DATADIR_OUT}"
    file_out="$datadir_out"/"${file_in//dtseries/$OUTSUFFIX.ptseries}"
    file_out="${file_out//dscalar/$OUTSUFFIX.pscalar}"

    if [ ! -d $datadir_out ]
    then
      echo creating $datadir_out
      mkdir -p $datadir_out
    fi

    echo creating $file_out

    FILE_OUT=$file_out
    
    wb_command -cifti-parcellate $f $LABELFILE COLUMN $FILE_OUT
    
    # copy the dense label file that has been used into the output directory
    labelfile_out="$datadir_out"/`basename "$LABELFILE"`
    cp $LABELFILE $labelfile_out
  
  done
fi

