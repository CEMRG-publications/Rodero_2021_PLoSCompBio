#!/bin/bash

if [ "$#" -ne 2 ]
then
  echo "Wrong number of input parameters"
  echo "Expecting one file with .neubc ending (input) and one without the extension (output name)."
  exit
fi

filename=$1
outname=$2
extension="${filename##*.}"
filename="${filename%.*}"


FILE=$filename".neubc"   # neubc file
OUT=$outname".surf"     # surf file

# --------------------------------
echo "Processing $FILE"

awk 'NR==1{ print }' $FILE > $OUT
awk 'NR>1{ print "Tr " $1 " " $2 " " $3 }' $FILE >> $OUT

echo "wrote output to $OUT "
