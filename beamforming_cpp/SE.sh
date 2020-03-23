#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit
fi

FOLDER_IN=$1
FOLDER_OUT=${FOLDER_IN}_6c_123567

if [ ! -d $FOLDER_OUT ]; then
    mkdir $FOLDER_OUT
fi

for file in ${FOLDER_IN}/*.wav
do
    file_name="$(basename -- $file)"
    ./beamforming ${FOLDER_IN}/${file_name} ${FOLDER_OUT}/${file_name}
done
