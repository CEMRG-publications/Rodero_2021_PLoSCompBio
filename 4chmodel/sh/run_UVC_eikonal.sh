#!/bin/bash

# bash /home/crg17/Desktop/scripts/4chmodel/sh/run_UVC_eikonal.sh -c 00

# <---- From Main_GlRF.cpp
clear

current_case=""

usage(){
    echo "Script to create a rule-based fibre field in the biventricular mesh."
    echo "Parameters:"
    echo "-c/--case: Case number with 2 digits."
    echo "-h/--help: Parameters usage."
}

while [ "$1" != "" ]; do
    case $1 in
        -c | --case )           shift
                                current_case=$1
								;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo "Command not found. Use -h or --help for more info."
                                exit 1
    esac
    shift
done

if [ -z $current_case ] ; then
    echo "Arguments missing."
    usage
    exit 1
fi

PATH2BIV="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$current_case/meshing/1000um/BiV"

cmd="/usr/bin/python2 /home/crg17/Desktop/scripts/4chmodel/Python/CARP/model_arch_ek.py --ID=$PATH2BIV/UVC --basename $PATH2BIV/BiV --uvc --mode biv --np 20 --tags=/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/forall/etags.sh"

echo $cmd
eval $cmd

UVC_PATH2BIV="$PATH2BIV/UVC"

cmd="/home/common/bin/GlVTKConvert -m $PATH2BIV/BiV -n $UVC_PATH2BIV/UVC/COORDS_Z.dat -n $UVC_PATH2BIV/UVC/COORDS_V.dat -n $UVC_PATH2BIV/UVC/COORDS_RHO.dat -n $UVC_PATH2BIV/UVC/COORDS_PHI.dat -o $PATH2BIV/BiV_UVC" 

echo $cmd
eval $cmd
