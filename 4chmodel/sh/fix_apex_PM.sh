#!/bin/bash

clear


CASE=""

usage(){
    echo "Script to extract the biventricular surfaces from a healthy case."
    echo "Parameters:"
    echo "-c/--case: Case number with 2 digits."
    echo "-h/--help: Parameters usage."
}

while [[ "$#" -gt 0 ]]; do # While there are options
    case "$1" in
        # Cases with the input as -option x
        -c) CASE="$2"; shift 2;;
        -h) usage; exit;;

        # Cases with the input as --option=x
        --case=*) CASE="${1#*=}"; shift 1;;
        --help) usage; exit;;

        -*) echo "Command $1" >&2 " not found. Use -h or --help for more info."; exit 1;; # Command not found
        *) handle_argument "$1"; shift 1;; # For spaces and stuff
    esac
done

# Mandatory option:
#----------------------------------------------
if [ $CASE = "" ] ; then
    echo "Arguments missing or incorrect."
    usage
    exit 1
fi
#----------------------------------------------
path24ch="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$CASE/meshing/1000um"
path2biv=$path24ch"/BiV"
cmd="cp /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$CASE/meshing/1000um/BiV/UVC_PM/model/biv/BiV.biv.uvcapex.vtx  /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$CASE/meshing/1000um/BiV/UVC_PM/model/biv/BiV.biv.uvcapex_original.vtx"

echo $cmd
eval $cmd

cmd="cp /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$CASE/meshing/1000um/BiV/UVC/model/biv/BiV.biv.uvcapex.vtx  /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$CASE/meshing/1000um/BiV/UVC_PM/model/biv/BiV.biv.uvcapex.vtx"
echo $cmd
eval $cmd

cmd="rm -rf /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$CASE/meshing/1000um/BiV/UVC_PM/sols"
echo $cmd
eval $cmd

cmd="rm -rf /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$CASE/meshing/1000um/BiV/UVC_PM/UVC"
echo $cmd
eval $cmd

cmd="/usr/bin/python2 /home/crg17/Desktop/scripts/4chmodel/Python/CARP/model_arch_ek.py --overwrite-behaviour overwrite --ID=$path2biv/UVC_PM --basename $path2biv/BiV --uvc --mode biv --np 20 --tags=/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/forall/etags.sh"

echo $cmd
eval $cmd

UVC_PATH2BIV="$path2biv/UVC_PM"

cmd="/home/common/bin/GlVTKConvert -m $path2biv/BiV -n $UVC_PATH2BIV/UVC/COORDS_Z.dat -n $UVC_PATH2BIV/UVC/COORDS_V.dat -n $UVC_PATH2BIV/UVC/COORDS_RHO.dat -n $UVC_PATH2BIV/UVC/COORDS_PHI.dat -o $path2biv/UVC_PM/UVC/BiV_UVC_PM" 

echo $cmd
eval $cmd