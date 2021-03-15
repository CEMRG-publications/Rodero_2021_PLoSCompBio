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

# The script starts:
# --------------------------------------------------------------------------------------------------------------------

# <---- From config-heart.yaml

path24ch="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$CASE/meshing/1000um"
path2biv=$path24ch"/BiV"
meshtool_path="/home/common/meshtool_old/meshtool"

cmd="cp "$path24ch"/FEC/h_case"$CASE"_FEC_w0_h0.elem "$path24ch"/h_case"$CASE".elem"
echo $cmd
eval $cmd

cmd="cp "$path2biv"/FEC/BiV_FEC_w0_h0.elem "$path2biv"/BiV.elem"
echo $cmd
eval $cmd

MESHFILE="h_case"$CASE

OUTFILE_1="Ao_RV_base"
OUTFILE_2="AV_base"

OP="2:5/1:9"


cmd="$meshtool_path extract surface -msh=$path24ch/$MESHFILE
               -surf=$path24ch/$OUTFILE_1,$path24ch/$OUTFILE_2
               -ofmt=vtk
               -op=$OP
               "

echo $cmd
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/cpp/bin/Main_merge_4PM.o "$CASE" PM"

echo $cmd
eval $cmd

cmd=$meshtool_path" map -submsh="$path2biv"/BiV -files="$path24ch"/PM_base.surf.vtx -outdir="$path2biv

echo $cmd
eval $cmd

cmd=$meshtool_path" map -submsh="$path2biv"/BiV -files="$path24ch"/PM_base.surf -outdir="$path2biv

echo $cmd
eval $cmd

cmd=$meshtool_path" map -submsh="$path2biv"/BiV -files="$path24ch"/PM_base.surf.neubc -outdir="$path2biv

echo $cmd
eval $cmd

cmd="mv "$path2biv"/PM_base.surf.vtx "$path2biv"/BiV.PMbase.surf.vtx; mv "$path2biv"/PM_base.surf "$path2biv"/BiV.PMbase.surf; mv "$path2biv"/PM_base.surf.neubc "$path2biv"/BiV.PMbase.surf.neubc"

echo $cmd
eval $cmd

cmd="cp "$path2biv"/BiV.PMbase.surf.vtx "$path2biv"/BiV.base.surf.vtx; mv "$path2biv"/BiV.PMbase.surf "$path2biv"/BiV.base.surf; mv "$path2biv"/BiV.PMbase.surf.neubc "$path2biv"/BiV.base.surf.neubc"

echo $cmd
eval $cmd

cmd="/usr/bin/python2 /home/crg17/Desktop/scripts/4chmodel/Python/CARP/model_arch_ek.py --overwrite-behaviour overwrite --ID=$path2biv/UVC_PM --basename $path2biv/BiV --uvc --mode biv --np 20 --tags=/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/forall/etags.sh"

echo $cmd
eval $cmd

UVC_PATH2BIV="$path2biv/UVC_PM"

cmd="/home/common/bin/GlVTKConvert -m $path2biv/BiV -n $UVC_PATH2BIV/UVC/COORDS_Z.dat -n $UVC_PATH2BIV/UVC/COORDS_V.dat -n $UVC_PATH2BIV/UVC/COORDS_RHO.dat -n $UVC_PATH2BIV/UVC/COORDS_PHI.dat -o $path2biv/UVC_PM/UVC/BiV_UVC_PM" 

echo $cmd
eval $cmd


