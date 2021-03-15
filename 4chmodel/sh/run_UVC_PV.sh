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

OUTFILE_1="Ao_RV_base_ava"
OUTFILE_2="AV_base_ava"
OUTFILE_3="PV_base_ava"
OUTFILE_4="Ao_base_ava"
OUTFILE_5="LV_MV_base_ava"
OUTFILE_6="RV_TV_base_ava"

OP="2:5/1:9/2:10/1:5/1:7/2:8"


cmd="$meshtool_path extract surface -msh=$path24ch/$MESHFILE
               -surf=$path24ch/$OUTFILE_1,$path24ch/$OUTFILE_2,$path24ch/$OUTFILE_3,$path24ch/$OUTFILE_4,$path24ch/$OUTFILE_5,$path24ch/$OUTFILE_6
               -ofmt=vtk
               -op=$OP
               "

echo $cmd
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/cpp/bin/Main_merge_4PM.o "$CASE" ava"

echo $cmd
eval $cmd

cmd=$meshtool_path" map -submsh="$path2biv"/BiV -files="$path24ch"/PM_base_ava.surf.vtx -outdir="$path2biv

echo $cmd
eval $cmd

cmd=$meshtool_path" map -submsh="$path2biv"/BiV -files="$path24ch"/PM_base_ava.surf -outdir="$path2biv

echo $cmd
eval $cmd

cmd=$meshtool_path" map -submsh="$path2biv"/BiV -files="$path24ch"/PM_base_ava.surf.neubc -outdir="$path2biv

echo $cmd
eval $cmd

cmd="mv "$path2biv"/PM_base_ava.surf.vtx "$path2biv"/BiV.PM_base_ava.surf.vtx; mv "$path2biv"/PM_base_ava.surf "$path2biv"/BiV.PM_base_ava.surf; mv "$path2biv"/PM_base_ava.surf.neubc "$path2biv"/BiV.PM_base_ava.surf.neubc"

echo $cmd
eval $cmd

cmd="cp "$path2biv"/BiV.PM_base_ava.surf.vtx "$path2biv"/BiV.base.surf.vtx; mv "$path2biv"/BiV.PM_base_ava.surf "$path2biv"/BiV.base.surf; mv "$path2biv"/BiV.PM_base_ava.surf.neubc "$path2biv"/BiV.base.surf.neubc"

echo $cmd
eval $cmd

cmd="/usr/bin/python2 /home/crg17/Desktop/scripts/4chmodel/Python/CARP/model_arch_ek.py --overwrite-behaviour delete --ID=$path2biv/UVC_ava --basename $path2biv/BiV --uvc --mode biv --np 20 --tags=/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/forall/etags.sh"

echo $cmd
eval $cmd

UVC_PATH2BIV="$path2biv/UVC_ava"

cmd="/home/common/bin/GlVTKConvert -m $path2biv/BiV -n $UVC_PATH2BIV/UVC/COORDS_Z.dat -n $UVC_PATH2BIV/UVC/COORDS_V.dat -n $UVC_PATH2BIV/UVC/COORDS_RHO.dat -n $UVC_PATH2BIV/UVC/COORDS_PHI.dat -F bin -o $UVC_PATH2BIV/UVC/BiV_UVC_ava" 

echo $cmd
eval $cmd


