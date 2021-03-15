#!/bin/bash

# <---- From config-heart.yaml
#~ /home/crg17/Desktop/scripts/4chmodel/sh/surfaceExtraction.sh
# ----> To epiendoExtraction.sh

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

#path="/media/crg17/Seagate\ Backup\ Plus\ Drive/CT_cases/h_case$CASE/meshing/1000um/"
path="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$CASE"/meshing/1000um/"
meshtool_path="/home/common/meshtool_old/meshtool"

MESHFILE="h_case"$CASE

OUTFILE_1="LV" # LV minus RV, MV and AV
OUTFILE_2="RV" # RV minus LV, TV and PV
OUTFILE_3="LA" # LA minus MV, PVs BCs and PVs
OUTFILE_4="RA" # RA minus TV, VCs BCs and VCs
OUTFILE_5="LV_base" # Intersec. LV with MV
OUTFILE_6="RV_base" # Intersec. RV with TV
OUTFILE_7="Ao_base"
OUTFILE_8="PA_base"

# OP1="1-2,7,9\;2-1,8,10\;3-7,11,12,13,14,15,18,19,20,21,22\;4-8,16,17,23,24"
# OP2="1:7\;2:8\;1:5\;2:10"

OP1="1-2,7,9/2-1,8,10/3-7,11,12,13,14,15,18,19,20,21,22/4-8,16,17,23,24"
OP2="1:7/2:8/1:5/2:10"


cmd="$meshtool_path extract surface -msh=$path$MESHFILE
               -surf=$path$outFolder$OUTFILE_1,$path$outFolder$OUTFILE_2,$path$outFolder$OUTFILE_3,$path$outFolder$OUTFILE_4
               -ifmt=carp_txt
               -ofmt=vtk
               -op=$OP1
               "

echo $cmd
eval $cmd

cmd="$meshtool_path extract surface -msh=$path$MESHFILE
               -surf=$path$outFolder$OUTFILE_5,$path$outFolder$OUTFILE_6,$path$outFolder$OUTFILE_7,$path$outFolder$OUTFILE_8
               -ifmt=carp_txt
               -ofmt=carp_txt
               -op=$OP2
               "

echo $cmd
eval $cmd

