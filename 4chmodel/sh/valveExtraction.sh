#/bin/bash

#~ The idx are verteces that are NOT in the label we want. Is like the reverse of connected component.
#~ /home/crg17/Desktop/scripts/4chmodel/sh/septumExtraction.sh

# <---- From epiendoExtraction.sh
# -----> To valvePlaneSurface.sh

clear

CASE=""

usage(){
    echo "Script to separate the endo from the epi."
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

if [ -z $CASE ] ; then
    echo "Arguments missing or incorrect."
    usage
    exit 1
fi
#----------------------------------------------

# The script starts:
# --------------------------------------------------------------------------------------------------------------------

current_case=$CASE

#path="/media/crg17/Seagate\ Backup\ Plus\ Drive/CT_cases/h_case"$current_case"/meshing/1000um/"
path="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/meshing/1000um/"
meshtool_path="/home/common/meshtool_old/meshtool"

MESHFILE="h_case"$current_case
OUTFILE_1="MV"
OUTFILE_2="TV"
OUTFILE_3="AV"
OUTFILE_4="PV"
OUTFILE_5="App"
OUTFILE_6="RIPV"
OUTFILE_7="LIPV"
OUTFILE_8="LSPV"
OUTFILE_9="RSPV"
OUTFILE_10="SVC"
OUTFILE_11="IVC"

cmd="$meshtool_path extract surface -msh=$path$MESHFILE
               -surf=$path$OUTFILE_1,$path$OUTFILE_2,$path$OUTFILE_3,$path$OUTFILE_4,$path$OUTFILE_5,$path$OUTFILE_6,$path$OUTFILE_7,$path$OUTFILE_8,$path$OUTFILE_9,$path$OUTFILE_10,$path$OUTFILE_11
               -ifmt=carp_txt
               -ofmt=vtk
               -op=7-1,3/8-2,4/9-1,5/10-2,6/11-3,12,13,14,15,18,19,20,21,22/12-3,11,13,14,15,18,19,20,21,22/13-3,12,11,14,15,18,19,20,21,22/14-3,12,13,11,15,18,19,20,21,22/15-3,12,13,14,11,18,19,20,21,22/16-4,17,23,24/17-4,16,23,24
    "
echo $cmd
eval $cmd
