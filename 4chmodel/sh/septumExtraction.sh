#/bin/bash

#~ The idx are verteces that are NOT in the label we want. Is like the reverse of connected component.
#~ /home/crg17/Desktop/scripts/4chmodel/sh/septumExtraction.sh

# <---- From epiendoExtraction.sh
# -----> Check in Paraview and to valveExtraction.sh

clear

CASE=""
SEPTUMENDO=""
SEPTUMEPI=""
LVENDO=""

usage(){
    echo "Script to separate the endo from the epi."
    echo "Parameters:"
    echo "-c/--case: Case number with 2 digits."
    echo "-LVendo"
    echo "-septumendo"
    echo "-septumepi"
    echo "-h/--help: Parameters usage."
}

while [[ "$#" -gt 0 ]]; do # While there are options
    case "$1" in
        # Cases with the input as -option x
        -c) CASE="$2"; shift 2;;
	-septumendo) SEPTUMENDO="$2"; shift 2;;
	-septumepi) SEPTUMEPI="$2"; shift 2;;
        -LVendo) LVENDO="$2"; shift 2;;

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

if [ -z $CASE ] || [ -z $LVENDO ] || [ -z $SEPTUMENDO ] || [ -z $SEPTUMEPI ]; then
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

## --------------------------------------------------------------------------------------------------------------------

# Extract LV epi
MESHFILE="LV_epi_temp.surfmesh"
OUTFILE="LV_epi.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$SEPTUMEPI
               -ifmt=vtk
               -ofmt=vtk
               -submsh=$path$OUTFILE
               "
echo $cmd
eval $cmd

cmd="$meshtool_path convert -imsh=$path$OUTFILE 
               -ifmt=vtk 
               -omsh=$path$OUTFILE 
               -ofmt=carp_txt
               "
echo $cmd
eval $cmd

# --------------------------------------------------------------------------------------------------------------------

# Extract LV endo
MESHFILE="LV_endo_temp.surfmesh"
OUTFILE="LV_endo.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$SEPTUMENDO
               -ifmt=vtk
               -ofmt=vtk
               -submsh=$path$OUTFILE
               "
echo $cmd
eval $cmd

cmd="$meshtool_path convert -imsh=$path$OUTFILE 
               -ifmt=vtk 
               -omsh=$path$OUTFILE 
               -ofmt=carp_txt
               "
echo $cmd
eval $cmd

# --------------------------------------------------------------------------------------------------------------------

# Extract LV septum
MESHFILE="LV_endo_temp.surfmesh"
OUTFILE="septum.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$LVENDO
               -ifmt=vtk
               -ofmt=vtk
               -submsh=$path$OUTFILE
               "
echo $cmd
eval $cmd

cmd="$meshtool_path convert -imsh=$path$OUTFILE 
               -ifmt=vtk 
               -omsh=$path$OUTFILE 
               -ofmt=carp_txt
               "
echo $cmd
eval $cmd
# /home/crg17/Desktop/scripts/4chmodel/sh/septumExtraction.sh

