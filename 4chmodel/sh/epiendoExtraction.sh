#!/bin/bash
# /home/crg17/Desktop/scripts/4chmodel/sh/epiendoExtraction.sh
# The idx are verteces that are NOT in the label we want. Is like the reverse of connected component.

# <---- From surfaceExtraction.sh
# ----> Check endocardia and epicardia and to septum extraction.sh

clear

CASE=""
LVENDO=""
RVENDO=""
LAENDO=""
RAENDO=""
LVEPI=""
RVEPI=""
LAEPI=""
RAEPI=""

usage(){
    echo "Script to separate the endo from the epi."
    echo "Parameters:"
    echo "-c/--case: Case number with 2 digits."
    echo "-LVendo : Point in the endo of the LV. Analogous to RVendo, LAendo, RAendo, LVepi, RVepi, LAepi, RAepi"
    echo "-h/--help: Parameters usage."
}

while [[ "$#" -gt 0 ]]; do # While there are options
    case "$1" in
        # Cases with the input as -option x
        -c) CASE="$2"; shift 2;;
        -LVendo) LVENDO="$2"; shift 2;;
        -RVendo) RVENDO="$2"; shift 2;;
        -LAendo) LAENDO="$2"; shift 2;;
        -RAendo) RAENDO="$2"; shift 2;;
        -LVepi) LVEPI="$2"; shift 2;;
        -RVepi) RVEPI="$2"; shift 2;;
        -LAepi) LAEPI="$2"; shift 2;;
        -RAepi) RAEPI="$2"; shift 2;;
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

if [ -z $CASE ] || [ -z $LVENDO ] || [ -z $RVENDO ] || [ -z $LAENDO ] || [ -z $RAENDO ] || [ -z $LVEPI ] || [ -z $RVEPI ] || [ -z $LAEPI ] || [ -z $RAEPI ] ; then
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

# # --------------------------------------------------------------------------------------------------------------------

# # Extract LV epi and endo
MESHFILE="LV.surfmesh"
OUTFILE="LV_epi_temp.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$LVENDO
               -ifmt=vtk
               -ofmt=vtk
               -submsh=$path$OUTFILE
               "
echo $cmd
eval $cmd


OUTFILE="LV_endo_temp.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$LVEPI
               -ifmt=vtk
               -ofmt=vtk
               -submsh=$path$OUTFILE
               "
echo $cmd
eval $cmd

# --------------------------------------------------------------------------------------------------------------------

# Extract RV epi and endo
MESHFILE="RV.surfmesh"
OUTFILE="RV_epi.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$RVENDO
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

OUTFILE="RV_endo.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$RVEPI
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

# Extract LA epi and endo
MESHFILE="LA.surfmesh"
OUTFILE="LA_epi.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$LAENDO
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

OUTFILE="LA_endo.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$LAEPI
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

# Extract RA epi and endo
MESHFILE="RA.surfmesh"
OUTFILE="RA_epi.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$RAENDO
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


MESHFILE="RA.surfmesh"
OUTFILE="RA_endo.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$RAEPI
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

# /home/crg17/Desktop/scripts/4chmodel/sh/epiendoExtraction.sh

