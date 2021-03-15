#!/bin/bash
#  /home/crg17/Desktop/scripts/4chmodel/sh/valvePlanesSurface.sh

# <---- From valveExtraction.sh
# -----> To Main.cpp

clear
#/bin/bash



clear

CASE=""
MVLV=""
MVLA=""
TVRV=""
TVRA=""
AVLV=""
AVAO=""
PV=""
APP=""
RIPV=""
LIPV=""
LSPV=""
RSPV=""
SVC=""
IVC=""

usage(){
    echo "Script to separate the endo from the epi."
    echo "Parameters:"
    echo "-c/--case: Case number with 2 digits."
    echo "All the seeds neded are: -MVLV -MVLA -TVRV -TVRA -AVLV -AVAo -PV -App -RIPV -LIPV -LSPV -RSPV -SVC -IVC."
    echo "-h/--help: Parameters usage."
}

while [[ "$#" -gt 0 ]]; do # While there are options
    case "$1" in
        # Cases with the input as -option x
        -c) CASE="$2"; shift 2;;
        -MVLV) MVLV="$2"; shift 2;;
        -MVLA) MVLA="$2"; shift 2;;
        -TVRV) TVRV="$2"; shift 2;;
        -TVRA) TVRA="$2"; shift 2;;
        -AVAo) AVAO="$2"; shift 2;;
        -AVLV) AVLV="$2"; shift 2;;
        -PV) PV="$2"; shift 2;;
        -App) APP="$2"; shift 2;;
        -RIPV) RIPV="$2"; shift 2;;
        -RSPV) RSPV="$2"; shift 2;;
        -LSPV) LSPV="$2"; shift 2;;
        -LIPV) LIPV="$2"; shift 2;;        
        -SVC) SVC="$2"; shift 2;;        
        -IVC) IVC="$2"; shift 2;;        
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

if [ -z $CASE ] || [ -z $MVLV ]  || [ -z $MVLA ]  || [ -z $TVRV ]  || [ -z $TVRA ]  || [ -z $AVAO ]  || [ -z $AVLV ]  || [ -z $PV ]  || [ -z $APP ] || [ -z $RIPV ] || [ -z $LIPV ] || [ -z $LSPV ] || [ -z $RSPV ] || [ -z $SVC ] || [ -z $IVC ] ; then
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

# # # --------------------------------------------------------------------------------------------------------------------

# Extract MV surfaces
MESHFILE="MV.surfmesh"
OUTFILE="MV_LV.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$MVLA
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

OUTFILE="MV_LA.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
               -idx=$MVLV
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

 # Extract TV surfaces
MESHFILE="TV.surfmesh"
OUTFILE="TV_RV.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
                -idx=$TVRA
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


OUTFILE="TV_RA.surfmesh"

cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
                -idx=$TVRV
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

#  # --------------------------------------------------------------------------------------------------------------------

# #  Extract AV surface
 MESHFILE="AV.surfmesh"
 OUTFILE="AV_LV.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$AVAO
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

 MESHFILE="AV.surfmesh"
 OUTFILE="AV_Ao.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$AVLV
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


# # #  --------------------------------------------------------------------------------------------------------------------

# #  Extract PV surface
 MESHFILE="PV.surfmesh"
 OUTFILE="PV_RV.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$PV
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

# # #  --------------------------------------------------------------------------------------------------------------------

# #  Extract PVe surfaces

 MESHFILE="App.surfmesh"
 OUTFILE="App_LA.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$APP
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

# # -------------------------
 MESHFILE="RIPV.surfmesh"
 OUTFILE="RIPV_LA.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$RIPV
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

# # -------------------------
 MESHFILE="LIPV.surfmesh"
 OUTFILE="LIPV_LA.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$LIPV
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

 # # -------------------------
 MESHFILE="LSPV.surfmesh"
 OUTFILE="LSPV_LA.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$LSPV
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
  # -------------------------
 MESHFILE="RSPV.surfmesh"
 OUTFILE="RSPV_LA.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$RSPV
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

# #  Extract SVC surface
 MESHFILE="SVC.surfmesh"
 OUTFILE="SVC_RA.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$SVC
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

#  # Extract IVC surface
 MESHFILE="IVC.surfmesh"
 OUTFILE="IVC_RA.surfmesh"

 cmd="$meshtool_path extract unreachable -msh=$path$MESHFILE
             -idx=$IVC
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


#  /home/crg17/Desktop/scripts/4chmodel/sh/valvePlanesSurface.sh

