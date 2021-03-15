#/bin/bash
#  /home/crg17/Desktop/scripts/4chmodel/sh/bivExtraction.sh
# Use  tr -d '\r' < bivExtraction.sh > bivExtraction_corr.sh; rm bivExtraction.sh; mv bivExtraction_corr.sh bivExtraction.sh; chmod +x bivExtraction.sh the first time
# <--- From pickApex.cpp
# ----> To Main_map.cpp



clear

current_case=""

usage(){
    echo "Script to extract the biventricular mesh from a healthy case. Ventricle tags are expected to be 1 and 2."
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

if [ "$current_case" = "" ]; then
    echo "Arguments missing."
    usage
    exit 1
fi

# --------------------------------------------------------------------------------------------------------------------
path="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/meshing/1000um/"

if [ ! -d $path ]; then
 echo "The directory for that case doesn't exist. Try writing two digits."
 exit 1
fi

if [ ! -d $path:"Biv/" ]; then
	cmd="mkdir -p $path""BiV/"
	echo $cmd
	eval $cmd
fi

outpath="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/meshing/1000um/BiV/"

MESHFILE="h_case"$current_case
#MESHFILE="HF_case"$current_case

OUTFILE="BiV"

#~ We want the ventricles

               
cmd="/home/common/meshtool/meshtool extract mesh -msh=$path$MESHFILE
               -tags=1,2
               -ifmt=carp_txt
               -ofmt=carp_txt
               -submsh=$outpath$OUTFILE
               "
echo $cmd
eval $cmd
