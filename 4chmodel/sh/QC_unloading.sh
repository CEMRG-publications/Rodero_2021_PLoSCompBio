#!/bin/bash

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
#path2mesh="/media/crg17/Seagate\ Expansion\ Drive/h_case"$current_case"/simulations"
path2mesh="/media/crg17/Seagate\ Backup\ Plus\ Drive/CT_cases/h_case"$current_case"/simulations/cohort"
meshtool="/home/common/meshtool_new/meshtool/meshtool"

cmd="scp -r crg17@login.archer.ac.uk:/work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/"$current_case*" "$path2mesh"/."
echo $cmd
eval $cmd

path2mesh=$path2mesh"/"$current_case"HC_unloading_tpeak_120"

cmd="cp "$path2mesh"/../../../meshing/"$current_case"HC/heart_case"*" "$path2mesh"/."
echo $cmd
eval $cmd

cmd=$meshtool" convert -ifmt=carp_bin -ofmt=carp_txt -imsh="$path2mesh"/heart_case"$current_case"HC_800um -omsh="$path2mesh"/"$current_case"_unloaded"
echo $cmd
eval $cmd

cmd="cp "$path2mesh"/reference.pts "$path2mesh"/"$current_case"_unloaded.pts"
echo $cmd
eval $cmd

cmd=$meshtool" convert -ifmt=carp_txt -ofmt=vtk_bin -omsh="$path2mesh"/"$current_case"_unloaded -imsh="$path2mesh"/"$current_case"_unloaded"
echo $cmd
eval $cmd

cmd=$meshtool" convert -ifmt=carp_txt -ofmt=carp_bin -omsh="$path2mesh"/heart_case"$current_case"HC_800um_unloaded -imsh="$path2mesh"/"$current_case"_unloaded"
echo $cmd
eval $cmd

cmd="scp "$path2mesh"/heart_case"$current_case"HC_800um_unloaded.bpts crg17@login.archer.ac.uk:/work/e348/e348/shared/carp-meshes/h4ckcl/meshes/"$current_case"HC/."
echo $cmd
eval $cmd
