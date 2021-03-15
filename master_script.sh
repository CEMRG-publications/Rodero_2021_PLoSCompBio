#/bin/bash

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
path2scripts="/home/crg17/Desktop/scripts"

cmd="cp -r "$path2scripts"/4chmodel "$path2scripts"/4chmodel_"$current_case
echo $cmd
eval $cmd

path2scripts=$path2scripts"/4chmodel_"$current_case

cmd=$path2scripts"/cpp/bin/Main_map_independent.o "$current_case" chambers"
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_map_independent.o "$current_case" valves"
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_map_independent.o "$current_case" triangles"
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_map_independent.o "$current_case" fibres"
echo $cmd
eval $cmd

cmd=$path2scripts"/sh/bivExtraction.sh -c "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_map_biv.o "$current_case
echo $cmd
eval $cmd

cmd="cd /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/meshing/1000um/BiV/"
echo $cmd
eval $cmd

cmd="mkdir -p fibres; cd ./fibres"
echo $cmd
eval $cmd

for fibre_exp in 'apba' 'epi' 'endoLV' 'endoRV'
do 
cmd=$path2scripts"/Python/run_fibers.py --experiment "$fibre_exp" --current_case "$current_case" --np 20 --overwrite-behaviour overwrite"
echo $cmd
eval $cmd
done

cmd=$path2scripts"/sh/fibres.sh -c "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_GlRF_independent.o "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/sh/run_UVC_eikonal.sh -c "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_FEC_independent.o $current_case 5 33 1"
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_FEC_independent.o $current_case 5 70 0"
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_FEC_independent.o $current_case 5 100 0"
echo $cmd
eval $cmd

cmd=$path2scripts"/sh/create_bottom_third.sh -c "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/Python/CARP/run_4ch_healthy.py --electrode_name bottom_third --electrode_path /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/simulations --experiment EP --mesh_name h_case"$current_case" --mesh_path /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/meshing/1000um --TestID "$current_case"_EP --FEC_ratio=7 --np 20 --sim_folder /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/simulations --fast"
echo $cmd
eval $cmd

cmd=$path2scripts"/sh/run_UVC_PM.sh -c "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/Python/eldata_UVC_ek.py "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/Python/create_PM_PV.py "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_premechanics.o "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/Main_premechanics_epi_no_isle.o "$current_case
echo $cmd
eval $cmd


cmd=$path2scripts"/cpp/bin/Main_premechanics_septum_no_isle.o "$current_case
echo $cmd
eval $cmd

cmd=$path2scripts"/cpp/bin/QC_unloading_compiled.o "$current_case" h_case"$current_case"_FEC_w5_h70.elem"
echo $cmd
eval $cmd

cmd="rm -rf "$path2scripts