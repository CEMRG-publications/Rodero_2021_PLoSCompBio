#!/bin/bash

# bash /home/crg17/Desktop/scripts/4chmodel/sh/fibres.sh -c 01

# <---- From run_fibers.py
# ----> To Main_GlRF.cpps


clear

f_endo="80"
f_epi="-60"
s_endo="-65"
s_epi="25"
# Land et al., 2017; Land and Niederer, 2018
current_case=""
HF=""

usage(){
    echo "Script to create a rule-based fibre field in the biventricular mesh."
    echo "Parameters:"
    echo "-c/--case: Case number with 2 digits."
	echo "-f_endo: Fibre direction in the endocardium. Default is 80."
	echo "-f_epi: Fibre direction in the epicardium. Default is -60."
	echo "-s_endo: Sheet direction in the endocardium. Default is -65."
	echo "-s_epi: Sheet direction in the epicardium. Default is 25."
    echo "-h/--help: Parameters usage."
}

while [ "$1" != "" ]; do
    case $1 in
        -c | --case )           shift
                                current_case=$1
								;;
		-f_endo)				shift
								f_endo=$1
								;;
		-f_epi) 				shift
								f_epi=$1
								;;
		-s_endo)				shift
								s_endo=$1
								;;
		-s_epi) 				shift
								f_epi=$1
								;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo "Command not found. Use -h or --help for more info."
                                exit 1
    esac
    shift
done

if [ -z $current_case ] ; then
    echo "Arguments missing."
    usage
    exit 1
fi

BiV_folder="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/meshing/1000um/BiV"

common_name=$BiV_folder"/fibres"


# -----------------------------------------------------------
FILE=$common_name"/apba/phie.igb"
OUTFILE=$common_name"/apba/phie.dat"

cmd="/home/common/bin/igbextract $FILE -o ascii -O $OUTFILE"

echo $cmd
eval $cmd

cmd="sed -i s/'\s'/'\n'/g $OUTFILE"

echo $cmd
eval $cmd

# -----------------------------------------------------------
FILE=$common_name"/epi/phie.igb"
OUTFILE=$common_name"/epi/phie.dat"

cmd="/home/common/bin/igbextract $FILE -o ascii -O $OUTFILE"

echo $cmd
eval $cmd

cmd="sed -i s/'\s'/'\n'/g $OUTFILE"

echo $cmd
eval $cmd

# -----------------------------------------------------------
FILE=$common_name"/endoLV/phie.igb"
OUTFILE=$common_name"/endoLV/phie.dat"

cmd="/home/common/bin/igbextract $FILE -o ascii -O $OUTFILE"

echo $cmd
eval $cmd

cmd="sed -i s/'\s'/'\n'/g $OUTFILE"

echo $cmd
eval $cmd

# -----------------------------------------------------------
FILE=$common_name"/endoRV/phie.igb"
OUTFILE=$common_name"/endoRV/phie.dat"

cmd="/home/common/bin/igbextract $FILE -o ascii -O $OUTFILE"

echo $cmd
eval $cmd

cmd="sed -i s/'\s'/'\n'/g $OUTFILE"

echo $cmd
eval $cmd

# # -----------------------------------------------------------

cmd="/home/common/bin/GlRuleFibers 
			-m "$BiV_folder"/BiV
			--type biv 
			-a "$common_name"/apba/phie.dat
			-e "$common_name"/epi/phie.dat
			-l "$common_name"/endoLV/phie.dat 
			-r "$common_name"/endoRV/phie.dat
			--alpha_endo $f_endo	
			--alpha_epi $f_epi
			--beta_endo $s_endo
			--beta_epi $s_epi
			-o "$BiV_folder"/fibres/fibres_bayer_"$f_epi"_"$f_endo".lon
			"

echo $cmd
eval $cmd



