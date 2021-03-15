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
path2mesh="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/meshing/1000um"

if [ ! -d $path2mesh ]; then
 echo "The directory for that case doesn't exist. Try writing two digits."
 exit 1
fi


outpath="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$current_case"/simulations"

if [ ! -d $outpath ]; then
	cmd="mkdir -p $outpath"
	echo $cmd
	eval $cmd
fi

MESHFILE="h_case"$current_case

cmd="cp "$path2mesh"/FEC/"$MESHFILE"_FEC_w5_h33.elem "$path2mesh"/"$MESHFILE".elem"
echo $cmd
eval $cmd

cmd="/home/common/meshtool/meshtool extract surface -msh="$path2mesh"/"$MESHFILE" -op=25,26 -ifmt=carp_txt -ofmt=carp_txt -surf="$outpath"/bottom_third"
echo $cmd 
eval $cmd

cmd="mv "$outpath"/bottom_third.surf.vtx "$outpath"/x.vtx"
echo $cmd
eval $cmd

cmd="rm "$outpath"/bottom"*
echo $cmd
eval $cmd

cmd="mv "$outpath"/x.vtx "$outpath"/bottom_third.vtx"
echo $cmd
eval $cmd

cmd="cp "$path2mesh"/FEC/"$MESHFILE"_FEC_w5_h70.elem "$path2mesh"/"$MESHFILE".elem"
echo $cmd
eval $cmd