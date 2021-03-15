CASE=""

usage(){
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

path2case="/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case"$CASE
path2mesh=$path2case"/meshing/1000um"
archerpath=$path2case"/meshing/"$CASE"HC"
path2cav=$path2mesh"/cavities"
path2BC=$path2mesh"/BC"

cmd="mkdir -p "$archerpath

echo $cmd
eval $cmd

cmd="cp "$path2case"/simulations/"$CASE"_EP/vm_act_seq.dat "$archerpath"/"$CASE"HC_vm_act_seq.dat"

echo $cmd
eval $cmd

cmd="cp "$path2mesh"/PM_$CASE""_ava.dat "$archerpath"/PM_"$CASE"HC.dat"

echo $cmd
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2BC"/epicardium_ventricles.neubc "$archerpath"/epicardium_ventricles"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2BC"/IVC_BC.neubc "$archerpath"/IVC"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2BC"/SVC_BC.neubc "$archerpath"/SVC"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2BC"/LA_App_BC.neubc "$archerpath"/appendage"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2BC"/LIPV_BC.neubc "$archerpath"/LI_pulm_vein"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2BC"/LSPV_BC.neubc "$archerpath"/LS_pulm_vein"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2BC"/RIPV_BC.neubc "$archerpath"/RI_pulm_vein"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2BC"/RSPV_BC.neubc "$archerpath"/RS_pulm_vein"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2cav"/LV_endo_closed.neubc "$archerpath"/LV_endo"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2cav"/RV_endo_closed.neubc "$archerpath"/RV_endo"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2cav"/LA_endo_closed.neubc "$archerpath"/LA_endo"

echo $cmd 
eval $cmd

cmd="/home/crg17/Desktop/scripts/4chmodel/sh/neubc2surf.sh "$path2cav"/RA_endo_closed.neubc "$archerpath"/RA_endo"

echo $cmd 
eval $cmd


echo "Finished! Now convert manually the mesh to vtk_bin with the correct elem file and one/two .bpts corresponding to the EDV and ESV configurations."