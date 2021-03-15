#/bin/bash

clear

current_case=""
jobname=""
time_days=0
time_hours="00"
time_minutes="00"
time_seconds="00"
ncores=64
nnodes=8
simtype=""
alpha_method=0
mech_mass_damping=0.1
mech_stiffness_damping=0.1
mech_rho_inf=0.0

function Usage(){
    echo "Script to run the simulations with exactly the same parameters as defined in the template files on Tom2."
    echo "Parameters:"
    echo "--alpha_method: Activate or desactivate the alpha method. Default is desactivated."
    echo "-c/--case: Case number with 2 digits."
    echo "-j/--jobname: Name of the job submitted."
    echo "-mdamp/--mech_mass_damping: mechanical mass damping factor in 1 / ms"
    echo "-sdamp/--mech_stiffness_damping: mechanical stiffness damping factor in ms."
    echo "-ncores/--num_cores: Number of cores per node to use. Default and maximum is 64."
    echo "-nnodes/--num_nodes: Number of nodes to use. Default and maximum is 8."
    echo "-rhoinf/--mech_rho_inf: Spectral radius rho_inf for the generalized alpha time integrator"
    echo "-simtype/--simulation: Type of simulation. Options are unloading or wk3."
    echo "-td/--time_days: Days (with one digit) in wall time when the simulation will be killed. Default is 0."
    echo "-th/--time_hours: Hours (with two digits) in wall time when the simulation will be killed. Default is 00."
    echo "-tm/--time_minutes: Hours (with two digits) in wall time when the simulation will be killed. Default is 00."
    echo "-ts/--time_seconds: Seconds (with two digits) in wall time when the simulation will be killed. Default is 00."
    echo "-h/--help: Parameters usage."
}

function Warnings(){
    if [[ $current_case == "" ]]; then
    echo "Case number missing (-c/--case)."
    exit 1
    fi

    if [[ $simtype -ne "unloading" ]]; then
        if [[ $simtype -ne "wk3" ]]; then
        echo "Invalid simulation type. Options are unloading or wk3."
        exit 1
        fi
    fi

    nchar="${#current_case}"
    if [ $nchar -lt 2 ]; then
    echo "Case number (-c/--case) needed with at least two digits."
    exit 1
    fi

    nchar="${#time_days}"
    if [ $nchar -ne 1 ]; then
    echo "Days in wall time (-td/--time_days) needed with one digit."
    exit 1
    fi

    nchar="${#time_hours}"
    if [ $nchar -ne 2 ]; then
    echo "Hours in wall time (-th/--time_hours) needed with two digits."
    exit 1
    fi

    nchar="${#time_minutes}"
    if [ $nchar -ne 2 ]; then
    echo "Minutes in wall time (-tm/--time_minutes) needed with two digits."
    exit 1
    fi

    nchar="${#time_seconds}"
    if [ $nchar -ne 2 ]; then
    echo "Seconds in wall time (-ts/--time_seconds) needed with two digits."
    exit 1
    fi

    if [ "$time_hours" -ge 24 ]; then
    echo "Hours in wall time should be less than 24."
    exit 1
    fi

    if [ "$time_minutes" -ge 60 ]; then
    echo "Minutes in wall time should be less than 60."
    exit 1
    fi

    if [ "$time_seconds" -ge 60 ]; then
    echo "Seconds in wall time should be less than 60."
    exit 1
    fi

    if [ "$ncores" -gt 64 ]; then
    echo "Number of cores must not be greater than 64."
    exit 1
    fi

    if [ "$nnodes" -gt 8 ]; then
    echo "Number of nodes must not be greater than 8."
    exit 1
    fi
    
    }

function Substitute_line {
    cmd="sed -i '$1"s"/.*/$2/' $simulation_path"/"$simulation_name"
    eval $cmd
}


while [ "$1" != "" ]; do
    case $1 in
        --alpha_method )        shift
                                alpha_method=$1
                                ;;
        -c | --case )           shift
                                current_case=$1
                                ;;
        -j | --jobname )        shift
                                jobname=$1
                                ;;
        -mdamp | --mech_mass_damping )  shift
                                        mech_mass_damping=$1
                                        ;;
        -sdamp | --mech_stiffness_damping )     shift
                                                mech_stiffness_damping=$1
                                                ;;
       -ncores | --num_cores )  shift
                                ncores=$1
                                ;;
       -nnodes | --num_nodes )  shift
                                nnodes=$1
                                ;;
        -rhoinf | --mech_rho_inf )  shift
                                    mech_rho_inf=$1
                                    ;;
       -simtype | --simulation ) shift
                                 simtype=$1
                                 ;;
       -td | --time_days )      shift
                                time_days=$1
                                ;;
        -th | --time_hours )    shift
                                time_hours=$1
                                ;;
        -tm | --time_minutes )  shift
                                time_minutes=$1
                                ;;
        -ts | --time_seconds )  shift
                                time_seconds=$1
                                ;;
        -h | --help )           Usage
                                exit
                                ;;
        * )                     echo "Command not found. Use -h or --help for more info."
                                exit 1
    esac
    shift
done



if [[ $jobname == "" ]]
    then
    jobname=$current_case"_"$simtype
fi

Warnings

template_path="/scratch/crg17/scripts"
template_name="template_"$simtype".in"
simulation_path="/scratch/crg17/simulations"
simulation_name=$jobname".in"

# We copy the template to a new file substituting the line we want

cmd="cp "$template_path"/"$template_name" "$simulation_path"/"$simulation_name
eval $cmd

cmd="chmod +w "$simulation_path"/"$simulation_name
eval $cmd

# Substitute JOBNAME
line_num=4
line_str="#SBATCH -J "$jobname

Substitute_line $line_num "$line_str"

line_num=19
line_str="JOBNAME="$jobname

Substitute_line $line_num "$line_str"

# Substitute WALLTIME

line_num=5
line_str="#SBATCH -t "$time_days"-"$time_hours":"$time_minutes":"$time_seconds

Substitute_line $line_num "$line_str"

# Substitute --nodes

line_num=6
line_str="#SBATCH --nodes="$nnodes

Substitute_line $line_num "$line_str"

# Substitute --ntasks-per-node

line_num=7
line_str="#SBATCH --ntasks-per-node="$ncores

Substitute_line $line_num "$line_str"

# Substitute NPROC

NPROC=$((ncores*nnodes))

line_num=11
line_str="NPROC="$NPROC

Substitute_line $line_num "$line_str"

# Substitute case_num

line_num=17
line_str="case_num="$current_case"HC"

Substitute_line $line_num "$line_str"

# Substitute alpha_method

line_num=18
line_str="alpha_method="$alpha_method

Substitute_line $line_num "$line_str"

# Substitute mech_mass_damping

line_num=20
line_str="mech_mass_damping="$mech_mass_damping

Substitute_line $line_num "$line_str"

# Substitute mech_stiffness_damping

line_num=21
line_str="mech_stiffness_damping="$mech_stiffness_damping

Substitute_line $line_num "$line_str"

# Substitute mech_rho_inf

line_num=22
line_str="mech_rho_inf="$mech_rho_inf

Substitute_line $line_num "$line_str"
