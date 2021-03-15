unload_proc='168'
cycle_proc='480'

#======== Running unloading =====#
heart_num=21
heart=$heart_num"HC"
par_modified='purkinje'
value='on'

#for value in '0.027' '0.033'; 
#do
##### UNLOADING #####
#cmd='python /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/run_cycle.py --np $unload_proc --runtime 24:00:00  --case $heart --experiment unloading --build CPU0819 --$par_modified $value --sim_name LSA_unload_$par_modified"_"$value --mesh_name_opt heart_case$heart"_800um_default"'

#cmd='python /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/run_cycle.py --np $unload_proc --runtime 24:00:00  --case $heart --experiment unloading --build CPU0819 --stiffness_pericardium $value --stiffness_LSPV $value --stiffness_RSPV $value --stiffness_SVC $value --sim_name LSA_unload_$par_modified"_"$value --mesh_name_opt heart_case$heart"_800um_default"'

#cmd='python /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/run_cycle.py --np $unload_proc --runtime 24:00:00  --case $heart --experiment unloading --build CPU0819 --c_atria $value --c_veins $value --sim_name LSA_unload_$par_modified"_"$value --mesh_name_opt heart_case$heart"_800um_default"'

#cmd='python /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/run_cycle.py --np $unload_proc --runtime 24:00:00  --case $heart --experiment unloading --build CPU0819 --act_seq_name LSA_$par_modified"_"$value --sim_name LSA_unload_$par_modified"_"$value --mesh_name_opt heart_case$heart"_800um_"$par_modified"_"$value'

##### CYCLE #####
#cmd='python /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/run_cycle.py --np $cycle_proc --runtime 24:00:00 --checkpoint 100 --case $heart --experiment wk3 --duration 800 --build CPU0819 --$par_modified $value --sim_name LSA_wk3_$par_modified"_"$value --mesh_name_opt LSA_$par_modified"_"$value'

#cmd='python /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/run_cycle.py --np $cycle_proc --runtime 24:00:00 --checkpoint 100 --case $heart --experiment wk3 --duration 800 --build CPU0819 --stiffness_pericardium $value --stiffness_LSPV $value --stiffness_RSPV $value --stiffness_SVC $value --sim_name LSA_wk3_$par_modified"_"$value --mesh_name_opt LSA_$par_modified"_"$value'

#cmd='python /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/run_cycle.py --np $cycle_proc --runtime 24:00:00  --checkpoint 100 --case $heart --experiment wk3 --duration 800 --build CPU0819 --c_atria $value --c_veins $value --sim_name LSA_wk3_$par_modified"_"$value --mesh_name_opt LSA_$par_modified"_"$value'

cmd='python /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/run_cycle.py --np $cycle_proc --runtime 24:00:00  --checkpoint 100 --case $heart --experiment wk3 --duration 800 --build CPU0819 --act_seq_name LSA_$par_modified"_"$value --sim_name LSA_wk3_$par_modified"_"$value --mesh_name_opt LSA_$par_modified"_"$value --alpha_method True'

##### KEEP PLAYING #####

#cmd=$cmd' --restore /work/e348/e348/crg17/benchmarks/mechanics/h4cKCL/'$heart'_wk350_tpeak100/checkpoint.400.0 --suffix restored'

echo $cmd
eval $cmd 
#done
