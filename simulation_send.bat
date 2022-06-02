#!/bin/bash

#SBATCH --job-name=test_simulation
#SBATCH --mail-user=noz42@kamis.me ##youremail here
#SBATCH --mail-type=NONE

#SBATCH --account=lsa1 
#SBATCH --partition=standard

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=0-41:11:02
#SBATCH --export=ALL
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.log
####


#### PBS -l qos=flux

#if [ -s "$PBS_NODEFILE" ] ; then
#    echo "Running on"
#    uniq -c $PBS_NODEFILE
#fi

#if [ -d "$PBS_O_WORKDIR" ] ; then
#    cd $PBS_O_WORKDIR
#    echo "Running from $PBS_O_WORKDIR"
#fi



#  simulation_loop.sh
#  
#
#  
#
#Modified by bolaji eniwaye starting 01/30/2017
# lsa_flux # michalz_flux
#get current directory
# 11/17/2018 so originally this was the code was the one used to try to show phase correlation but it didnt do what I want it to do. Originally it had zero dc spread but now I adding a dc spread of .5 to see if it changes anything




#this just prints 1 and 2 which have nothing assigned to them
# so the dollar sign just references a variable 
echo $1
echo $2
cwd=$(pwd)
cwd1=${cwd:1:${#cwd}}

#compile code

#g++ -std=c++11 -o "${iter}run.out" "${cwd}/richneuromod_multi_stim_2_mg_nmda_model_book_tot_fix_log_norm_dc_min_noise_variable_syn_net_ruin.cpp"

g++ -std=c++11 -o "${iter}run.out" "${cwd}/anesthetic_sim.cpp"


#g++ -std=c++11 -o run2.out "${cwd}/mpc_v2.cpp"

#qsub -v "gks=$g,xn=$h,xc=$i,CurOffSet=$k,runs=$j,wi=$l,we=$m,dcbaseline=$n" "${cwd}/simulation_loopdistcur_dynrun.pbs"


daydate=$(date '+%d')
datepeice[0]=$(date '+%m') #month
datepeice[1]=$(date '+%d') #day
datepeice[2]=$(date '+%y') #year

timestart[0]=$(date '+%d')
timestart[1]=$(date '+%H')
timestart[2]=$(date '+%M')
timestart[3]=$(date '+%S')



                #make excecution directory
               # newdir=${cwd}"/gks_"${g}"_xn_"${h}"_xc_"${i}"_run_"${j}"_distcur_CurOffSet_"${k}"_wi_"${l}"_we_"${m}"_dcbaseline_"${n}"_DcSpread_"${o}
               
#newdir=${cwd}"/GKS_"${g}"_WE_"${h}"_WI_"${i}"_AE_"${j}"_AI_"${k}"_TauS_I_"${l}"_STIM_"${n}"_RUN_"${m}"_AMPAMULT_"${o}"_dpblk"
#newdir=${cwd}"/GKS_"${GKS}"_WE_"${WE}"_WI_"${WI}"_AE_"${AE}"_AI_"${AI}"_TauS_I_"${TAU}"_STIM_"${STIM}"_RUN_"${RUN}"_AMPAMULT_"${AMPA_MULT}"_dpblk"
#newdir=${cwd}"/GKS_"${GKS}"_WE_"${WE}"_WI_"${WI}"_AE_"${AE}"_AI_"${AI}"_TauS_I_"${TAU}"_STIM_"${STIM}"_RUN_"${RUN}"_AMPAMULT_"${AMPA_MULT}"_dpblk"
newdir=${cwd}"/gks_"${gks}"_pee_"${pee}"_pii_"${pii}"_pie_"${pie}"_pei_"${pei}"_wee_"${wee_mult}"_wii_"${wii_mult}"_wie_"${wie_mult}"_wei_"${wei_mult}"_CIe_"${CIe}"_CIi_"${CIi}"_stim_"${stim}"_w_mult_"${w_mult}"_nmda_strength_wee_"${wee_nmda_mult}"_run_"${run}"_nmda_strength_wei_"${wei_nmda_mult}"_log_mu_"${log_norm_mu}"_log_std_"${log_norm_std}"_noise_freq_"${noise_freq}"_noise_strength_"${noise_strength}"_noise_dur_"${noise_dur}"_net_"${net}"_ruin_time_"${run_time}
#newdir=${cwd}"/gks_"${k}"_pee_"${a}"_pii_"${b}"_pie_"${c}"_pei_"${d}"_wee_"${e}"_wii_"${f}"_wie_"${G[pp]}"_wei_"${h}"_CIe_"${i}"_CIi_"${j}"_stim_"${l}"_w_mult_"${n}"_nmda_strength_"${O[pp]}"_run_"${m}"_nmda_scan"
		 mkdir "$newdir"
                cp  "network_dyn.hpp" "$newdir"
               # cp  "curr.txt" "$newdir"
                #cp  "mcond.txt" "$newdir"                 

 
                #cd to execution directory
                cd "$newdir"
               
                #run code //////// gks probE probI rEEm rIIm wei wie tauCaSlow # been using, .4 .5 , wieghts .1 .1
	       
              # "${cwd}/run1.out"  $g 10.0 $h 0.0 0.25 1.0 0.02 0.1 $i $k # .06 .3 .5 .06 0.00005 0.001300 0.00190 0.00046 $i $k $g $h 
             #  "${cwd}/run1.out"  $g 10.0 20020 $o 1.0 1.0 ${m} ${l} ${i} ${k} ${n} # .06 .3 .5 .06 0.00005 0.001300 0.00190 0.00046 $i $k $g $h
           #  "${cwd}/${iter}run.out"  $g 10.0 20020 $o 1.0 1.0 ${m} ${l} ${i} ${k} ${n}
#"${cwd}/${iter}run.out"  1 1 1 1 ${h} ${i} ${i} ${h} -0.05 -0.1 ${g} ${j} ${k} ${l} 2400 ${n} ${o}
echo "pee pii pie pei wee wii wie wei Ie Ii gks stim w_mult wee_nmda_mult wei_nmda_multi log_norm_mu log_norm_std"

echo "$pee" "$pii" "$pie" "$pei" "$wee_mult" "$wii_mult" "$wie_mult"  "$wei_mult" "$CIe" "$CIi" "$gks" "$stim" "$w_mult" "$wee_nmda_mult" "$wei_nmda_mult" "$log_norm_mu" "$log_norm_std"

echo " new  pee pii pie pei wee wii wie wei Ie Ii gks stim w_mult wee_nmda_mult wei_nmda_mult log_norm_mu log_norm_std" 


echo "${pee}" "${pii}" "${pie}" "${pei}" "${wee_mult}" "${wii_mult}" "${wie_mult}"  "${wei_mult}" "${CIe}" "${CIi}" "${gks}" "${stim}" "${w_mult}" "${wee_nmda_mult}" "${wei_nmda_mult}" "${log_norm_mu}" "${log_norm_std}"

"${cwd}/${iter}run.out"  ${pee} ${pii} ${pie} ${pei} ${wee_mult} ${wii_mult} ${wie_mult} ${wei_mult} ${CIe} ${CIi} ${gks} ${stim} ${w_mult} ${wee_nmda_mult} ${wei_nmda_mult} ${log_norm_mu} ${log_norm_std} ${noise_freq} ${noise_strength} ${noise_dur} ${net} ${run_time}

              # 
              # #  $g .2 .2 .1 .1 $k $k $h $i #$g .2 .2 .1 .1 .1 .1 $h  #probe,probi .2      #"${cwd}/run.out" $g 0.5 0.5 $h $i
               # "${cwd}/run2.out"  "raster_dat.txt"
		cd "$cwd"
                 
               # newdir1=${cwd}"/gks_"${g}"_net_"${h}"_xc_"${i}"_run_"${j}"_selectcur_CurOffSet"{k} #"/AmpaMult_"${g}"_NoiseMag_"${h}"_ExcCur_"${i}"_InhCur_"${k}"_run_"${j}
#                cd "$cwd" && tar -czf ${newdir}.tar.gz "${newdir}"
#                tar -C "$cwd" -cvzf ${newdir1}.tar.gz "${newdir}"
#                tar -zcvf ${newdir}.tar.gz "${newdir}"
               # zip -r ${newdir}.zip ${newdir1}/* # co 2/7
               # rm -rf "$newdir"#

#make a new directory that i can put all of the individual runs in 
#rundir=${cwd}"/rundata/kncl_"${datepeice[0]}"_"${datepeice[1]}"_"${datepeice[2]}

#savedir="/scratch/michalz_flux/eniwbola/"${datepeice[0]}"_"${datepeice[1]}"_"${datepeice[2]}
savedir="/scratch/michalz_root/michalz/eniwbola/"${datepeice[0]}"_"${datepeice[1]}"_"${datepeice[2]}"_one_per_low_input_p_12_less_r_nmda_1_8_prob_0p1_hz_gks_net_"${net}"_"${run_time}"_spike_rand_prob_cor_con_p5_variabil_2_per_s2"
#-----------------------------------------------had this before to make a group day
#rundatadir=${cwd}"/rundata"
#mkdir "$rundatadir"

mkdir "$savedir"

#rundir=${cwd}"/rundata/Stable_bump_learn_"${datepeice[0]}"_"${datepeice[1]}"_"${datepeice[2]}
#mkdir "$rundir"
#mv  ${cwd}"/gks"*  ${rundir}/ #change
mv   ${newdir}    ${savedir}
#----------------------------------------------------------#
#


## rclone stuf it usaully works but messed up one of my sims
#if [ 40 == 45 ]
#then

#print "okay Isync is 1" 

#sleep 2m

#NewgDir="rdrive:Flux/FluxData/AllData/AmpaVar_"${datepeice[0]}"_"${datepeice[1]}"_"${datepeice[2]}
#CurDatDir="rdrive:/Flux/FluxData/CurrentData"

#sleep 20s

#rclone -vv purge ${NewgDir}

#sleep 20s

#rclone -vv purge ${CurDatDir}

#sleep 30s

#rclone -vv mkdir ${CurDatDir}

#sleep 5s

#rclone -vv mkdir ${NewgDir}

#print "newgdir",${NewgDir}
#print "rundir ",${rundir}


#sleep 1m

##rclone sync ${rundir} "rdrive:/Flux/FluxData/CurrentData"


#rclone -vv copy ${rundir}  ${CurDatDir}

#sleep 1m

#print "past2"
##rclone copy  ${rundir} "rdrive:Flux/FluxData/AllData/AmpaVar_"${datepeice[0]}"_"${datepeice[1]}"_"${datepeice[2]}

#rclone -vv copy  ${rundir} ${NewgDir}

#print "past2"

#fi

#get the parameter range in a text file so that i can iterate over the string
#in matlab
#cat paramsrange.txt


#echo -n "AmpaMult " >> "paramsrange.txt"
#for b in "${AM[@]}"
#do
#echo -n "${b} " >> "paramsrange.txt"
#echo -n " " >> "paramsrange.txt"
#done

#echo " " >> "paramsrange.txt"


 
#echo -n "NoiseMag " >> "paramsrange.txt"
#for a in "${NM[@]}"
#do
#echo -n "${a} " >> "paramsrange.txt"
#echo -n " " >> "paramsrange.txt"
#done

#echo " " >> "paramsrange.txt "


#echo -n "ExcCur " >> "paramsrange.txt"

#for c in "${EC[@]}" 
#do
#echo -n "${c} " >> "paramsrange.txt"
#echo -n " " >> "paramsrange.txt"
#done

#echo " " >> "paramsrange.txt"

#echo -n "InhCur " >> "paramsrange.txt"
#for d in "${IC[@]}"
#do
#echo -n "${d}" >> "paramsrange.txt"
#echo -n "  "   >> "paramsrange.txt"
#done



#mv  paramsrange.txt  ${rundir}/


timeend[0]=$(date '+%d')
timeend[1]=$(date '+%H')
timeend[2]=$(date '+%M')
timeend[3]=$(date '+%S')

dura[0]=$(expr ${timeend[0]} - ${timestart[0]})
dura[1]=$(expr ${timeend[1]} - ${timestart[1]})
dura[2]=$(expr ${timeend[2]} - ${timestart[2]})
dura[3]=$(expr ${timeend[3]} - ${timestart[3]})


echo days "${dura[0]}"
echo hours "${dura[1]}"
echo minutes "${dura[2]}"
echo seconds "${dura[3]}"

