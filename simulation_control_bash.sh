#!/bin/sh




#get current directory


#echo $1
#echo $2


#compile code
#icpc "${cwd}/network_dyn.cpp" -O2 -o run1.out
#icpc "${cwd}/AMD_v4.cpp" -O2 -o run2.out

A=(0.01) #pconee
B=(0.02) #pconii
C=(0.02) #pconie
D=(0.01)   #pconei
#E=(0.004 0.0045 0.005 0.0055 0.006 0.0065)   #( 0.01 0.04   ) #( 0.005 0.0075 0.001   ) #gsynee 
E=(0.0065 0.006 0.0055 0.005 0.0045 0.004)  # gsynee revers
F=(0.00190) # gsynii
G=(0.0005 0.0008 0.001175 0.0013 0.0015 0.002) #G=(0.001175 ) # gsynie
#G=( 0.002  0.0015   0.001175  0.0010  0.0008 0.0005) #gsynie revers
H=(0.05) # gsynei
I=(1.0) #iE
J=(0.01) # iI
#K=(0.2 0.6 1.0 1.2)  # gksexc  
K=(0.6 0.8 1.0 1.2 1.4) #gks_  
#L=(0.0 1.0 1.4) # stim_val 
L=(1.0) # stim_val 
M=(1) #(1 2 3 4 5) #run
N=(4) #( 2 3 5 7)
#O=(0.00001 0.0005 0.001 0.005 0.01 0.05) #nmda
O=(0.05 0.01 0.005 0.001 0.0005 0.00001) #nmda reverse

PP=(0 1 2 3 4 5)

A=(0.1) #pconee
B=(0.1) #pconii
C=(0.1) #pconie
D=(0.1)   #pconei
E=(1) #(1 1 1 1 1 1) #wee_mult
F=(1) #wii_mult
G=(1) #(0.1 0.4 0.7 1.0 1.4 1.8) #wie_mult
H=(1) #wei_mult
I=(0.8) #iE
J=(0.01) # iI
K=(1.0)  # gks_
L=(0.0) # stim_val
M=(1) #(1 2 3 4 5) #run
N=(1) #( 2 3 5 7) #w_mult
O=(1) #(1.8 1.4 1.0 0.7 0.4 0.1) # wee_nmda_mult
Q=(1) # wei_nmda_mult
R=(-20.0) # log_norm_mu
S=(10.0) # log_norm_std
PP=(0) #PP=(0 1 2 3 4 5)


A=(0.1) #pconee
B=(0.1) #pconii
C=(0.1) #pconie
D=(0.1)   #pconeii
E=(1.2 1.2 1.2 1.2) #(1 1 1 1 1 1) #wee_mult
F=(1) #wii_mult
G=(0.8 1.0 1.2 1.4) #(0.1 0.4 0.7 1.0 1.4 1.8) #wie_mult
H=(1) #wei_mult
I=(0.0) #iE
J=(0.0) # iI
K=(0.4 0.6 0.8 1.0 1.2)  # gks_
L=(0.0) # stim_val
M=(1) #(1 2 3 4 5) #run
N=(1) #( 2 3 5 7) #w_mult
O=(1.5 0.9 0.46 0.37) #(1.8 1.4 1.0 0.7 0.4 0.1) # wee_nmda_mult
Q=(1) # wei_nmda_mult
R=(-20.0) # log_norm_mu
S=(9.4) # log_norm_std
T=(0.1) #frequency
U=(4.0) #strength
V=(2.0) # duration
W=(3) # net
X=(20000) # run time
PP=(0 1 2 3) #PP=(0 1 2 3 4 5)




#.1 .1 .1 .1 1.2 1 .8 1 0.0 0.0 1.2 0.0 1 1 1.5 -20.0 9.4  .1 4.0 2.0 3 2111


echo Pconee
echo ${A}
echo Pconii
echo ${B}
echo Pconie
echo ${C}
echo Pconei
echo ${D}
echo WEE
echo ${E}
echo WII
echo ${F}
echo WIE
echo ${G}
echo WEI
echo ${H}
echo IE
echo $I
echo II
echo $J
echo GKS
echo $K
echo STIM
echo $L
echo NMDA
echo $O

module load  rclone
cwd=$(pwd)


ite=0
daydate=$(date '+%d')
datepeice[0]=$(date '+%m') #month
datepeice[1]=$(date '+%d') #day
datepeice[2]=$(date '+%y') #year

timestart[0]=$(date '+%d')
timestart[1]=$(date '+%H')
timestart[2]=$(date '+%M')
timestart[3]=$(date '+%S')


rundatadir=${cwd}"/rundata"
mkdir "$rundatadir"

rundir=${cwd}"/rundata/Rich_Neuro_Mod_"${datepeice[0]}"_"${datepeice[1]}"_"${datepeice[2]}
mkdir "$rundir"

for a in "${A[@]}" ; do for b in "${B[@]}"; do for c in "${C[@]}"; do for d in "${D[@]}"; do for pp in "${PP[@]}";

do for f in "${F[@]}"; do for h  in "${H[@]}"; do for i in "${I[@]}"; do for j in "${J[@]}";
do for k in "${K[@]}"; do for l in "${L[@]}"; do for m  in "${M[@]}"; do for n  in "${N[@]}"; do for q in "${Q[@]}"; 
do for rr in "${R[@]}"; do for s in "${S[@]}"; do for t in "${T[@]}"; do for u in "${U[@]}"; do for v in "${V[@]}"; do for w in "${W[@]}"; do for x in "${X[@]}"; do
            #make excecution directory
            
	  #  #newdir=${cwd}"/gks_"${k}"_pee_"${a}"_pii_"${b}"_pie_"${c}"_pei_"${d}"_wee_"${e}"_wii_"${f}"_wie_"${G[pp]}"_wei_"${h}"_CIe_"${i}"_CIi_"${j}"_stim_"${l}"_w_mult_"${n}"_nmda_strength_"${O[pp]}"_run_"${m}"_nmda_scan"

	    
           # mkdir "$newdir"
            #cd to execution directory
          #  cd "$newdir"
         #   cp "${cwd}/network_dyn.hpp"    "${newdir}/network_dyn.hpp"
           
            
            #"${cwd}/AnSyPr" $k 3.0 $l 0.0 1.0 1.0 0.02 0.02
             #"${cwd}/AnSyPr" $k 3.0 $l 1.0 1.0 1.0 0.04 0.02
             #"${cwd}/AnSyPr" 0.01 0.02 0.02 0.01 0.005 0.00190 0.001175 0.05 1.0 0.1 1.2 1.4
        #     echo "pee pii pie pei wee wii wie wei iE iI gks stim run"
       # echo  $a  $b $c $d "${E[pp]}" $f "${G[pp]}" $h $i $j $k $l $m $n "${O[pp]}"
       #echo  "${E[pp]}" "${G[pp]}" "${O[pp]}"
        

	#"${cwd}/AnSyPr3" $a  $b $c $d $e $f "${G[pp]}" $h $i $j $k $l $n "${O[pp]}"      
	#sbatch   --export="gks=$k,pee=$a,pii=$b,pie=$c,pei=$d,wee_mult=${E[pp]},wii_mult=${G[pp]},wie_mult=${G[pp]},wei_mult=${E[pp]},CIe=$i,CIi=$j,stim=$l,w_mult=$n,wee_nmda_mult=${O[pp]},run=$m,iter=$ite,wei_nmda_mult=${O[pp]},log_norm_mu=$rr,log_norm_std=$s,noise_freq=$t,noise_strength=$u,noise_dur=$v,net=$w,run_time=$x" "${cwd}/richneuromod_log_norm_run_noise_variable_net_ruin.batch"
    sbatch   --export="gks=$k,pee=$a,pii=$b,pie=$c,pei=$d,wee_mult=${E[pp]},wii_mult=${G[pp]},wie_mult=${G[pp]},wei_mult=${E[pp]},CIe=$i,CIi=$j,stim=$l,w_mult=$n,wee_nmda_mult=${O[pp]},run=$m,iter=$ite,wei_nmda_mult=${O[pp]},log_norm_mu=$rr,log_norm_std=$s,noise_freq=$t,noise_strength=$u,noise_dur=$v,net=$w,run_time=$x" "${cwd}/simulation_send.batch"
#sbatch "${cwd}/richneuromod_run.batch"    --export  "gks=$k,pee=$a,pii=$b,pie=${c},pei=$d,wee=$e,wii=$f,wie=$h,CIe=$i,CIi=$j,stim=$l,w_mult=$n,nmdaee=${O[pp]},run=${m}"   

let "ite=ite+1"

sleep 1


echo "gks pee pii pie pei wee wii  wie     wei CIe CIi stim w_mult wee_nmda_mult run iter wee_nmda_mult  log_norm_mu   log_norm_std "
echo "$k  $a  $b  $c  $d  ${E[pp]} $f   ${G[pp]} $h  $i  $j  $l   $n     ${O[pp]}     $m $ite $q             $rr            $s           "

done; done; done; done; done; done; done; done; done; done; done; done; done; done; done; done; done; done; done; done; done;
