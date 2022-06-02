#!/bin/sh


#get current directory 


#echo $1 
#echo $2 


# bedit 212 if b (directory size)  is equal to base 20 and baseit is zero , run new paramters from last matlab calculation 
	     # so this should really say that if b<= 20 then generate 30 -b from matlab , then send to flux 
		 # so I need to figure out how to stop this so check the flux dir3ectory 
		 #  and wait until its done  
		 # so after checking the flux and its done, revaluate and do the same thing and keep going 

#so so first get rid of the 10 done file and just ask if there is nothing in the flux directory 
# bedit we should  change this to above 30 because matlab will get rid of the ones greater than 20 and create 10 and discard the rest, make sure this is the case

# so think there are two stages, generating the paramters in matlab and sending them out 
# so to solve this 
#


# The process is send out the 10 to  make it go back to 30, then run matlab to remove the 10 and gneration paramters 
# so I need to check what is running in flux, if nothing then check how many are in the directory, 
# if less than 30  call matlab make 30- x based on the prvious 
# and then call the other function to run the last 30 -x 
#  do this until it is 30  
# in each scenerero put the parameters in the in the file and when there are 30 write these paramters to to the updating paramter file
# so now I need to take out the current write to the up dating parameter file 

# then change maltab 
# then add the param chanage 
 
#so all i need to do now is say , -> if less than 30 matlab generates 30-x and gen sends them out 
# So only take from the top highest  to generate the 30 -x , and generate but if not 
# only if it actually is 30 should matlab  send to rejects to make 20 
	# so steps make sure that the call to run is  


# so how do i do this? so 
# first make if and say if 30 then then send to rejects 


cwd=$(pwd)


base_param_nums=20  #20  # these are the paramters that you hope that are there after matlabshavs off
total_param_nums=30			  
my_time=5000  

ite=0
base_ite=0
ze=0
ons=1

module load matlab
while(($ite<5))
do
	base_ite=0
	c_ite=0
	let "ite=ite+1"
	echo you are on iter
	echo $ite
	Pl_ite=0
	check_file=0		


	running_tot=$(squeue -u eniwbola -h -t running -r | wc -l)   
	pending_tot=$(squeue -u eniwbola -h -t pending -r | wc -l)  
	d=$(squeue -u eniwbola -h -t pending -r | wc -l)                                                                                                                                                 echo d
	pend_run_tot=$((running_tot+pending_tot))  #$((d+e))

	while true	
	do	


 
		running_tot=$(squeue -u eniwbola -h -t running -r | wc -l)   
		pending_tot=$(squeue -u eniwbola -h -t pending -r | wc -l)  
		d=$(squeue -u eniwbola -h -t pending -r | wc -l)                                                                                                                                                 echo d
		pend_run_tot=$((running_tot+pending_tot))  #$((d+e))
		#10end_file.txt

		#rm -r *__* 
		#rm "-r *__*"
		dir=$(pwd)


		echo dir
		echo $dir
		b=$(ls -ld gks* | grep "^d" | wc -l)
		echo gks_grep
		echo $b
		echo end_gks_grep	

		echo first_set
		total_param_nums=30	
		echo pend_run_tot
		echo "$pend_run_tot"
		echo running_tot
		echo "$running_tot"

		echo total_param_nums
		echo $total_param_nums

		echo gks_grep
		echo $b


		echo base_param_nums
		echo $base_param_nums
		if [ "$b" -ge "$total_param_nums" ] && [ "$pend_run_tot" -eq "$ze" ]; 
		then

			module load matlab

			dirs="'$dir'"

			echo running matlab

			#matlab -nodisplay -r 	 "classifcation_strength_neurons_fix_fig_gen_alg_5_5; exit"
			matlab -nodisplay -r 	 "classifcation_strength_neurons_fix_fig_gen_al5; exit"
			#matlab -nodisplay -r 	 "classifcation_strength_analysis_subset_neurons_fix_fig_gen_al5; exit"

			sleep 4
			break	
		

# this is going to be the case when you have less than 30 but not 20

		elif [ \( "$b" -lt "$total_param_nums" -a "$pend_run_tot" -eq "$ze" \) ]; then


			#echo second_set
			total_param_nums=30	
			echo pend_run_tot
			echo "$pend_run_tot"
			echo running_tot
			echo "$running_tot"

			#echo total_param_nums
			#echo $total_param_nums

			echo b
			echo $b

			#echo ze
			#echo $ze

			echo base_param_nums
			echo $base_param_nums

			if [ \( "$b" -gt "$base_param_nums" -a "$pend_run_tot" -eq "$ze"  -a  "$b" -le "$total_param_nums" \) ]; then
				# at this point is its greater than 20 but less than 30 run matlab to shave off to get to 20
				matlab -nodisplay -r 	 "gen_alg; exit"#"classifcation_strength_neurons_fix_fig_gen_al5; exit"
				#matlab -nodisplay -r 	 "classifcation_strength_analysis_subset_neurons_fix_fig_gen_al5; exit"

				#run_FILE="new_param.txt"  #"remaining_param.txt"

			elif [ \( "$b" -eq "$base_param_nums" -a "$pend_run_tot" -eq "$ze" \) ]; then # bedit 212 if b (directory size)  is equal to base 20 and baseit is zero , run new paramters from last matlab calculation
				# at this point if directors is equat to 20
                                                        # so this should really say that if b<= 20 then generate 30 -b from matlab , then send to flux
																	 # so I need to figure out how to stop this so check the flux dir3ectory
																	 #  and wait until its done  
																	 # so after checking the flux and its done, revaluate and do the same thing and keep goin
				if [ -s new_param.txt ];
				then
				
				#matlab -nodisplay -r 	 "classifcation_strength_analysis_subset_neurons_fix_fig_gen_al5; exit"
				matlab -nodisplay -r 	 "gen_alg; exit" #"classifcation_strength_neurons_fix_fig_gen_al5; exit"
				#run_FILE="remaining_param.txt"
				else
				#matlab -nodisplay -r 	 "classifcation_strength_analysis_subset_neurons_fix_fig_gen_al5; exit"
				matlab -nodisplay -r 	 "gen_alg; exit" #"classifcation_strength_neurons_fix_fig_gen_al5; exit"
				#run_FILE="new_param.txt"
				fi


			fi	
	
	
			echo running_new_params
			# so this an
			pee_new=()

			while IFS=", " read -r -a a
			do
		
					anew=${a[0]}

		
					t=${anew[0]}

					pee_new+=($(echo $t |cut -d " " -f 1))

				
	

				A+=($(echo $t |cut -d " " -f 1))
				B+=($(echo $t |cut -d " " -f 2))
				C+=($(echo $t |cut -d " " -f 3))
				D+=($(echo $t |cut -d " " -f 4))
				E+=($(echo $t |cut -d " " -f 5))
				F+=($(echo $t |cut -d " " -f 6))
				G+=($(echo $t |cut -d " " -f 7))
				H+=($(echo $t |cut -d " " -f 8))
				I+=($(echo $t |cut -d " " -f 9))
				J+=($(echo $t |cut -d " " -f 10))
				K+=($(echo $t |cut -d " " -f 11))
				L+=($(echo $t |cut -d " " -f 12))
				M+=($(echo $t |cut -d " " -f 13))
				N+=($(echo $t |cut -d " " -f 14))
				O+=($(echo $t |cut -d " " -f 15))
				Q+=($(echo $t |cut -d " " -f 16))
				RN=($(echo $t |cut -d " " -f 17))
				R+=($(expr $RN*$neg | bc))
				S+=($(echo $t |cut -d " " -f 18))
				T+=($(echo $t |cut -d " " -f 19))
				U+=($(echo $t |cut -d " " -f 20))
				V+=($(echo $t |cut -d " " -f 21))
				W+=($(echo $t |cut -d " " -f 22))
				X+=($(echo $t |cut -d " " -f 23))



				let "base_ite=base_ite+1"
		
				let "c_ite=c_ite+1"

			done < "new_param.txt" #run_FILE


			#echo ${pee_new[@]}
			let "base_ite=base_ite+1"
			echo c_ite
			echo $c_ite
			echo base_ite
			echo $base_ite
			echo ite
			echo $ite
			PL=(1 2 3 4 5 6 7 8 9 10)
			#for pl in "${PL[@]}"; do
			#for pl in {1..$base_ite}; do
			for (( c=1; c<=$c_ite; c++ )); do

				R_fix=(-20.0)
				echo "gks pee pii pie pei wee_mult wii_mult wie_mult wei_mult CIe Cii stim w_mult wee_nmda_mult run iter wei_nmda_mult log_norm_mu log_nor_std noise_freq noise_strength noise__dur net run_time"
				echo "${K[pl]}" "${A[pl]}" "${B[pl]}" "${C[pl]}"  "${D[pl]}"  "${E[pl]}" "${G[pl]}" "${G[pl]}"  "${E[pl]}" "${I[pl]}" "${J[pl]}" "${L[pl]}" "${N[pl]}" "${O[pl]}" "${M[pl]}" "$Pl_ite" "${O[pl]}""$R_fix" "${S[pl]}" "${T[pl]}" "${U[pl]}" "${V[pl]}" "${W[pl]}" "$my_time" 

				let "Pl_ite=Pl_ite+1"
				sbatch   --export="gks=${K[pl]},pee=${A[pl]},pii=${B[pl]},pie=${C[pl]},pei=${D[pl]},wee_mult=${E[pl]},wii_mult=${G[pl]},wie_mult=${G[pl]},wei_mult=${E[pl]},CIe=${I[pl]},CIi=${J[pl]},stim=${L[pl]},w_mult=${N[pl]},wee_nmda_mult=${O[pl]},run=${M[pl]},iter=$Pl_ite,wei_nmda_mult=${O[pl]},log_norm_mu=$R_fix,log_norm_std=${S[pl]},noise_freq=${T[pl]},noise_strength=${U[pl]},noise_dur=${V[pl]},net=${W[pl]},run_time=$my_time" "${cwd}/richneuromod_log_norm_run_noise_variable_net_ruin_1st_twen.batch"

				sleep 1
				echo pl
				echo $pl


			done

			c_ite=0

	# send to bash

		else 
			echo hey im waiting to run matlab
			sleep 11

		fi
	done		


done




sleep 1
