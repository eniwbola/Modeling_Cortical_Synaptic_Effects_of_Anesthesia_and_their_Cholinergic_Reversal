# Modeling Cortical Synaptic Effects of Anesthesia and their Cholinergic Reversal

Modeling and Analysis code for study "Modeling Cortical Synaptic Effects of Anesthesia and their Cholinergic Reversal"

## anesthetic_sim.cpp
cpp file for simulating model/all simulation spike trains.

## simulation_send.batch
batch file that sends file to cluster, this can removed if doing this locally

## simulation_control_bash.sh
This is used to simulated specific paramters for different intial condictions  and network structure


## gen_alg.m
This file evaluates the new simualations by comparing them to the experimental data and then generatesa new set of paramters

## gen_alg_bash_control.sh

Shell file that controls loop for genetic algorrithm. The file first calls the matlab gen_alg file which evaluates the current generation of simulations and removes the highest cost simulations as well as generaties a new candidate set of paramters. After this is done the shell file then calls the simulation send abtch file to simulate those paramters.

## BaseTS 
Raw experimental spike data. Data collected in study leading to prior publication:
Vizuete JA, Pillay S, Diba K, Ropella KM, Hudetz AG. Monosynaptic functional connectivity in cerebral cortex during wakefulness and under graded levels of anesthesia. Front Integr Neurosci. 2012 Oct 12;6:90. doi: 10.3389/fnint.2012.00090. PMID: 23091451; PMCID: PMC3469825.


## *_dat_low files
statistical Analysis of simulation outputs




