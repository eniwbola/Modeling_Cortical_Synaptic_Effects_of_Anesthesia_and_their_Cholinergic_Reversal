# Modeling Cortical Synaptic Effects of Anesthesia and their Cholinergic Reversal

Modeling and Analysis code for study "Modeling Cortical Synaptic Effects of Anesthesia and their Cholinergic Reversal

## anesthetic_sim.cpp
cpp file for simulating model/all simulation spike trains.

## simulation_send.batch
batch file that send file to clusther, this can rwmoved if doing this locatlly

## simulation_control_bash.sh
This is used to simulated specific paramters for different intial condictions  and network structure


## gen_alg.m
  This file evaluates the new simualations by comparing them to the experimental data and then generatesa new set of paramters

## gen_alg_bash_control.sh

Shell file that controls loop for genetic algorrithm. The file first calls the matlab gen_alg file which evaluates the current generation of simulations and removes the highest cost simulations as well as generaties a new candidate set of paramters. After this is done the shell file then calls the simulation send abtch file to simulate those paramters.

## BaseTS 
Raw experimental spike data

## *_dat_low files
statistical Analysis of simulation outputs




