#!/bin/bash

inputfile="../output/bilayer_stack/coord"
directory="../output/bilayer_stack_evaporate"
inputfile2="../output/bilayer_stack_evaporate/coord"
directory2="../output/bilayer_stack_evaporate_relax"
mkdir -p $directory
mkdir -p $directory2

time ../code/peptoid -inputfile $inputfile -base $directory -ictype 0 -cycles 100000 -wholebilayertranslatefreq 0.1 -shiftbilayergapfreq 0.1 -reset 1 -normalforceperunitarea 0.001

time ../code/peptoid -inputfile $inputfile2 -base $directory2 -ictype 0 -cycles 500000 -wholebilayertranslatefreq 0.1 -shiftbilayergapfreq 0.1 -reset 1



