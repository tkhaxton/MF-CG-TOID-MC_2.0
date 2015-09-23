#!/bin/bash

inputfile="../output/monolayer_search/c10000.best.config"
directory="../output/bilayer_search"
mkdir -p $directory

time ../code/bilayer_from_monolayer -inputfile $inputfile -base $directory -chemistry 1 1 24 2 1 2 7 3 7 -leafxoffsetfrac 0.8 -leafyoffsetfrac 0.0 -leafspacing 12.0 -temperature 10000 -cycles 10000



