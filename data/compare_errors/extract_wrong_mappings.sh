#!/bin/bash
gaf=$1
../k8-Linux ../paftools.js mapeval -Q0 $gaf | sed -n "/^E/p" -  > wrong_mappings_$gaf
../k8-Linux ../paftools.js mapeval $gaf | awk '{ print $6 }' | tail -n 1 - > count_mappings_$gaf
