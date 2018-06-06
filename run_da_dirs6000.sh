#!/bin/bash

nr=0
for seed in `cat seed4.txt`
do
		
		nr=$(($nr+1))
		./brownian res300k_$nr statis300k_$nr $seed 300000 >out_300k_$nr.txt &
		
		nr_runnings=`ps -e | grep "brownian" | wc -l`
		while [ $nr_runnings -ge 30 ]
		do
		  sleep 3
		  nr_runnings=`ps -e | grep "brownian" | wc -l`
		done
done
