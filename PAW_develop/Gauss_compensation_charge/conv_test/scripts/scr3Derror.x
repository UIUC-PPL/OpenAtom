#!/bin/bash
rm diff.dat
#for nr in `seq 4 2 14`
for nr in 14
	do	rm EHarshortselfscr$nr.dat
		for i in `seq 1.0 0.5 6.0`
 		do 	./model_PAW.x PAW.in $nr 10 20 0 $i > tmp.txt
			EHarselfscrA=`grep "E_Har_short_self_3D_scr" tmp.txt|head -1|awk '{print $3}'`
			EHarselfscrG=`grep "E_Har_short_self_3D_scr" tmp.txt|head -1|awk '{print $5}'`
			rm tmp.txt
			echo $i   $EHarselfscrA $EHarselfscrA >> EHarshortselfscr$nr.dat
	done
done

