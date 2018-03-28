#!/bin/bash
for nr in `seq 4 2 14`
#for nr in 14
	do rm diff.dat
#	rm EHarselfscr$nr.dat
	for i in `seq 1.0 0.2 6.2`
 		do ./model_PAW.x PAW.in $nr 10 20 0 $i > tmp.txt
			EHarselfscrA=`grep "E_Har_self_0D_scr" tmp.txt|head -1|awk '{print $3}'`
			EHarselfscrG=`grep "E_Har_self_0D_scr" tmp.txt|head -1|awk '{print $5}'`
			E0=`grep "E_Har_self_0D" tmp.txt|head -1|awk '{print$3}'`
			diffA=`echo '('$EHarselfscrA-$E0')'/$E0*100|bc -l`
			diffgridA=`echo '('$EHarselfscrG-$E0')'/$E0*100|bc -l`
			rm tmp.txt
			echo $i   $diffA $diffgridA >> diff.dat
	done
done

