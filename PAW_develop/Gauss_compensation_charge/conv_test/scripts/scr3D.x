#!/bin/bash
rm diff.dat
#for nr in `seq 4 2 14`
for nr in 14
#	rm EHarselfscr$nr.dat
  do   for i in `seq 1.0 0.5 6.0`
 		do 	./model_PAW.x PAW.in $nr 10 20 0 $i > tmp.txt
			EHarselfscrA=`grep "E_Har_short_self_3D_scr" tmp.txt|head -1|awk '{print $3}'`
			EHarselfscrG=`grep "E_Har_short_self_3D_scr" tmp.txt|head -1|awk '{print $5}'`
			E0=`grep "E_Har_short_self_3D" tmp.txt|head -1|awk '{print$3}'`
			diffA=`echo '('$EHarselfscrA-$E0')'/$E0*100|bc -l`
			diffgridA=`echo '('$EHarselfscrG-$E0')'/$E0*100|bc -l`
			rm tmp.txt
			echo $i   $diffA $diffgridA >> diff.dat
	done
done

