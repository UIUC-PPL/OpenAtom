#!/bin/bash
for beta in `seq 2.0 0.5 4.5`
	do 
	rm betaconv$beta.dat
	for i in `seq 4 2 14`
 		do ./model_PAW.x PAW.in $i 10 20 0 $beta > tmp.txt
			EHarselfscrA=`grep "E_Har_self_0D_scr" tmp.txt|head -1|awk '{print $3}'`
			EHarselfscrG=`grep "E_Har_self_0D_scr" tmp.txt|head -1|awk '{print $5}'`
			E0=`grep "E_Har_self_0D" tmp.txt|head -1|awk '{print$3}'`
			diffA=`echo '('$EHarselfscrA-$EHarselfscrG')'/$EHarselfscrG|bc -l`
			diffgridA=`echo '('$EHarselfscrG-$E0')'/$E0*100|bc -l`
			rm tmp.txt
			echo $i   $diffA >> betaconv$beta.dat
	done
done

