#!/bin/bash
rm res.dat
#for i in `seq 2 2 20`
# do i2=`echo $i*2-2|bc -l`
#	./model_PAW.x PAW.in 14 $i $i2 > tmp.txt
for i in `seq 2 1 15`
 do ./model_PAW.x PAW.in $i 10 20 3 > tmp.txt
# for i in `seq 0 1 10`
#  do ./model_PAW.x PAW.in 14 10 20 $i > tmp.txt
	EeN0D=`grep "E_eN_0D" tmp.txt|awk '{print $5}'`
	EeN0DA=`grep "E_eN_0D" tmp.txt|awk '{print $3}'`
	EHar0DA=`grep "E_Har_0D" tmp.txt|awk '{print $3}'`
	EHar0D=`grep "E_Har_0D" tmp.txt|awk '{print $5}'`
	EHarself0DA=`grep "E_Har_self_0D" tmp.txt|awk '{print $3}'`
	EHarself0D=`grep "E_Har_self_0D" tmp.txt|awk '{print $5}'`
	rm tmp.txt
	echo $i   $EeN0D $EeN0DA $EHar0D $EHar0DA $EHarself0DA $EHarself0D >> res.dat
done
