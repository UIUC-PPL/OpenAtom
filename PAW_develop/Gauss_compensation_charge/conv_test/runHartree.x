#!/bin/bash
rm res.dat
echo "rorder E3D   E0D  E3D_A  E0D_A   f3D   f0D  f3D_A   f0D_A" > res.dat
#for i in `seq 2 2 20`
# do i2=`echo $i*2-2|bc -l`
#	./model_PAW.x PAW.in 14 $i $i2 > tmp.txt
for i in `seq 2 2 14`
 do ./model_PAW.x PAW.in $i 10 20 3 > tmp.txt
# for i in `seq 0 1 10`
#  do ./model_PAW.x PAW.in 14 10 20 $i > tmp.txt
	Har0D=`grep "E_Har_0D" tmp.txt|awk '{print $5}'`
	Har0DA=`grep "E_Har_0D" tmp.txt|awk '{print $3}'`
	Harself0D=`grep "E_Har_self_0D" tmp.txt|awk '{print $5}'`
	Harself0DA=`grep "E_Har_self_0D" tmp.txt|awk '{print $3}'`
	rm tmp.txt
	echo $i $Har0D $Har0DA $Harself0D $Harself0DA  >> res.dat
done

