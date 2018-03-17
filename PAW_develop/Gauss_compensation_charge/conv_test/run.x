#!/bin/bash
rm res.dat
echo "rorder E3D   E0D  E3D_A  E0D_A   f3D   f0D  f3D_A   f0D_A" > res.dat
#for i in `seq 2 2 20`
# do i2=`echo $i*2-2|bc -l`
#	./model_PAW.x PAW.in 14 $i $i2 > tmp.txt
for i in `seq 2 1 15`
 do ./model_PAW.x PAW.in $i 10 20 3 > tmp.txt
# for i in `seq 0 1 10`
#  do ./model_PAW.x PAW.in 14 10 20 $i > tmp.txt
	E3D=`grep "E_tot_3D" tmp.txt|awk '{print $5}'`
	E0D=`grep "E_tot_0D" tmp.txt|awk '{print $5}'`
	E3DA=`grep "E_tot_3D" tmp.txt|awk '{print $3}'`
	E0DA=`grep "E_tot_0D" tmp.txt|awk '{print $3}'`
	f3D=`grep "fxg  \[0" tmp.txt|head -1|awk '{print $5}'`
	f0D=`grep "fx0g \[0" tmp.txt|head -1|awk '{print $5}'`
	f3DA=`grep "fx   \[0" tmp.txt|head -1|awk '{print $5}'`
	f0DA=`grep "fx0  \[0" tmp.txt|head -1|awk '{print $5}'`
	rm tmp.txt
	echo $i   $E3D   $E0D  $E3DA $E0DA  $f3D   $f0D $f3DA $f0DA >> res.dat
done

