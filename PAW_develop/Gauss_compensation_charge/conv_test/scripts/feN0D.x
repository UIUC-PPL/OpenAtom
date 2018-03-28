#!/bin/bash
rm feN.dat
#for i in `seq 2 2 20`
# do i2=`echo $i*2-2|bc -l`
#	./model_PAW.x PAW.in 14 $i $i2 > tmp.txt
for i in `seq 2 1 15`
 do ./model_PAW_feN.x PAW.in $i 10 20 0 2.0 > tmp.txt
# for i in `seq 0 1 10`
#  do ./model_PAW.x PAW.in 14 10 20 $i > tmp.txt
	f0D=`grep "fx0g \[0" tmp.txt|head -1|awk '{print $5}'`
	f0DA=`grep "fx0  \[0" tmp.txt|head -1|awk '{print $5}'`
	rm tmp.txt
	echo $i  $f0D $f0DA >> feN.dat
done
