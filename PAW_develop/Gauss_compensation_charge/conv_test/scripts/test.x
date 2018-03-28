#!/bin/bash
rm test.dat
#for i in `seq 2 2 20`
# do i2=`echo $i*2-2|bc -l`
#	./model_PAW.x PAW.in 14 $i $i2 > tmp.txt
for i in `seq 2 1 15`
 do ./model_PAW.x PAW.in $i 10 20 3 > tmp.txt
# for i in `seq 0 1 10`
#  do ./model_PAW.x PAW.in 14 10 20 $i > tmp.txt
	test=`grep "test: " tmp.txt|head -1|awk '{print $2}'`
	testA=`grep "test: " tmp.txt|head -1|awk '{print $3}'`
	rm tmp.txt
	echo $i   $test   $testA >> test.dat
done

