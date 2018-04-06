#!/bin/bash
rm *.dat
#for beta in 2.0
for beta in 2.0 3.0 4.0
do  for n in `seq 10 10 100`
	do  ./testortho.x 1 1 $n $beta > tmp.txt
		result=`grep result4 tmp.txt|awk '{print $2}'`
		result3=`grep result3 tmp.txt|awk '{print $2}'`
		rm tmp.txt
		echo $beta   $n $result $result3 >> res$beta.dat
	done
done
