#!/bin/bash
#for beta in 3.0
for beta in 2.0 3.0 4.0
do	rm analytic_res$beta.dat
    for n in `seq 10 10 1000`
	do  ./testortho.x 1 1 $n $beta > tmp.txt
		result=`grep result4 tmp.txt|awk '{print $2}'`
		result3=`grep result3 tmp.txt|awk '{print $2}'`
		result5=`grep result5 tmp.txt|awk '{print $2}'`
		result6=`grep result6 tmp.txt|awk '{print $2}'`
		rm tmp.txt
		echo $beta   $n $result $result3 $result5 $result6 >> analytic_res$beta.dat
	done
done
