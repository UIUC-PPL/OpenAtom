#!/bin/bash
rm res.dat
echo "rorder E3D E0D f3D f0D" > res.dat
for i in `seq 2 2 20`
 do i2=`echo $i*2|bc -l`
	./model_PAW.x PAW.in 14 $i $i2 > tmp.txt
	E3D=`grep "E_tot_3D" tmp.txt|awk '{print $5}'`
	E0D=`grep "E_tot_0D" tmp.txt|awk '{print $5}'`
	f3D=`grep "fxg" tmp.txt|head -1|awk '{print $5}'`
	f0D=`grep "fx0g" tmp.txt|head -1|awk '{print $5}'`
	rm tmp.txt
	echo $i   $E3D   $E0D   $f3D   $f0D >> res.dat
done

