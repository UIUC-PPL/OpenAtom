#!/bin/bash
for nr in 20 30
    do  rm testntchknr$nr.dat
		for nt in 2 4 6 8 10 12 14 16 
			do 	nphi=`echo 2*$nt|bc -l`
				./model_PAW.x PAW.in $nr $nt $nphi  0 3.0 > tmp.txt
				answer=`grep result4 tmp.txt|head -1|awk '{print $2}'`
				Hself=`grep E_Har_self_0D_scr tmp.txt|awk '{print $5}'`
				eN0D=`grep E_eN_0D tmp.txt|head -1|awk '{print $5}'`
				eN0Dself=`grep E_eNself_0D tmp.txt|head -1|awk '{print $5}'`
				Elong3D=`grep E_long_3D tmp.txt|head -1|awk '{print $5}'`
				fxg=`grep -e "fxg  " tmp.txt|head -1|awk '{print $5}'`
#				H0D=`grep E_Har_0D tmp.txt|head -1|awk '{print $5}'`
#				H0Dself=`grep E_Har_self_0D tmp.txt|head -1|awk '{print $5}'`
				rm tmp.txt
				echo  $nt $eN0D $eN0Dself $Elong3D $fxg>> testntchknr$nr.dat
		done
done
