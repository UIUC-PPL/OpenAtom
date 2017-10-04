#!/bin/bash

echo "Making STATES directory..."
mkdir STATES

numk=$1
numspin=$2

numkc=`expr $numk - 1`
numspinc=`expr $numspin - 1`

echo "number of k points $numk"

cd STATES

for s in `seq 0 $numspinc`
do
  for i in `seq 0 $numkc`
  do
    mkdir Spin.${s}_Kpt.${i}_Bead.0_Temper.0
  done
done  

cd ..

echo "making directories for shifted wavefunctions"

cd STATES

for s in `seq 0 $numspinc`
do
  for i in `seq 0 $numkc`
  do
    mkdir Spin.${s}_Kpt.0${i}_Bead.0_Temper.0
  done
done

