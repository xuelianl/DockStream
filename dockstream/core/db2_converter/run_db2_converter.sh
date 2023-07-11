#!/bin/bash

line=$1
max_conf=$2
sampletp=$3
checkstereo=$4
useff=$5

name=`echo $line|awk '{print $2}'|awk -F ':' '{print $1}'`
echo $line > $name.smi
/pubhome/xli02/git/db2_converter/ligand/generate/build_ligand.sh $name $max_conf $sampletp $checkstereo $useff &> conformator_$name.log
