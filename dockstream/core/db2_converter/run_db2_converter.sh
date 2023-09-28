#!/bin/bash

line=$1
max_conf=$2
workingpath=$3
method=$4
checkstereo=$5
useff=$6
sampletp=$7

source /pubhome/xli02/opt/miniconda/etc/profile.d/conda.sh
conda activate basic
name=`echo $line|awk '{print $2}'|awk -F ':' '{print $1}'`
echo $line > $name.smi
# /pubhome/xli02/scripts/db2_converter/ligand/generate/build_ligand.sh $name $max_conf $sampletp $checkstereo $useff &> conformator_$name.log
command="build_ligand -i $name.smi -m $max_conf --workingpath $workingpath --outputpath $workingpath --method $method"
if [ $checkstereo == 'True' ];
then
    command+=" -c"
fi

if [ $useff == 'True' ];
then
    command+=" -f"
fi

if [ $sampletp == 'True' ];
then
    command+=" -s"
fi

# echo $command
$command
conda deactivate &> db2_${method}_$name.log
