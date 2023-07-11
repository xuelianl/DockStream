#!/bin/bash

dockfiles_indock_dir=$1
# indock=$2
# db2s=$3
dock_path=$2


pipeline(){
    cp -r $dockfiles_indock_dir/dockfiles .
    cd docking
    cp $dockfiles_indock_dir/INDOCK .
    # readlink -f $db2s > split_database_index
    $dock_path/docking/DOCK/bin/dock64 INDOCK

    # conda info --envs
    # source /pubhome/xli02/opt/miniconda/etc/profile.d/conda.sh
    # conda deactivate
    # conda info --envs
    cd ..
    ls docking -d > dirlist
    python2 $dock_path/analysis/extract_all_blazing_fast.py ./dirlist extract_all.txt 1000.0 # 1000 is maxenergy
    # modify getposes_blazing_faster.py: add ':' at line 36
    python2 $dock_path/analysis/getposes_blazing_faster.py '' extract_all.sort.uniq.txt 5000 poses.mol2 test.mol2.gz
}

pipeline &> dock.log
