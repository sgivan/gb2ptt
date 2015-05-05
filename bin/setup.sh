#!/bin/bash

echo 'splitting input file'
ds=$(~/projects/gb2ptt/bin/seqsplit.pl -f gbkfile -F genbank -O genbank -n 1)
string=''

echo 'creating directory for each input file'
cnt=0
for dir in *.gbk
do
    echo $dir
    newdir=$(echo $dir | sed 's/.gbk//')
    echo $newdir
    mkdir -p $newdir
    mv -f $dir $newdir
    cd $newdir
    echo "gb2ppt.pl --infile $dir"
    $(~/projects/gb2ptt/bin/gb2ptt.pl --infile $dir)
    echo "fasta_format.pl -f $dir -i genbank -O $dir.fna"
    $( ~/bin/fasta_format.pl -f $dir -i genbank -O $dir.fna)
    curdir=$(pwd)
    ((cnt++))
    if [[ $cnt == 1 ]]; then
        string=$curdir
    else
        string="$string,$curdir"
    fi
    cd ..
done

echo $string
