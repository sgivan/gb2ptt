#!/bin/bash

toolpath=$( echo $0 | sed 's/\/setup.sh//' )

echo "splitting $1"
ds=$( $toolpath/seqsplit.pl -f $1 -F genbank -O genbank -n 1 )
string=''

echo 'creating directory for each input file'
cnt=0
for dir in *.gbk
do
    echo $dir
    newdir=$( echo $dir | sed 's/.gbk//' )
    echo $newdir
    mkdir -p $newdir
    mv -f $dir $newdir
    cd $newdir
    echo "gb2ppt.pl --infile $dir"
    $( $toolpath/gb2ptt.pl --infile $dir )
    echo "fasta_format.pl -f $dir -i genbank -O $dir.fna"
    $( $toolpath/fasta_format.pl -f $dir -i genbank -O $dir.fna )
    curdir=$( pwd )
    ((cnt++))
    if [[ $cnt == 1 ]]; then
        string=$curdir
    else
        string="$string,$curdir"
    fi
    cd ..
done

echo $string
