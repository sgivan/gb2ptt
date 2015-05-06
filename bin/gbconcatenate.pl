#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  gbconcatenate.pl
#
#        USAGE:  ./gbconcatenate.pl  
#
#  DESCRIPTION:  Script to concatenate multiple GenBank sequences into a single
#                   GenBank sequence file.
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  05/06/15 09:52:31
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use warnings;
use strict;
use Bio::SeqIO;

my ($debug,$verbose,$help,$infile,$joinseq);

my $result = GetOptions(
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
    "infile:s"  =>  \$infile,
    "joinseq:s" =>  \$joinseq,
);

if ($help) {
    help();
    exit(0);
}

sub help {

    say <<HELP;

--debug
--verbose
--help
--infile

HELP

}

$joinseq = 'NNNNNNNNNNNNTCCANNNNNNNNNNNN';
my $joinlength = length($joinseq);

my $seqio = Bio::SeqIO->new(
    -format     =>  'genbank',
    -file       =>  $infile,
    -verbose    =>  -1,
);

my $outseq = Bio::SeqIO->new(
    -format     =>  'genbank',
    -file       =>  '>concatseq.gbk',
);

my $concatseq;
my $cnt = 0;
while (my $seq = $seqio->next_seq()) {
    ++$cnt;

    if ($cnt == 1) {
        $concatseq = $seq;
        next;
    }

    my $inlength = $concatseq->length();
    my @seqfeatures = $seq->get_SeqFeatures();
    if ($debug) {
        say "inlength = '$inlength'";
        say $seq->desc() . " has " . scalar(@seqfeatures) . " features" if ($debug);
    }

    $concatseq->seq($concatseq->seq() . $joinseq . $seq->seq());
    for my $seqfeature (@seqfeatures) {
        $seqfeature->start($seqfeature->start() + $inlength + $joinlength);
        $seqfeature->end($seqfeature->end() + $inlength + $joinlength);
        $concatseq->add_SeqFeature($seqfeature);
    }
}


$outseq->write_seq($concatseq);

