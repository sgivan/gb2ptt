#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  gb2ptt.pl
#
#        USAGE:  ./gb2ptt.pl  
#
#  DESCRIPTION:  Script to convert a genbank flat file to ptt
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  05/04/15 15:23:32
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use warnings;
use strict;
use Bio::SeqIO;

my ($debug,$verbose,$help,$infile);

my $result = GetOptions(
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
    "infile:s"  =>  \$infile,
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

$infile = 'test' unless ($infile);

my $seqio = Bio::SeqIO->new(
                                -file   =>  $infile,
                                -format =>  'genbank',
                            );

open(my $PTT, ">", $infile . ".ptt");
open(my $RNT, ">", $infile . ".rnt");

my %tags = ();
while (my $seq = $seqio->next_seq()) {
    my $header1 = $seq->desc() || 'unknown' . " - 1.." . $seq->length();

    if ($debug) {
        say "\$seq isa '" . ref($seq) . "'" if ($debug);
        say $seq->desc() . " - 1.." . $seq->length();
        say "CDS: " . scalar($seq->get_SeqFeatures('CDS'));
        exit();
    }

    my @CDS = $seq->get_SeqFeatures('CDS');
    my @RNA = $seq->get_SeqFeatures('tRNA');
    push(@RNA,$seq->get_SeqFeatures('rRNA'));

    my $featcnt = 0;
    
    say $PTT $header1;
    say $PTT scalar(@CDS) . " proteins";
    say $RNT $header1;
    say $RNT scalar(@RNA) . " RNAs";
    say $PTT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct";
    say $RNT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct";

    for my $feature (@CDS) {

        next if $feature->has_tag('pseudo');

        my $tag = $feature->primary_tag();
        ++$tags{$tag};


        last if ($debug && ++$featcnt > 10);

        my $start = $feature->start();
        my $stop = $feature->end();
        my $strand = $feature->strand();
        my $length = $feature->length();
        my @pid = $feature->has_tag('protein_id') ? $feature->get_tag_values('protein_id') : $feature->get_tag_values('locus_tag');
        my @gene = $feature->has_tag('gene') ? $feature->get_tag_values('gene') : '-';
        my @synonym = $feature->get_tag_values('locus_tag');
        my $code = '-';
        my $cog = '-';
        my @description = $feature->get_tag_values('product');

        if ($strand > 0) {
            $strand = '+';
        } elsif ($strand < 0) {
            $strand = '-';
        }

#            my @misc = $feature->primary_tag('misc_feature');
#
#            for my $misc_feature (@misc) {
#                say "\$misc_feature isa '" . ref($misc_feature) . "'";
#                if ($misc_feature->has_tag('note')) {
#                    my @notes = $misc_feature->get_tag_values('note');
#                    for my $note (@notes) {
#                        if ($note =~ /(COG\d\d\d\d)/) {
#                            $cog = $1;
#                        }
#                    }
#                }
#            }

        say $PTT $start . ".." . "$stop\t$strand\t$length\t$pid[0]\t$gene[0]\t$synonym[0]\t$code\t$cog\t$description[0]";
    }

    #foreach my $feature (sort { $a->start() <=> $b->start() } (@tRNA, @rRNA)) {
    for my $feature (sort { $a->start() <=> $b->start() } @RNA) {

            my $start = $feature->start();
            my $stop = $feature->end();
            my $strand = '+';
            my $length = $feature->length();
            my $pid = $feature->has_tag('db_xref') ? ($feature->get_tag_values('db_xref'))[0] : ($feature->get_tag_values('locus_tag'))[0];
            my @gene = $feature->has_tag('gene') ? $feature->get_tag_values('gene') : '-';
            my @synonym = $feature->get_tag_values('locus_tag');
            my $code = '-';
            my $cog = '-';
            my @description = $feature->get_tag_values('product');

            say $RNT $start . ".." . "$stop\t$strand\t$length\t$pid\t$gene[0]\t$synonym[0]\t$code\t$cog\t$description[0]";
        }
}

