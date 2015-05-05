#!/usr/bin/env perl
# $Id: seqsplit.pl,v 1.7 2010/02/25 19:35:27 givans Exp $

use Bio::SeqIO;
use Getopt::Std;
use warnings;
use strict;
use vars qw/$opt_f $opt_F $opt_d $opt_o $opt_O $opt_n $opt_x $opt_h/;

my $descr = "
#
# This is a program to split a file with many fasta sequences into several
# files with a user-defined number of sequences in each output file 
# (except for possibly the last file.
#
# usage:  seqsplit.pl <options>
";

if (!$ARGV[0]) {
  print "$descr\n";
  exit(0);
}

getopts('f:F:d:o:O:n:x:h');
my ($file,$infile_format,$file_out,$dir,$out,$outfile_format,$cnt,$cntall,$outcnt,$splitnum,$seqin,$seqout) = ();

if ($opt_h) {

print <<HELP;
$descr

Command line options:

-f	name of input file
-F	input file format [default=fasta]
-d	name of output directory to place files
-o	output file name
-O	output file format
-n	maximum number of sequences in each file
-x	stop after parsing this number of sequences from input file
-h	print this help menu


HELP

exit(0);
}

#
# Read input file name
#
if (!$opt_f) {
  print "What is the name of the input file? ";
  chomp($file = <STDIN>);
  die "'$file' isn't a valid filename\n" unless (&namecheck($file));
  die "file doesn't exist\n" unless (&file_e($file));
} else {
  $file = $opt_f;
  $infile_format = $opt_F || 'fasta';
}

#
# Read number of sequences per output file
if ($opt_n) {
  $splitnum = $opt_n;
} else {
  print "How many sequences do you want in each output file? ";
  chomp($splitnum = <STDIN>);
}
die "you must enter a value for the number of sequences per output file\n" unless ($splitnum =~ /[\d]/);

#
# Read output filename root
#
if ($opt_o) {
  die "'$opt_o' isn't a valid file name\n" unless (&namecheck($opt_o));
  $file_out = $opt_o;
} else {
  $file_out = $file;
}
$outfile_format = $opt_O || 'fasta';
#
# Read output directory
#
if ($opt_d) {
    if (&namecheck($opt_d)) {
            if (&file_e($opt_d)) {
                $dir = $opt_d;

                if ($file_out =~ /^.+\/([\w\.\_\-]+)$/) {
                    $file_out = $1;
                }
                $file_out = "$dir" . "/" . "$file_out";
            } else {
                die "'$opt_d' doesn't exist\n";
            }
    } else {
        die "'$opt_d' isn't a valid file name\n";
    }
}


$seqin = Bio::SeqIO->new(
						-file	=>	$file,
#						-format	=>	$infile_format,
						-format	=>	'genbank',
				);
				
				
$seqout = Bio::SeqIO->new(
    #-format		=>	$outfile_format,
							-format		=>	'genbank',
							# print sequence ID on quality lines also:
                            # -quality_header	=>	1, # only affects fastq output
						);
					

while (my $seq = $seqin->next_seq()) {
  ++$cnt;
  ++$cntall;
  
   if ($cnt <= $splitnum && $cnt != 1) {
# 
   } else {
        $cnt = 1;
        ++$outcnt;
        #my $new_file_out = "$file_out" . "_$outcnt";
        #my $new_file_out = $seq->accession() . "." . $seq->version() . ".gbk";
        my $new_file_out = $seq->accession() . ".gbk";
        my $newFH;
        open($newFH, ">$new_file_out") or die "can't open '$new_file_out': $!";
        $seqout->_fh($newFH);

   }
   $seqout->write_seq($seq);

  if ($opt_x) {
    exit if ($cntall == $opt_x);
  }
}


sub namecheck {
  my $name = shift;

#  print "checking '$name'\n";
  if ($name =~ /[\w\d\.\-_]/) {
    return 1;
  } else {
    return 0;
  }
}

sub file_e {
  my $file = shift;

  if (-e $file) {
    return 1;
  } else {
    return 0;
  }
}

END {
  if ($cntall && $splitnum) {
    print "$cntall sequences split into individual files of $splitnum sequences each\n";
  }
  exit(0);
}

