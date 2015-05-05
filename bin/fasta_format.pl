#!/usr/bin/env perl
# $Id: fasta_format.pl,v 1.21 2007/03/24 00:05:49 givans Exp $
# $Log: fasta_format.pl,v $
# Revision 1.21  2007/03/24 00:05:49  givans
# added swissprot and kegg as possible input formats.
#
# Revision 1.20  2004/10/05 22:15:54  givans
# also remove ascii <= 31 from descriptions
#
# Revision 1.19  2004/10/05 21:35:08  givans
# added section in cleanDescr() to remove ascii characters > 127
#
# Revision 1.18  2004/10/05 19:47:40  givans
# added tab as acceptable file format
#
# Revision 1.17  2004/10/02 00:05:11  givans
# added embl as an accepted format
#
# Revision 1.16  2004/09/28 22:17:44  givans
# added to help message to explain -P # usage
#
# Revision 1.15  2004/09/28 22:16:19  givans
# added option to prefix each id with a incremented number - makes each ID unique
#
# Revision 1.14  2004/09/23 20:19:02  givans
# $opt_R was in a separate use vars statement so I moved it into the main use vars
# added $opt_v flag for verbose output to terminal
#
# Revision 1.13  2004/09/23 19:54:46  givans
# accepts phd as output format
#
# Revision 1.12  2004/09/23 19:54:17  givans
# accepts 'phd' as file input format
#
# Revision 1.11  2004/07/26 21:25:33  givans
# added some characters to replace when -I is used
#
# Revision 1.10  2004/07/26 21:19:16  givans
# *** empty log message ***
#
# Revision 1.9  2004/04/01 20:58:37  givans
# Added ability to provide a regular expression to use in replacing the sequence
# identifier with a portion of the current sequence identifier
#
# Revision 1.8  2004/03/12 19:48:57  givans
# Added -I option to clean problem characters fro sequence identifiers
#
# Revision 1.7  2004/03/11 20:06:06  givans
# Added some comments
# Chaged how $rootprefix is made
#
# Revision 1.6  2004/03/11 17:25:15  givans
# Added an eval{} around $seq->sequence(), which was causing problems.
# Also added a -w option, which will print warnings to a file
#
# Revision 1.5  2004/03/06 00:17:49  givans
# Added some new capabilities
# -f all will process all files in current directory
# -F appends file name root to all sequence names
#
# Revision 1.4  2004/01/15 19:24:38  givans
# added -D option to clean the description line
# added Log line to header
#
#
#use lib '/local/lib/perl5';
use strict;
use warnings;
use Carp;
use Bio::SeqIO;
use Getopt::Std;
use vars qw($opt_f $opt_i $opt_o $opt_O $opt_a $opt_h $opt_d $opt_D $opt_P $opt_S $opt_n $opt_F $opt_w $opt_A $opt_I $opt_R $opt_v);

getopts('f:i:o:O:ahdDP:S:nFw:AIR:v');

my $usage = "fasta_format.pl <options>";

if ($opt_h) {
print <<HELP;

This script converts an input DNA or protein file in a standard
format into another standard format.  Standard formats are fasta,
gcg, staden and genbank.

 usage:  $usage

Command line options

-f	input file name [use 'all' to process all files in current directory]
-i	input file format (defaults to 'fasta')
-o	desired output format (defaults to 'fasta')
-O	output file name
-A	append sequences to this file
-a	convert complex id to simple accession number
-I	clean nasty characters from sequence ID's
-D	clean the description - remove 'bad' characters
-P	prefix all sequence names with ( # prefixes unique number to each ID)
-F	prefix all sequence names with root of file name
-S	suffix all sequence names with
-R	regular expression to use for replacing sequence ID
-n	remove description (only use seq id)
-h	print this help menu
-d	debugging enabled:  verbose output
-w	print warnings to this file
-v	verbose output to terminal

HELP
exit(0);
}

#################################
#  Variable Declarations	#
#################################

my($file,@files,$format_in,$format_out,$file_in,$file_out,$prefix,$suffix,$warnfile,$RE) = "";

#################################
#	Collect User Input	#
#################################

if (!$opt_f) {
    print "What file do you want to convert to fasta format?: ";
    $file_in = <STDIN>;
    chomp $file;
} else {
    $file_in = $opt_f;
}

if ($file_in eq 'all') {
  opendir(DIR, '.') or croak("can't open this directory");
  @files = readdir(DIR);
  closedir(DIR);
} else {
  push(@files,$file_in);
}

# foreach my $file_tmp (@files) {
#   print "file root: '$file_tmp'\n";
# }
# exit(0);

if ($opt_i) {
  if ($opt_i =~ /[\w]/) {
#    if ($opt_i eq 'raw' || $opt_i eq 'fasta' || $opt_i eq 'gcg' || $opt_i eq 'staden' || $opt_i eq 'genbank' || $opt_i eq 'phd' || $opt_i eq 'embl' || $opt_i eq 'tab' || $opt_i eq 'swiss' || $opt_i eq 'kegg') {
      $format_in = $opt_i;
#    } else {
#      croak "'$opt_i' is not accepted";
#    }
  } else {
    croak "'$opt_i' isn't a valid format name";
  }
} else {
  $format_in = 'fasta';
}

if ($opt_o) {
  if ($opt_o =~ /[\w]/) {
#    if ($opt_o eq 'raw' || $opt_o eq 'fasta' || $opt_o eq 'gcg' || $opt_o eq 'staden' || $opt_o eq 'genbank' || $opt_o eq 'phd' || $opt_o eq 'embl' || $opt_o eq 'tab') {
      $format_out = $opt_o;
#    } else {
#      croak "'$opt_o' is not accepted";
#    }
  } else {
    croak "'$opt_o' isn't a valid format name";
  }
} else {
  $format_out = 'fasta';
}

if ($opt_O) {
  $file_out = $opt_O;
}# else {
#  $file_out = "$file" . ".out";
#}

if ($opt_P) {
  $prefix = $opt_P;
}

if ($opt_S) {
  $suffix = $opt_S;
}

if ($opt_w) {
  $warnfile = $opt_w;
  open(WARN,">$warnfile") or croak("can't open warnings file '$warnfile': $!");
}

if ($opt_R) {
  $RE = $opt_R;
}

print "$0 -f $file -i $format_in -O $file_out -o $format_out\n" if ($opt_d);

#################################
#  Perform Actions		#
#################################
my $cnt;
foreach my $file (@files) {

  next unless ($file !~ /^\./);
  my ($output_file,$rootprefix);
  ++$cnt;
  if (!$file_out) {
    $output_file = "$file" . ".out";
  } else {
    $output_file = $file_out;
  }

  if ($opt_F) {
    print "appending root file of filename to prefix\n";
     $rootprefix = $file;
#     $rootprefix =~ s/\..+//;
     if ($rootprefix =~ /^(.+)\.\w+?$/) {
       $rootprefix = $1 . "_";
     } else {
       $rootprefix = "untitled$cnt" . "_";
     }
   }

  my $seqio_in = Bio::SeqIO->new(	-file	=>	$file,
					-format	=>	$format_in,
				);
  if ($opt_A) {
    $output_file = ">>$output_file";
  } else {
    $output_file = ">$output_file";
  }
  my $seqio_out = Bio::SeqIO->new(	-file	=>	"$output_file",
					-format	=>	$format_out,
				 );
  my $seq_cnt = 0;

  while (my $seq = $seqio_in->next_seq()) {
    my $sequence = $seq->seq();### getting current sequence
    ++$seq_cnt;

    $sequence = uc($sequence);### converting sequence to uppercase

    eval { $seq->seq($sequence); };### setting sequence to the uppercase version
    ### wrapped in an eval to catch errors but continue processing
    if ($@) {### if $sequence is empty
      warn($seq->id() . " has no sequence");
      print WARN $seq->id() . " has no sequence\n" if ($opt_w);
      next;
    }

    if ($opt_R) {
      if ($seq->id =~ /$RE/) {
	$seq->id($1);
      }
    }

    if ($opt_a) {
      my $simpleID = $seq->primary_id();
      if ($simpleID =~ /ref\|([\w\._]+?)\|/) {
	$simpleID = $1;
      } elsif ($simpleID =~ /gi\|(\w+?)\|/) {
	$simpleID = $1;
      } elsif ($simpleID =~ /lcl\|(\w+?)\|/) {
	$simpleID = $1;
      } else {
	$simpleID = 'unavailable';
      }
      print "simpleID = '$simpleID'\n" if ($opt_d);
      $seq->id($simpleID);
    }

    if ($opt_I) {
      cleanID($seq);
    }

    if ($prefix || $suffix || $rootprefix) {
      my $tempID = $seq->id();
      my $prefix_temp;
      if ($prefix) {
	if ($prefix eq '#') {
	  $prefix_temp = "$seq_cnt" . "_";
	} else {
	  $prefix_temp = $prefix;
	}
	$tempID = "$prefix_temp$tempID";
      }
      $tempID = "$rootprefix$tempID" if ($rootprefix);
      $tempID .= $suffix if ($suffix);
      $seq->id($tempID);
    }

    if ($opt_D) {
      cleanDescr($seq);
    }

    if ($opt_n) {
      $seq->description('');
    }

    $seqio_out->write_seq($seq);
  }
}

if ($opt_w) {
  close(WARN);
}

print "$cnt sequence files processed\n" if ($opt_v);

sub cleanDescr {
  my $seq = shift;
  my $mod = 0;
  my $desc = $seq->description();
  my $temp;

#  $temp = 1 if ($desc =~ /68415\.m02594/);

  my @chars = split //,$desc;
  for (my $i = 0; $i < scalar(@chars); ++$i) {
    if (ord($chars[$i]) >= 128 || ord($chars[$i]) <= 31) {
      print "char: '$chars[$i]' ord: '", ord($chars[$i]), "'\n" if ($opt_v);
      $chars[$i] = " ";
      $mod = 1;
    }
  }
  $desc = join "", @chars;
  
  if ($desc) {
    if ($desc =~ /[>#;]/) {
      $desc =~ s/[>#;]/ /g;
      $mod = 1;
    }

    if ($mod) {
      $seq->description($desc);
    }
  }

}

sub cleanID {
  my $seq = shift;
  my $mod = 0;

  my $id = $seq->id();
  if ($id =~ /[\.\*\#\(\)\~\\\/\-\'\"]/) {
    $id =~ s/[\.\*\#\(\)\~\\\/\-\'\"]/_/g;
    $mod = 1;
  }

  if ($mod) {
    $seq->id($id);
  }
}

