#!/usr/bin/perl
# Extracts ONLY the character sequence from manually-created file
# INPUTS
# ARGV[0]: "organism_seqfile.txt" = file containing the sequence part of the Genbank record (from begin to end of sequence lines, exclude ORIGIN and //) 

use strict;
# use warnings;
# declare and initialize variables
my $seq = '';
my @seq = ( );

# print "\nFirst Argument: "; print $ARGV[0];

#-------------------------------------------------------
# Input file
my $infile = $ARGV[0];
#-------------------------------------------------------

# Read from file
 unless ( open(INFILE, "$infile") ) {
     print "Cannot open file \"$infile\" for reading!!\n\n";
     exit;
 }
 @seq = <INFILE>; 
 close(INFILE);

# Clean-up
$seq = join('',@seq);
$seq =~ s/[\s0-9]//g; # remove whitespace and line numbers from sequence

print $seq;

exit;                        
