#!/usr/bin/perl -w
# process_fa.pl --- change the format of the first line for fasta files
# Author: Fei Zhan <fei@fei-laptop>
# Created: 3 Dec 2012
# Version: 0.01
# Modified by Yanghua at 26 Aug 2013
# Modified by: Alexis Marceau
# Modified: 29 Oct 2021
# Being used in RNAseq to lncRNA pipeline development

use warnings;
use strict;

print "Condition name, please be consistent with previous answers:\n";
$Condition_name = <>;
chomp($Condition_name);

my $filein = '$Condition_name.intergenic_loci.fa';
my $fileout = '$Condition_name.intergenic_loci_2.fa';

my $usage = "Uasge:\n\t perl process_fa.pl -i [fasta] -o [fasta]\n\n";
foreach my $i (0 .. scalar(@ARGV)-1) {
   if($ARGV[$i] eq '-i') {
      $filein = $ARGV[++$i];
   }elsif($ARGV[$i] eq '-o') {
      $fileout = $ARGV[++$i];
   }
}  
if(@ARGV ==0) {
    die $usage;
}

open(IN,$filein) or die "Cannot open $filein.\n";
open(OUT,'>',$fileout) or die "Cannot open $fileout.\n";

while(<IN>) {
    chomp;
    my $line = $_;
    if($line =~ m/^>/) {
	my($loci, $name) = (split(" ",$line))[1,2];
	print OUT ">$name $loci\n";
    }else{
	print OUT $line,"\n";
    }
}

close IN;
close OUT;

