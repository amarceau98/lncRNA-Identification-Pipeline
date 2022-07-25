#!/usr/bin/perl -w
# get_intergenicLoci_list.pl --- get transcripts list of common intergenic loci in both of cmp_ens.gtf and cmp_ref.gtf based on class_code "u"
# Author: Yanghua
# Created: 20 Aug 2013
# Modified by: Alexis Marceau
# Modified: 23 May 2022
# Being used in RNAseq to lncRNA pipeline development
# Version: 0.01

use warnings;
use strict;

my $loci = 'loci_replace_ref1';
my $gtf = 'loci_replace_ref2';
my $fileout = 'list_replace';

foreach my $i (0 .. scalar(@ARGV)-1) {
  if($ARGV[$i] eq '-g1') {
    $gtf = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-g2') {
    $loci = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-o') {
    $fileout = $ARGV[++$i];
  }
}
if($loci eq '' | $gtf eq '' | $fileout eq '') {
  die "Usage:\n\t perl get_intergenic.pl -g1 [ens_gtf] -g2 [ref_gtf] -o [output]\n\n";
}

open(LOCI,$loci) or die "Cannot open $loci.\n\nUsage:\n\t perl get_intergenic.pl -g1 [ens_gtf] -g2 [ref_gtf] -o [output]\n\n";
open(GTF,$gtf) or die "Cannot open $gtf.\n\nUsage:\n\t perl get_intergenic.pl -g1 [ens_gtf] -g2 [ref_gtf] -o [output]\n";
open(OUT,'>',$fileout) or die "Cannot open $fileout.\n\nUsage:\n\t perl get_intergenic.pl -g1 [ens_gtf] -g2 [ref_gtf] -o [output]\n";

my %hash_loci_u;
my %hash_loci_nou;
while(<LOCI>) {
  chomp;
  if($_ =~ /.+transcript_id "(.+?)";.+class_code "(.+?)";.+/) {
    if($2 eq "u") {
	  if(!exists $hash_loci_nou{$1}){
        $hash_loci_u{$1} = 1;
	  }
    }else{
	  if(exists $hash_loci_u{$1}){
	    delete $hash_loci_u{$1};
	  }
	  $hash_loci_nou{$1} = 1
  }
}
}
close LOCI;

my %hash_gtf_u;
my %hash_gtf_nou;
while(<GTF>) {
  chomp;
  if($_ =~ /.+transcript_id "(.+?)";.+class_code "(.+?)";.+/) {
    if($2 eq "u") {
	  if(!exists $hash_gtf_nou{$1}){
        $hash_gtf_u{$1} = 1;
	  }
    }else{
	  if(exists $hash_gtf_u{$1}){
	    delete $hash_gtf_u{$1};
	  }
	  $hash_gtf_nou{$1} = 1
  }
}
}
foreach my $key (sort keys %hash_gtf_u){
    if(exists $hash_loci_u{$key}){
	  print OUT $key,"\n";
	}
}
close GTF;
close OUT;











__END__

=head1 NAME

get_intergenic.pl - Describe the usage of script briefly

=head1 SYNOPSIS

get_intergenic.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for get_intergenic.pl, 

=head1 AUTHOR

Fei Zhan, E<lt>fei@fei-laptopE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by Fei Zhan

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
