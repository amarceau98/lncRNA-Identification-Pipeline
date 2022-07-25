#!/usr/bin/perl -w
# get_intergenic.pl --- Get intergenic transcript information
# Author: Fei Zhan <fei@fei-laptop>
# Created: 04 Dec 2012
# Modified by: Alexis Marceau
# Modified: 23 May 2022
# Being used in RNAseq to lncRNA pipeline development
# Version: 0.01

use warnings;
use strict;

my $loci = 'loci_replace';
my $gtf = 'gtf_replace';
my $fileout = 'out_replace';

foreach my $i (0 .. scalar(@ARGV)-1) {
  if($ARGV[$i] eq '-g') {
    $gtf = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-l') {
    $loci = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-o') {
    $fileout = $ARGV[++$i];
  }
}
if($loci eq '' | $gtf eq '' | $fileout eq '') {
  die "Usage:\n\t perl get_intergenic.pl -g [gtf] -l [loci] -o [output]\n\n";
}

open(LOCI,$loci) or die "Cannot open $loci.\n\nUsage:\n\t perl get_intergenic.pl -g [gtf] -l [loci] -o [output]\n\n";
open(GTF,$gtf) or die "Cannot open $gtf.\n\nUsage:\n\t perl get_intergenic.pl -g [gtf] -l [loci] -o [output]\n";
open(OUT,'>',$fileout) or die "Cannot open $fileout.\n\nUsage:\n\t perl get_intergenic.pl -g [gtf] -l [loci] -o [output]\n";

my %hash_loci;
while(<LOCI>) {
  chomp;
  $hash_loci{$_}++;
}
close LOCI;

while(<GTF>) {
  chomp;
  my($attributes) = (split("\t",$_))[8];
  if($attributes =~ m/transcript_id "(.+?)";/) {
    my $loci = $1;
    if(exists $hash_loci{$loci}) {
      print OUT $_,"\n";
    }
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
