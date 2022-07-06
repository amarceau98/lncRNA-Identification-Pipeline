#!/usr/bin/perl -w
###########################################################################
# summary_gtf.pl --- Process gtf files
# <transcript_id> <gene_id> <chr> <start> <end> <len> <# of exons> <len_exon1;len_exon2;...>
###########################################################################
# Author: Fei Zhan <fei@fei-laptop>
# Created: 19 Dec 2012
# Version: 0.01
# Modified by Yanghua at 14 Sep. 2013
# Modified by: Alexis Marceau
# Modified: 29 Oct 2021
# Being used in RNAseq to lncRNA pipeline development

use warnings;
use strict;

my $filein = "gtf_in";
my $fileout = "file_out";

my $min_length = 200;
my $min_exons = 2;

my $usage = "\nUsage:\n\t perl summary_gtf.pl -i [filein.gtf] -o [fileout.txt]\n\n";
foreach my $i (0 ..scalar(@ARGV)-1) {
  if($ARGV[$i] eq '-i') {
    $filein = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-o') {
    $fileout = $ARGV[++$i];
  }
}
if(@ARGV ==0) {
    die $usage;
}

open(IN,'<',$filein) or die "Cannot open $filein.\n$usage";
open(OUT,'>',$fileout) or die "Cannot open $fileout.\n$usage";

my $headline = join("\t", ("TCONS", "XLOC", "chr", "strand", "start", "end", "num_exons", "length", "starts", "ends"));
print OUT $headline,"\n";

my %hash_gtf;
while(<IN>) {
  chomp;
  my @terms = split("\t",$_);
  my($chr,$start,$end,$strand,$info) = @terms[0,3,4,6,8];
  my $gene_id = $1 if $info =~ /gene_id "(.+?)";/;
  my $transcript_id = $1 if $info =~ /transcript_id "(.+?)";/;
  my $len_exon = $end-$start+1;
  if(exists $hash_gtf{$transcript_id}) {
    my $ref = $hash_gtf{$transcript_id};
    if($transcript_id ne $$ref[0] or $gene_id ne $$ref[1] or $chr ne $$ref[2] or $strand ne $$ref[3]) {
      print $_;
      die "Error!\n";
    }
    push(@{$$ref[4]}, $start);
    push(@{$$ref[5]}, $end);
    push(@{$$ref[6]}, $len_exon);
  }else{
    my @lens = ($len_exon);
    my @starts = ($start);
    my @ends = ($end);
    my $data = [$transcript_id,$gene_id,$chr,$strand,\@starts,\@ends,\@lens];
    $hash_gtf{$transcript_id} = $data;
  }
}

foreach my $key (keys %hash_gtf) {
  my $ref = $hash_gtf{$key};
  my @starts = @{$$ref[4]};
  my @ends = @{$$ref[5]};
  my $num_exons = scalar(@{$$ref[6]});
  my $len = 0;
  foreach (@{$$ref[6]}) {
    $len += scalar($_);
  }
  if($len >= 200 and $num_exons >= 2) {
    print OUT join("\t", (@{$hash_gtf{$key}}[0,1,2,3], $starts[0], $ends[-1], $num_exons, $len, join(",", @starts), join(",", @ends))),"\n";
  }
#  print OUT join("\t", (@{$hash_gtf{$key}}[0,1,2,3], $len, $num_exons)),"\n";
}


close IN;
close OUT;



__END__

=head1 NAME

process_gtf.pl - Describe the usage of script briefly

=head1 SYNOPSIS

process_gtf.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for process_gtf.pl, 

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
