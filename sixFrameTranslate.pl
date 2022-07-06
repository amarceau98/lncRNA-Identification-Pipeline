#!/usr/bin/perl -w
# sixFrameTranslate.pl --- Six Frame Translation
# Author: Fei Zhan <fei@fei-laptop>
# Created: 11 Feb 2013
# Version: 0.01
# Modified by: Alexis Marceau
# Modified: 29 Oct 2021
# Being used in RNAseq to lncRNA pipeline development

use warnings;
use strict;

my $codon = 'codon.txt';
my $filein = 'file_in_replace';
(my $fileout = $filein) =~ s/\.fa$/\.ff/;
my $min_length = 6;

foreach my $i (0 .. scalar(@ARGV)-1) {
  if($ARGV[$i] eq '-i') {
    $filein = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-o') {
    $fileout = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-l') {
    $min_length = $ARGV[++$i];
  }
}


my $usage = "Usage:\n\tsixFrameTranslate.pl -i [input] -o [output] -l [min_length]\n\n";
if(@ARGV ==0) {
    die $usage;
}

my %hash_codon = readCodon($codon);

open(IN, $filein) or die "Cannot open $filein.\n$usage";
open(OUT,'>',$fileout) or die "Cannot open $fileout.\n$usage";
local $/ = ">";
my $pattern = qr/(.+?)\n(.+)/ms;
while(<IN>) {
  if($_ =~ /$pattern/) {
    my $name = $1;
    my $seq = $2;
    $seq =~ s/\s+//g;
    $seq =~ s/>$//; #remove tailing '>' if any
#    print OUT $name,"\n",$seq,"\n";
    my $frame_1 = translate(+1, $seq);
    my $frame_2 = translate(+2, $seq);
    my $frame_3 = translate(+3, $seq);
#    my $frame_r1 = translate(-1, reverse_complement($seq));
#    my $frame_r2 = translate(-2, reverse_complement($seq));
#    my $frame_r3 = translate(-3, reverse_complement($seq));
#    print OUT join("\n", ($frame_1,$frame_2,$frame_3)),"\n";
#    print OUT join("\n", ($frame_r1,$frame_r2,$frame_r3)),"\n";
    $name =~ s/ /\|/;
    print OUT ">$name|f1\n",$frame_1,"\n";
    print OUT ">$name|f2\n",$frame_2,"\n";
    print OUT ">$name|f3\n",$frame_3,"\n";
  }
}

close IN;
close OUT;


sub reverse_complement {
  my $dna = shift;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}


sub translate{
  my($frame, $seq) = @_;
  my @aa;
  for(my $pos = abs($frame)-1; $pos < length($seq)-2; $pos += 3) {
    my $codon = uc(substr($seq, $pos, 3));
    if(exists $hash_codon{$codon}) {
      push(@aa, $hash_codon{$codon});
    }else{
      push(@aa, 'X');
#      die "Unrecognized codon: $codon.\n";
    }
  }
  my $aa = join('', @aa);
  return $aa;

}






# Read in codon
sub readCodon {
  my($codon) = @_;
  my %codon;
  open(IN,$codon) or die "Cannot open $codon.\n";
  while(<IN>) {
    my($base,$aa) = split("\t",$_);
    if(exists $codon{$base}) {
      die "Error! Duplicate codon: $base.\n";
    }else{
      $codon{$base} = $aa;
    }
  }
  close IN;
  return %codon;
}





__END__

=head1 NAME

sixFrameTranslate.pl - Describe the usage of script briefly

=head1 SYNOPSIS

sixFrameTranslate.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for sixFrameTranslate.pl, 

=head1 AUTHOR

Fei Zhan, E<lt>fei@fei-laptopE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Fei Zhan

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
