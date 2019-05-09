#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dump::Color qw<dd>;
use Ntthal::XS;
my $Ntthal = Ntthal::XS->new();
my $seq1   = shift;
my $seq2   = shift;
usage() unless $seq1;

my ($result, $cmd, $exe);
chomp($exe = `which ntthal` || '# <ntthal not found> ');

for my $mode ('ANY', 'END1', 'END2', 'HAIRPIN') {
  # change the type, and run again
  $Ntthal->type($mode);
  
  $result = $Ntthal->ntthal( $seq1, $seq2 );
  dd $result;
  
  $cmd = sprintf "%s -s1 %s -s2 %s\n", $Ntthal->args('command'), $result->{s1}, $result->{s2};
  if ($mode eq 'HAIRPIN') {
    $cmd = sprintf "%s -s1 %s\n", $Ntthal->args('command'), $result->{s1};
  }
  print "\n$exe $cmd\n";
  if (-x $exe) {
    print `$exe $cmd`, "\n";
  }
}

my $Test = Ntthal::XS->new('mv' => 30, 'only' => 'tm' );
dd $Test->ntthal($seq1);

sub usage {
  print <<"END_OF_USAGE"

Usage:
  
  $0 <seq1> [ <seq2> ]
  
  This script is just a driver/testing script for the XS binding 
  of the ntthal C routines that are part of the primer3 code base.
  
  Sequences are hybridized as they are, so you need to provide 
  the strand that actually form the hybrid structure.
  
  If the second sequence is omitted, the perfect matching reverse 
  complement of the first will be generated to calculate the 
  hybrid's values. In that case, the maxLoop parameter will be set 
  to zero (0).

END_OF_USAGE
  ;
  exit;
};

