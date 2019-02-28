# Ntthal-XS

This Perl module uses [Inline::C](https://metacpan.org/pod/Inline::C) 
to bind the ntthal alignment from the [primer3](https://github.com/primer3-org/primer3) 
code and make it available as a Perl subroutine call. It uses an object-oriented 
[Mouse](https://metacpan.org/pod/Mouse) interface.

The main motivation to do this was speed of execution, especially to _reduce overhead_ 
generated by process creation (fork) and passing of sequences, retrieval of results 
via shell process communication. See _BENCHMARK_ for results.

## INSTALLATION

To install this module, clone the code and run the following commands:

	perl Makefile.PL
	make
	make test
	make install

Alternatively, you can download a [release](/releases) and run:

	cpanm Ntthal-XS-#.#.#.tar.gz

Here, the # signs are placeholders for the version number.

The package will also install a demo/driver script named 
[ntthal_xs.pl](/bin/ntthal_xs.pl), which should be executable 
after isntallation. The purpose is to test the functionality 
and serve as a template for own code.

## SYNOPSIS

```perl
use Ntthal::XS;
my $NT = Ntthal::XS->new();

# tm will accept only one sequence
print $NT->tm('TAATACGACTCACTATAGGG'); # 44.5215658215416

# change parameters
$NT->type('END1');
$NT->dv(0.8);

# ntthal will use two sequences (different strands!)
my $hybrid = $NT->ntthal('CCCGAAAAGTGCCACCTG', 'CAGGTGGCATTTTTTCGGG');

# $hybrid = {
#    dG => -13335.3624649896,     # { 0}
#    dH => -132500,               # { 1}
#    dna_conc => 50,              # { 2}
#    dntp => 0.8,                 # { 3}
#    dS => -384.216145526392,     # { 4}
#    dv => 0.8,                   # { 5}
#    maxLoop => 8,                # { 6}
#    mv => 50,                    # { 7}
#    s1 => "CCCGAAAAGTGCCACCTG",  # { 8}
#    s2 => "CAGGTGGCATTTTTTCGGG", # { 9}
#    t => 42.0422986972263,       # {10}
#    temp => 37,                  # {11}
#    type => "END1",              # {12}
#  }
 
# compare to the ntthal command line output
#
#  ntthal -a END1 -dv 0.8 -s1 CCCGAAAAGTGCCACCTG -s2 CAGGTGGCATTTTTTCGGG
#  Calculated thermodynamical parameters for dimer:	dS = -384.216	dH = -132500	dG = -13335.4	t = 42.0423
#  SEQ	        G-         
#  SEQ	CCCGAAAA  TGCCACCTG
#  STR	GGGCTTTT  ACGGTGGAC
#  STR	        TT         
```

## BENCHMARK

Here, a comparison is shown for the traditional **ntthal** command line call, a newer version that 
does not read all parameter files at every call (it uses pre-populated lists) named **ntthal_new** 
and the **xs** version implemented here. *Note:* The underlying algorithm was __not changed__, only the 
overhead calling it from Perl was reduced.

```
Benchmark: running ntthal, ntthal_new, xs for at least 3 CPU seconds...

     ntthal: 23 wallclock secs ( 0.46 usr  2.57 sys + 14.47 cusr  5.04 csys = 22.54 CPU) @ 106.70/s (n=2405)
 ntthal_new: 18 wallclock secs ( 0.34 usr  2.84 sys +  8.99 cusr  2.31 csys = 14.48 CPU) @ 389.36/s (n=5638)
         xs:  3 wallclock secs ( 3.13 usr +  0.00 sys =  3.13 CPU) @ 1490.42/s (n=4665)

             Rate     ntthal ntthal_new         xs
 ntthal      107/s         --       -73%       -93%
 ntthal_new  389/s       265%         --       -74%
 xs         1490/s      1297%       283%         --
```
That is about a 4-fold improvement on my system.

## LICENSE AND COPYRIGHT

Copyright (C) 2019 Frank Reinecke (frank.reinecke@qiagen.com)

This program is released under the following license: AGPL3
