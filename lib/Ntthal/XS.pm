package Ntthal::XS;
our $VERSION = '0.2.7';

#>>>
use Ntthal::XS::Inline C => 'DATA';
#<<<
use Mouse;
use Method::Signatures;



=head1 NAME

Ntthal::XS - Perl bindings to ntthal C code

=head1 VERSION

Version 0.2.7

=head1 SYNOPSIS

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
    #   dG => -13335.3624649896,          # { 0}
    #   dH => -132500,                    # { 1}
    #   dna_conc => 50,                   # { 2}
    #   dntp => 0.8,                      # { 3}
    #   dS => -384.216145526392,          # { 4}
    #   duplex0 => "        G-         ", # { 5}
    #   duplex1 => "CCCGAAAA  TGCCACCTG", # { 6}
    #   duplex2 => "GGGCTTTT  ACGGTGGAC", # { 7}
    #   duplex3 => "        TT         ", # { 8}
    #   dv => 0.8,                        # { 9}
    #   maxLoop => 8,                     # {10}
    #   mv => 50,                         # {11}
    #   s1 => "CCCGAAAAGTGCCACCTG",       # {12}
    #   s2 => "CAGGTGGCATTTTTTCGGG",      # {13}
    #   t => 42.0422986972263,            # {14}
    #   temp => 37,                       # {15}
    #   type => "END1",                   # {16}
    # }
    
    # compare to the ntthal command line output
    #
    #  ntthal -a END1 -dv 0.8 -s1 CCCGAAAAGTGCCACCTG -s2 CAGGTGGCATTTTTTCGGG
    #  Calculated thermodynamical parameters for dimer:	dS = -384.216	dH = -132500	dG = -13335.4	t = 42.0423
    #  SEQ	        G-         
    #  SEQ	CCCGAAAA  TGCCACCTG
    #  STR	GGGCTTTT  ACGGTGGAC
    #  STR	        TT         

=head1 ARGUMENTS

=head2 ALIGNMENT TYPE

=over 3

=item type - alignment type, END1, END2, ANY and HAIRPIN, by default ANY (when duplex)

=back

=head2 HYBRIDIZATION CONDITIONS

=over 3

=item mv - monovalent_conc = concentration of monovalent cations in mM, by default 50 mM

=item dv - divalent_conc = concentration of divalent cations in mM, by default 0 mM

=item dntp - concentration of deoxynycleotide triphosphate in mM, by default 0 mM

=item dna_conc - concentration of DNA strands in nM, by default 50 nM

=item maxLoop - the maximum size of secondary structures loops. 

=item temp - temperature (in Celsius) at which duplex is calculated, by default 37

=back

=head2 RETURN VALUES / INTERNALS

=over 3

=item only - do not return a hashref, but only the tm (^t) or dG (^[dg]) scalar value

=item typemap - used internally to convert named mode into integers

=item counter - used internally to count the number of calls

=back


=cut


# type - alignment type, END1, END2, ANY and HAIRPIN, by default ANY (when duplex)
has 'type' => (
  is => 'rw',
  isa => 'Str',
  default => 'ANY'
);

# mv - monovalent_conc = concentration of monovalent cations in mM, by default 50 mM
has 'mv' => (
  is => 'rw', 
  isa => 'Num',
  default => 50
);

# dv - divalent_conc = concentration of divalent cations in mM, by default 0 mM
has 'dv' => (
  is => 'rw', 
  isa => 'Num',
  default => 0
);
# dntp - concentration of deoxynycleotide triphosphate in mM, by default 0 mM
has 'dntp' => (
  is => 'rw', 
  isa => 'Num',
  default => 0.8
);

# dna_conc - concentration of DNA strands in nM, by default 50 nM
has 'dna_conc' => (
  is => 'rw', 
  isa => 'Num',
  default => 50
);

# maxLoop - the maximum size of secondary structures loops. 
## NOTE: The default is 8, which is different from ntthal, where the 
## default is 30 (which is maximum allowed length, currently).
has 'maxLoop' => (
  is => 'rw', 
  isa => 'Num',
  default => 8
);

# temp - temperature (in Celsius) at which duplex is calculated, by default 37
has 'temp' => (
  is => 'rw', 
  isa => 'Num',
  default => 37
);

# only - do not return a hashref, but only the tm (^t) or dG (^[dg]) scalar value
has 'only' => (
  is => 'rw',
  isa => 'Str',
  default => ''
);

# typemap - used internally to convert named mode into integers
has 'typemap' => (
  is => 'rw',
  isa => 'HashRef',
  default => sub { {
    ANY => 1,
    END1 => 2,
    END2 => 3,
    HAIRPIN => 4
  }
  }
);
# counter - used internally to count the number of calls
has 'counter' => (
  is => 'rw',
  isa => 'Int',
  default => 0
);

# args - generated an array (or named hash) for arguments we will be using
method args($mode = 'array') {
  if ($mode =~ /^n/i) {
    my $args;
    map { $args->{$_} = $self->$_ } qw<type mv dv maxLoop dntp dna_conc temp>;
    return $args;
  }
  elsif ($mode =~ /^c/i) {
    my @cmd;
    push @cmd, '-a', $self->type;
    push @cmd, '-mv', $self->mv;
    push @cmd, '-dv', $self->dv;
    push @cmd, '-maxloop', $self->maxLoop;
    push @cmd, '-n', $self->dntp;
    push @cmd, '-d', $self->dna_conc;
    push @cmd, '-t', $self->temp;
    return join(' ', @cmd);
  }
  # the default will generate an array with the same order 
  # that is expected by the wrapped C function 'align'
  my @out;
  if (not defined $self->typemap->{uc($self->type)}) {
    warn (join(' ', "Alignment mode", $self->type, "is undefined, falling back to ANY"));
    $self->type('ANY');
  }
  push @out, $self->typemap->{uc($self->type)}; # default is ANY
  map { push @out, $self->$_ } qw<mv dv maxLoop dntp dna_conc temp>;
  return @out;
};

# tm - convenience method for just one single input sequence
method tm($seq) {
  $self->counter($self->counter+1);
  $seq = cleanseq($seq);
  my @args = $self->args;
  $args[3] = 0;
  my @values = align($seq, revcomp($seq), @args); # this calls the C code routine
  return $values[0];
};

# ntthal - run the actual alignment
method ntthal($seq1, $seq2 = undef) {
  if (not defined $seq1) {
    warn "No sequence given to align!";
    return 0 if $self->only;
    return { error => 'no sequence given' };
  }
  $seq1 = cleanseq($seq1);
  my @args = $self->args;
  if (length($seq1) < 10) {
    warn "The first sequence is less than 10 bases long ($seq1)";
    return 0 if $self->only;
    return { error => 'short input' };
  }
  if (not defined $seq2) {
    if ($self->type ne 'HAIRPIN') {
      $args[3] = 0; # we know there will be no gaps in the best alignment
      $seq2 = revcomp($seq1);
    } else {
      $seq2 = $seq1;
    }
  }
  else {
    $seq2 = cleanseq($seq2);
  }
  $self->counter($self->counter+1);
  my @values = align($seq1, $seq2, @args); # this calls the C code routine
  
  push @values, 'NA' while (scalar @values < 4);
  $values[0] = 0 if $values[0] < 0;
  return $values[0] if $self->only =~ /^t/i;
  return $values[3] if $self->only =~ /^[dg]/i;
  
  my $result = \%{ $self->args('named') };
  $result->{s1} = $seq1;
  $result->{s2} = $seq2;
  map { $result->{$_} = shift @values } qw<t dH dS dG duplex0 duplex1 duplex2 duplex3>;
  return $result;
};

# cleanseq - remove whitespace and convert non ACGT characters to N
func cleanseq($seq) {
  chomp($seq);
  $seq =~ s/\s+//gsm;
  $seq =~ s/\d+//gsm;
  $seq = uc($seq);
  $seq =~ s/[^ACGT]/N/gsm;
  return $seq;
};

# revcomp - reverse complement the sequence
func revcomp($in) {
  $in =~ tr/ATGC/TACG/;
  return scalar reverse $in;
};





=head1 AUTHOR

Frank Reinecke, C<< <frank.reinecke at qiagen.com> >>

=head1 LICENSE AND COPYRIGHT

Copyright 2019 Frank Reinecke.

This program is released under the following license: AGPL3.

Please note the copyright for the wrapped C code (in the sources).

=cut

__PACKAGE__->meta->make_immutable();

__DATA__
__C__

/* Copyright (c) 1996 - 2018
 Whitehead Institute for Biomedical Research, Steve Rozen, Andreas Untergasser
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
 All rights reserved.

     This file is part the primer3 software suite.

     This software suite is free software;
     you can redistribute is and/or modify it under the terms
     of the GNU General Public License as published by the Free
     Software Foundation; either version 2 of the License, or (at
     your option) any later version.

     This software is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this file (file gpl-2.0.txt in the source
     distribution); if not, write to the Free Software
     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#ifndef _THAL_H
#define _THAL_H

#include <float.h> /* ! mul ei ole float.h-d includes DBL_MAX */
#include <math.h>
#include <limits.h>


#ifndef THAL_ERROR_SCORE
#define THAL_ERROR_SCORE -_INFINITY
#endif

/* The maximum length of _one_ of the two sequences being aligned in a
   thermodynamic alignment. In other words, the length of one sequence
   must be <= THAL_MAX_ALIGN, but the other sequence can be longer.
   The rationale behind this value (60) is that this is the maxium
   reasonable length for nearest neighbor models. It is the maxium
   length at which we can restrict our model to only two states of
   melting: fully intact duplex or completely dissociated single
   strands. */
#ifndef THAL_MAX_ALIGN
#define THAL_MAX_ALIGN 60
#endif

/* The maxium length of the other sequence in a thermodynamic
   alignment. This value can be increased, though alignments against
   very long sequences will be quite slow. As of 2012-05-18, we only
   potentially see sequences longer this when checking for mispriming
   in the template ('max_template_mispriming') in libprimer3.c, which
   is really designed to find sites of ectopic primer very close (a
   few kilobases) from the location of the cadidate primer. */
#ifndef THAL_MAX_SEQ
#define THAL_MAX_SEQ   10000
#endif

/*** BEGIN CONSTANTS ***/

extern const double _INFINITY;
extern const double ABSOLUTE_ZERO;
extern const int MAX_LOOP; /* the maximum size of loop that can be calculated;
                              for larger loops formula must be implemented */
extern const int MIN_LOOP;

/*** END CONSTANTS ***/

/* BEGIN TYPEDEFs */

typedef enum thal_alignment_type {
  thal_any = 1,
  thal_end1 = 2,
  thal_end2 = 3,
  thal_hairpin = 4,
} thal_alignment_type;


/* Structure for passing arguments to THermodynamic ALignment calculation */
typedef struct thal_args {
   int type; /* one of the
              1 THAL_ANY, (by default)
              2 THAL_END1,
              3 THAL_END2,
              4 THAL_HAIRPIN */
   int maxLoop;  /* maximum size of loop to consider; longer than 30 bp are not allowed */
   double mv; /* concentration of monovalent cations */
   double dv; /* concentration of divalent cations */
   double dntp; /* concentration of dNTP-s */
   double dna_conc; /* concentration of oligonucleotides */
   double temp; /* temperature from which hairpin structures will be calculated */
   int dimer; /* if non zero, dimer structure is calculated */
} thal_args;


/* The files from the directory primer3_config loaded as strings */
typedef struct thal_parameters {
  char *dangle_dh;
  char *dangle_ds;
  char *loops_dh;
  char *loops_ds;
  char *stack_dh;
  char *stack_ds;
  char *stackmm_dh;
  char *stackmm_ds;
  char *tetraloop_dh;
  char *tetraloop_ds;
  char *triloop_dh;
  char *triloop_ds;
  char *tstack_tm_inf_ds;
  char *tstack_dh;
  char *tstack2_dh;
  char *tstack2_ds;
} thal_parameters;


/* Structure for receiving results from the thermodynamic alignment calculation */

typedef struct thal_results {
   char msg[255];
   double temp;
   int align_end_1;
   int align_end_2;
   char *sec_struct;
   double S;
   double H;
   double dG;
   char duplex0[255];
   char duplex1[255];
   char duplex2[255];
   char duplex3[255];
} thal_results;


/* 
 * THL_FAST    = 0 - score only with optimized functions (fast)
 * THL_GENERAL = 1 - use general function without debug (slow)
 * THL_DEBUG_F = 2 - debug mode with fast, print alignments on STDERR
 * THL_DEBUG   = 3 - debug mode print alignments on STDERR
 * THL_STRUCT  = 4 - calculate secondary structures as string
 */
 
typedef enum thal_mode { 
  THL_FAST    = 0,
  THL_GENERAL = 1,
  THL_DEBUG_F = 2,
  THL_DEBUG   = 3, 
  THL_STRUCT  = 4
} thal_mode;


/*** END OF TYPEDEFS ***/

void set_thal_default_args(thal_args *a);
void set_thal_oligo_default_args(thal_args *a);

/* Read the thermodynamic values (parameters) from the parameter files
   in the directory specified by 'path'.  Return 0 on success and -1
   on error. The thermodynamic values are stored in multiple static
   variables. */
/* Here is an example of how this function is used in 
   primer3_boulder_main.c: */
#if 0
  if (get_thermodynamic_values(thermodynamic_params_path, &o)) {
    fprintf(stderr, "%s\n", o.msg);
    exit(-1);
  }
#endif
int  thal_set_null_parameters(thal_parameters *a);

int  thal_load_parameters(const char *path, thal_parameters *a, thal_results* o);

int  thal_free_parameters(thal_parameters *a);

int  get_thermodynamic_values(const thal_parameters *tp, thal_results *o);

void destroy_thal_structures();

/* Central method for finding the best alignment.  On error, o->temp
   is set to THAL_ERROR_SCORE and a message is put in o->msg.  The
   error might be caused by ENOMEM. To determine this it is necessary
   to check errno.
*/

void thal(const unsigned char *oligo1, 
          const unsigned char *oligo2, 
          const thal_args* a,
          const thal_mode mode, 
          thal_results* o);

#endif
/*
Copyright (c) 2018
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky.
All rights reserved.

    This file is part of the primer3 suite and libraries.

    The primer3 suite and libraries are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef THERMODYNAMIC_PARAMETERS_H
#define THERMODYNAMIC_PARAMETERS_H 1

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * Copy the default thermodynamic parameter strings to *a 
 */
int set_default_thal_parameters(thal_parameters *a);

#ifdef __cplusplus
}
#endif

#endif

/*
 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2009,2010,
               2011,2012
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
 All rights reserved.

       This file is part of primer3 software suite.

       This software suite is is free software;
       you can redistribute it and/or modify it under the terms
       of the GNU General Public License as published by the Free
       Software Foundation; either version 2 of the License, or (at
       your option) any later version.

       This software is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this software (file gpl-2.0.txt in the source
       distribution); if not, write to the Free Software
       Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

#if defined(__sun)
#include <ieeefp.h>
#endif


/*#define DEBUG*/
#ifndef MIN_HRPN_LOOP
#define MIN_HRPN_LOOP 3 /*  minimum size of hairpin loop */
#endif

#ifndef THAL_EXIT_ON_ERROR
#define THAL_EXIT_ON_ERROR 0
#endif

/* table where bp-s enthalpies, that retrieve to the most stable Tm, are saved */
#ifdef EnthalpyDPT
# undef EnthalpyDPT
#endif
#define EnthalpyDPT(i, j) enthalpyDPT[(j) + ((i-1)*len3) - (1)]

/* table where bp-s entropies, that retrieve to the most stable Tm, are saved */
#ifdef EntropyDPT
# undef EntropyDPT
#endif
#define EntropyDPT(i, j) entropyDPT[(j) + ((i-1)*len3) - (1)]

/* entropies of most stable hairpin terminal bp */
#ifndef SEND5
# define SEND5(i) send5[i]
#endif

/* enthalpies of most stable hairpin terminal bp */
#ifndef HEND5
# define HEND5(i) hend5[i]
#endif

#define CHECK_ERROR(COND,MSG) if (COND) { strcpy(o->msg, MSG); errno = 0; longjmp(_jmp_buf, 1); }
#define THAL_OOM_ERROR { strcpy(o->msg, "Out of memory"); errno = ENOMEM; longjmp(_jmp_buf, 1); }
#define THAL_IO_ERROR(f) { sprintf(o->msg, "Unable to open file %s", f); longjmp(_jmp_buf, 1); }

#define bpIndx(a, b) BPI[a][b] /* for traceing matrix BPI */
#define atPenaltyS(a, b) atpS[a][b]
#define atPenaltyH(a, b) atpH[a][b]

#define STR(X) #X
#define LONG_SEQ_ERR_STR(MAX_LEN) "Target sequence length > maximum allowed (" STR(MAX_LEN) ") in thermodynamic alignment"
#define XSTR(X) STR(X)

#define SMALL_NON_ZERO 0.000001
#define DBL_EQ(X,Y) (((X) - (Y)) < (SMALL_NON_ZERO) ? (1) : (2)) /* 1 when numbers are equal */

#ifdef INTEGER
# define isFinite(x) (x < _INFINITY / 2)
#else
# define isFinite(x) isfinite(x)
#endif

#define isPositive(x) ((x) > 0 ? (1) : (0))

/*** BEGIN CONSTANTS ***/
# ifdef INTEGER
const double _INFINITY = 999999.0;
# else
# ifdef INFINITY
const double _INFINITY = INFINITY;
# else
const double _INFINITY = 1.0 / 0.0;
# endif
# endif

static const double R = 1.9872; /* cal/Kmol */
static const double ILAS = (-300 / 310.15); /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
static const double ILAH = 0.0; /* Internal Loop EntHalpy Asymmetry correction */
static const double AT_H = 2200.0; /* AT penalty */
static const double AT_S = 6.9; /* AT penalty */
static const double MinEntropyCutoff = -2500.0; /* to filter out non-existing entropies */
static const double MinEntropy = -3224.0; /* initiation */
static const double G2 = 0.0; /* structures w higher G are considered to be unstabile */
const double ABSOLUTE_ZERO = 273.15;
const double TEMP_KELVIN = 310.15;
const int MAX_LOOP = 30; /* the maximum size of loop that can be calculated; for larger loops formula must be implemented */
const int MIN_LOOP = 0;
//static const char BASES[5] = {'A', 'C', 'G', 'T', 'N'}; /* bases to be considered - N is every symbol that is not A, G, C,$
//                                                  */
//static const char BASE_PAIRS[4][4] = {"A-T", "C-G", "G-C", "T-A" }; /* allowed basepairs */
/* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
static const int BPI[5][5] =  {
     {0, 0, 0, 1, 0}, /* A, C, G, T, N; */
     {0, 0, 1, 0, 0},
     {0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0},
     {0, 0, 0, 0, 0}};

/*** END OF CONSTANTS ***/

/*** BEGIN STRUCTs ***/

struct triloop {
  char loop[5];
  double value; };

struct tetraloop {
  char loop[6];
  double value; };

struct tracer /* structure for tracebacku - unimolecular str */ {
  int i;
  int j;
  int mtrx; /* [0 1] EntropyDPT/EnthalpyDPT*/
  struct tracer* next;
};

/*** END STRUCTs ***/

static int length_unsig_char(const unsigned char * str); /* returns length of unsigned char; to avoid warnings while compiling */

static unsigned char str2int(char c); /* converts DNA sequence to int; 0-A, 1-C, 2-G, 3-T, 4-whatever */

static double saltCorrectS (double mv, double dv, double dntp); /* part of calculating salt correction
                                                                   for Tm by SantaLucia et al */
static char* readParamFile(const char* dirname, const char* fname, thal_results* o); /* file of thermodynamic params */

/* get thermodynamic tables */
static double readDouble(char **str, thal_results* o);

static void readLoop(char **str, double *v1, double *v2, double *v3, thal_results *o);

static int readTLoop(char **str, char *s, double *v, int triloop, thal_results *o);

static void getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o);

/*static void verifyStackTable(double stack[5][5][5][5], char* type);*/ /* just for debugging; the method is turned off by default */

static void getStackint2(double stackEntropiesint2[5][5][5][5], double stackint2Enthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o);

static void getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5],
                      double dangleEnthalpies5[5][5][5], const thal_parameters *tp, thal_results* o);

static void getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o);

static void getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o);

static void getTriloop(struct triloop**, struct triloop**, int* num, const thal_parameters *tp, thal_results* o);

static void getTetraloop(struct tetraloop**, struct tetraloop**, int* num, const thal_parameters *tp, thal_results* o);

static void getLoop(double hairpinLoopEnntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropiess[30],
             double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], const thal_parameters *tp, thal_results* o);

static void tableStartATS(double atp_value, double atp[5][5]); /* creates table of entropy values for nucleotides
                                                                  to which AT-penlty must be applied */

static void tableStartATH(double atp_value, double atp[5][5]);

static int comp3loop(const void*, const void*); /* checks if sequnece consists of specific triloop */

static int comp4loop(const void*, const void*); /* checks if sequnece consists of specific tetraloop */

static void initMatrix(); /* initiates thermodynamic parameter tables of entropy and enthalpy for dimer */

static void initMatrix2(); /* initiates thermodynamic parameter tables of entropy and enthalpy for monomer */

static void fillMatrix(int maxLoop, thal_results* o); /* calc-s thermod values into dynamic progr table (dimer) */

static void fillMatrix2(int maxLoop, thal_results* o); /* calc-s thermod values into dynamic progr table (monomer) */

static void maxTM(int i, int j); /* finds max Tm while filling the dyn progr table using stacking S and stacking H (dimer) */

static void maxTM2(int i, int j); /* finds max Tm while filling the dyn progr table using stacking S and stacking H (monomer) */

/* calculates bulges and internal loops for dimer structures */
static void calc_bulge_internal(int ii, int jj, int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop);

/* calculates bulges and internal loops for monomer structures */
static void calc_bulge_internal2(int ii, int jj, int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop);

/* carries out Bulge and Internal loop and stack calculations to hairpin */
static void CBI(int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop);

/* finds monomer structure that has maximum Tm */
static void calc_hairpin(int i, int j, double* EntropyEnthalpy, int traceback);

static double Ss(int i, int j, int k); /* returns stack entropy */
static double Hs(int i, int j, int k); /* returns stack enthalpy */

/* calculate terminal entropy S and terminal enthalpy H starting reading from 5'end (Left hand/3' end - Right end) */
static void LSH(int i, int j, double* EntropyEnthalpy);
static void RSH(int i, int j, double* EntropyEnthalpy);

static void reverse(unsigned char *s);

static int max5(double, double, double, double, double);

/* Is sequence symmetrical */
static int symmetry_thermo(const unsigned char* seq);

/* traceback for dimers */
static void traceback(int i, int j, double RT, int* ps1, int* ps2, int maxLoop, thal_results* o);

/* traceback for hairpins */
static void tracebacku(int*, int, thal_results*);

/* prints ascii output of dimer structure */
char *drawDimer(int*, int*, double, double, double, const thal_mode mode, double, thal_results*);

/* prints ascii output of hairpin structure */
char *drawHairpin(int*, double, double, const thal_mode mode, double, thal_results*);

static void save_append_string(char** ret, int *space, thal_results *o, const char *str);

static void save_append_char(char** ret, int *space, thal_results *o, const char str);

static int equal(double a, double b);

static void strcatc(char*, char);

static void push(struct tracer**, int, int, int, thal_results*); /* to add elements to struct */

/* terminal bp for monomer structure */
static void calc_terminal_bp(double temp);

/* executed in calc_terminal_bp; to find structure that corresponds to max Tm for terminal bp */
static double END5_1(int,int); /* END5_1(X,1/2) - 1=Enthalpy, 2=Entropy*/
static double END5_2(int,int);
static double END5_3(int,int);
static double END5_4(int,int);

static double Hd5(int,int); /* returns thermodynamic value (H) for 5' dangling end */
static double Hd3(int,int); /* returns thermodynamic value (H) for 3' dangling end */
static double Sd5(int,int); /* returns thermodynamic value (S) for 5' dangling end */
static double Sd3(int,int); /* returns thermodynamic value (S) for 3' dangling end */
static double Ststack(int,int); /* returns entropy value for terminal stack */
static double Htstack(int,int); /* returns enthalpy value for terminal stack */

/* memory stuff */
static void* safe_calloc(size_t, size_t, thal_results* o);
static void* safe_malloc(size_t, thal_results* o);
static void* safe_realloc(void*, size_t, thal_results* o);
static double* safe_recalloc(double* ptr, int m, int n, thal_results* o);

static int numTriloops; /* hairpin triloop penalties */
static int numTetraloops; /* hairpin tetraloop penalties */
static double atpS[5][5]; /* AT penalty */
static double atpH[5][5]; /* AT penalty */
static double *send5, *hend5; /* calc 5'  */
/* w/o init not constant anymore, cause for unimolecular and bimolecular foldings there are different values */
static double dplx_init_H; /* initiation enthalpy; for duplex 200, for unimolecular structure 0 */
static double dplx_init_S; /* initiation entropy; for duplex -5.7, for unimoleculat structure 0 */
static double saltCorrection; /* value calculated by saltCorrectS, includes correction for monovalent and divalent cations */
static double RC; /* universal gas constant multiplied w DNA conc - for melting temperature */
static double SHleft; /* var that helps to find str w highest melting temperature */
static int bestI, bestJ; /* starting position of most stable str */
static double* enthalpyDPT; /* matrix for values of enthalpy */
static double* entropyDPT; /* matrix for values of entropy */
static unsigned char *oligo1, *oligo2; /* inserted oligo sequenced */
static unsigned char *numSeq1, *numSeq2; /* same as oligo1 and oligo2 but converted to numbers */
static int len1, len2, len3; /* length of sequense 1 and 2 *//* 17.02.2009 int temponly;*/ /* print only temperature of the predicted structure */
static double dangleEntropies3[5][5][5]; /* thermodynamic paramteres for 3' dangling ends */
static double dangleEnthalpies3[5][5][5]; /* ther params for 3' dangling ends */
static double dangleEntropies5[5][5][5];  /* ther params for 5' dangling ends */
static double dangleEnthalpies5[5][5][5]; /* ther params for 5' dangling ends */
static double stackEntropies[5][5][5][5]; /* ther params for perfect match pairs */
static double stackEnthalpies[5][5][5][5]; /* ther params for perfect match pairs */
static double stackint2Entropies[5][5][5][5]; /*ther params for perfect match and internal mm */
static double stackint2Enthalpies[5][5][5][5]; /* ther params for perfect match and internal mm*/
static double interiorLoopEntropies[30]; /* interior loop params according to length of the loop */
static double bulgeLoopEntropies[30]; /* bulge loop params according to length of the loop */
static double hairpinLoopEntropies[30]; /* hairpin loop params accordint to length of the loop */
static double interiorLoopEnthalpies[30]; /* same as interiorLoopEntropies but values of entropy */
static double bulgeLoopEnthalpies[30]; /* same as bulgeLoopEntropies but values of entropy */
static double hairpinLoopEnthalpies[30]; /* same as hairpinLoopEntropies but values of entropy */
static double tstackEntropies[5][5][5][5]; /* ther params for terminal mismatches */
static double tstackEnthalpies[5][5][5][5]; /* ther params for terminal mismatches */
static double tstack2Entropies[5][5][5][5]; /* ther params for internal terminal mismatches */
static double tstack2Enthalpies[5][5][5][5]; /* ther params for internal terminal mismatches */
static struct triloop* triloopEntropies = NULL; /* ther penalties for given triloop seq-s */
static struct triloop* triloopEnthalpies = NULL; /* ther penalties for given triloop seq-s */
static struct tetraloop* tetraloopEntropies = NULL; /* ther penalties for given tetraloop seq-s */
static struct tetraloop* tetraloopEnthalpies = NULL; /* ther penalties for given tetraloop seq-s */
static jmp_buf _jmp_buf;

/* Initialize the thermodynamic values (parameters) */
int  thal_set_null_parameters(thal_parameters *a) {
  a->dangle_dh = NULL;
  a->dangle_ds = NULL;
  a->loops_dh = NULL;
  a->loops_ds = NULL;
  a->stack_dh = NULL;
  a->stack_ds = NULL;
  a->stackmm_dh = NULL;
  a->stackmm_ds = NULL;
  a->tetraloop_dh = NULL;
  a->tetraloop_ds = NULL;
  a->triloop_dh = NULL;
  a->triloop_ds = NULL;
  a->tstack_tm_inf_ds = NULL;
  a->tstack_dh = NULL;
  a->tstack2_dh = NULL;
  a->tstack2_ds = NULL;
  return 0;
}

/* Free the thermodynamic values (parameters) */
int  thal_free_parameters(thal_parameters *a) {
  if (NULL != a->dangle_dh) {
    free(a->dangle_dh);
    a->dangle_dh = NULL;
  }
  if (NULL != a->dangle_ds) {
    free(a->dangle_ds);
    a->dangle_ds = NULL;
  }
  if (NULL != a->loops_dh) {
    free(a->loops_dh);
    a->loops_dh = NULL;
  }
  if (NULL != a->loops_ds) {
    free(a->loops_ds);
    a->loops_ds = NULL;
  }
  if (NULL != a->stack_dh) {
    free(a->stack_dh);
    a->stack_dh = NULL;
  }
  if (NULL != a->stack_ds) {
    free(a->stack_ds);
    a->stack_ds = NULL;
  }
  if (NULL != a->stackmm_dh) {
    free(a->stackmm_dh);
    a->stackmm_dh = NULL;
  }
  if (NULL != a->stackmm_ds) {
    free(a->stackmm_ds);
    a->stackmm_ds = NULL;
  }
  if (NULL != a->tetraloop_dh) {
    free(a->tetraloop_dh);
    a->tetraloop_dh = NULL;
  }
  if (NULL != a->tetraloop_ds) {
    free(a->tetraloop_ds);
    a->tetraloop_ds = NULL;
  }
  if (NULL != a->triloop_dh) {
    free(a->triloop_dh);
    a->triloop_dh = NULL;
  }
  if (NULL != a->triloop_ds) {
    free(a->triloop_ds);
    a->triloop_ds = NULL;
  }
  if (NULL != a->tstack_tm_inf_ds) {
    free(a->tstack_tm_inf_ds);
    a->tstack_tm_inf_ds = NULL;
  }
  if (NULL != a->tstack_dh) {
    free(a->tstack_dh);
    a->tstack_dh = NULL;
  }
  if (NULL != a->tstack2_dh) {
    free(a->tstack2_dh);
    a->tstack2_dh = NULL;
  }
  if (NULL != a->tstack2_ds) {
    free(a->tstack2_ds);
    a->tstack2_ds = NULL;
  }
  return 0;
}

/* Read the thermodynamic values (parameters) from the parameter files
   in the directory specified by 'path'.  Return 0 on success and -1
   on error. The thermodynamic values are stored in multiple static
   variables. */
int
get_thermodynamic_values(const thal_parameters *tp, thal_results *o)
{
  if (setjmp(_jmp_buf) != 0) {
     return -1;
  }
  getStack(stackEntropies, stackEnthalpies, tp, o);
  /* verifyStackTable(stackEntropies, "entropy");
     verifyStackTable(stackEnthalpies, "enthalpy"); */ /* this is for code debugging */
  getStackint2(stackint2Entropies, stackint2Enthalpies, tp, o);
  getDangle(dangleEntropies3, dangleEnthalpies3, dangleEntropies5, dangleEnthalpies5, tp, o);
  getLoop(hairpinLoopEntropies, interiorLoopEntropies, bulgeLoopEntropies, hairpinLoopEnthalpies,
          interiorLoopEnthalpies, bulgeLoopEnthalpies, tp, o);
  getTstack(tstackEntropies, tstackEnthalpies, tp, o);
  getTstack2(tstack2Entropies, tstack2Enthalpies, tp, o);
  getTriloop(&triloopEntropies, &triloopEnthalpies, &numTriloops, tp, o);
  getTetraloop(&tetraloopEntropies, &tetraloopEnthalpies, &numTetraloops, tp, o);
  /* getting the AT-penalties */
  tableStartATS(AT_S, atpS);
  tableStartATH(AT_H, atpH);

  return 0;
}

void
destroy_thal_structures()
{
  if (triloopEntropies != NULL) {
    free(triloopEntropies);
    triloopEntropies = NULL;
  }
  if (triloopEnthalpies != NULL) {
    free(triloopEnthalpies);
    triloopEnthalpies = NULL;
  }
  if (tetraloopEntropies != NULL) {
    free(tetraloopEntropies);
    tetraloopEntropies = NULL;
  }
  if (tetraloopEnthalpies != NULL) {
    free(tetraloopEnthalpies);
    tetraloopEnthalpies = NULL;
  }
}

/* central method: execute all sub-methods for calculating secondary
   structure for dimer or for monomer */
void
thal(const unsigned char *oligo_f,
     const unsigned char *oligo_r,
     const thal_args *a,
     const thal_mode mode,
     thal_results *o)
{
   double* SH;
   int i, j;
   int len_f, len_r;
   int k;
   int *bp;
   unsigned char *oligo2_rev = NULL;
   double mh, ms;
   double G1, bestG;

   send5 = hend5 = NULL;
   enthalpyDPT = entropyDPT = NULL;
   numSeq1 = numSeq2 = NULL;
   oligo1 = oligo2 = NULL;
   strcpy(o->msg, "");
   o->temp = THAL_ERROR_SCORE;
   errno = 0;

   if (setjmp(_jmp_buf) != 0) {
     o->temp = THAL_ERROR_SCORE;
     return;  /* If we get here, that means we returned via a
                 longjmp.  In this case errno might be ENOMEM,
                 but not necessarily. */
   }

   CHECK_ERROR(NULL == oligo_f, "NULL first sequence");
   CHECK_ERROR(NULL == oligo_r, "NULL second sequence");
   len_f = length_unsig_char(oligo_f);
   len_r = length_unsig_char(oligo_r);

   /*CHECK_ERROR(1==len_f, "Length 1 first sequence");
   CHECK_ERROR(1==len_r, "Length 1 second sequence"); */
   /* The following error messages will be seen by end users and will
      not be easy to understand. */
   CHECK_ERROR((len_f > THAL_MAX_ALIGN) && (len_r > THAL_MAX_ALIGN),
               "Both sequences longer than " XSTR(THAL_MAX_ALIGN)
               " for thermodynamic alignment");
   CHECK_ERROR((len_f > THAL_MAX_SEQ),
               LONG_SEQ_ERR_STR(THAL_MAX_SEQ) " (1)");
   CHECK_ERROR((len_r > THAL_MAX_SEQ),
               LONG_SEQ_ERR_STR(THAL_MAX_SEQ) " (2)");

   CHECK_ERROR(NULL == a,  "NULL 'in' pointer");
   if (NULL == o) return; /* Leave it to the caller to crash */
   CHECK_ERROR(a->type != thal_any
               && a->type != thal_end1
               && a->type != thal_end2
               && a->type != thal_hairpin,
               "Illegal type");
   o->align_end_1 = -1;
   o->align_end_2 = -1;
   if (oligo_f && '\0' == *oligo_f) {
      strcpy(o->msg, "Empty first sequence");
      o->temp = 0.0;
      return;
   }
   if (oligo_r && '\0' == *oligo_r) {
      strcpy(o->msg, "Empty second sequence");
      o->temp = 0.0;
      return;
   }
   if (0 == len_f) {
      o->temp = 0.0;
      return;
   }
   if (0 == len_r) {
      o->temp = 0.0;
      return;
   }
   if(a->type!=3) {
      oligo1 = (unsigned char*) safe_malloc((len_f + 1) * sizeof(unsigned char), o);
      oligo2 = (unsigned char*) safe_malloc((len_r + 1) * sizeof(unsigned char), o);
      strcpy((char*)oligo1,(const char*)oligo_f);
      strcpy((char*)oligo2,(const char*)oligo_r);
   } else  {
      oligo1 = (unsigned char*) safe_malloc((len_r + 1) * sizeof(unsigned char), o);
      oligo2 = (unsigned char*) safe_malloc((len_f + 1) * sizeof(unsigned char), o);
      strcpy((char*)oligo1,(const char*)oligo_r);
      strcpy((char*)oligo2,(const char*)oligo_f);
   }
   /*** INIT values for unimolecular and bimolecular structures ***/
   if (a->type==4) { /* unimolecular folding */
      len2 = length_unsig_char(oligo2);
      len3 = len2 -1;
      dplx_init_H = 0.0;
      dplx_init_S = -0.00000000001;
      RC=0;
   } else if(a->type!=4) {
      /* hybridization of two oligos */
      dplx_init_H = 200;
      dplx_init_S = -5.7;
      if(symmetry_thermo(oligo1) && symmetry_thermo(oligo2)) {
         RC = R  * log(a->dna_conc/1000000000.0);
      } else {
         RC = R  * log(a->dna_conc/4000000000.0);
      }
      if(a->type!=3) {
         oligo2_rev = (unsigned char*) safe_malloc((length_unsig_char(oligo_r) + 1) * sizeof(unsigned char), o);
         strcpy((char*)oligo2_rev,(const char*)oligo_r);
      } else {
         oligo2_rev = (unsigned char*) safe_malloc((length_unsig_char(oligo_f) + 1) * sizeof(unsigned char), o);
         strcpy((char*)oligo2_rev,(const char*)oligo_f);
      }
      reverse(oligo2_rev); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
      free(oligo2);
      oligo2=NULL;
      oligo2=&oligo2_rev[0];
   } else {
      strcpy(o->msg, "Wrong alignment type!");
      o->temp = THAL_ERROR_SCORE;
      errno=0;
#ifdef DEBUG
      fprintf(stderr, o->msg);
#endif
      return;
   }
   len1 = length_unsig_char(oligo1);
   len2 = length_unsig_char(oligo2);
   /* convert nucleotides to numbers */
   numSeq1 = (unsigned char*) safe_realloc(numSeq1, len1 + 2, o);
   numSeq2 = (unsigned char*) safe_realloc(numSeq2, len2 + 2, o);

   /*** Calc part of the salt correction ***/
   saltCorrection=saltCorrectS(a->mv,a->dv,a->dntp); /* salt correction for entropy, must be multiplied with N, which is
                                                   the total number of phosphates in the duplex divided by 2; 8bp dplx N=7 */

   if(a->type == 4){ /* monomer */
      /* terminal basepairs */
      send5 = (double*) safe_realloc(send5, (len1 + 1) * sizeof(double), o);
      hend5 = (double*) safe_realloc(hend5, (len1 + 1) * sizeof(double), o);
   }
   for(i = 0; i < len1; i++) oligo1[i] = toupper(oligo1[i]);
   for(i = 0; i < len2; i++) oligo2[i] = toupper(oligo2[i]);
   for(i = 1; i <= len1; ++i) numSeq1[i] = str2int(oligo1[i - 1]);
   for(i = 1; i <= len2; ++i) numSeq2[i] = str2int(oligo2[i - 1]);
   numSeq1[0] = numSeq1[len1 + 1] = numSeq2[0] = numSeq2[len2 + 1] = 4; /* mark as N-s */
   if (a->type==4) { /* calculate structure of monomer */
      enthalpyDPT = safe_recalloc(enthalpyDPT, len1, len2, o);
      entropyDPT = safe_recalloc(entropyDPT, len1, len2, o);
      initMatrix2();
      fillMatrix2(a->maxLoop, o);
      calc_terminal_bp(a->temp);
      mh = HEND5(len1);
      ms = SEND5(len1);
      o->align_end_1 = (int) mh;
      o->align_end_2 = (int) ms;
      bp = (int*) safe_calloc(len1, sizeof(int), o);
      for (k = 0; k < len1; ++k) bp[k] = 0;
      if(isFinite(mh)) {
        tracebacku(bp, a->maxLoop, o);
        /* traceback for unimolecular structure */
        o->sec_struct=drawHairpin(bp, mh, ms, mode,a->temp, o); /* if mode=THL_FAST or THL_DEBUG_F then return after printing basic therm data */
      } else if((mode != THL_FAST) && (mode != THL_DEBUG_F) && (mode != THL_STRUCT)) {
        fputs("No secondary structure could be calculated\n",stderr);
      }

      if(o->temp==-_INFINITY && (!strcmp(o->msg, ""))) o->temp=0.0;
      free(bp);
      free(enthalpyDPT);
      free(entropyDPT);
      free(numSeq1);
      free(numSeq2);
      free(send5);
      free(hend5);
      free(oligo1);
      free(oligo2);
      return;
   } else if(a->type!=4) { /* Hybridization of two moleculs */
      len3 = len2;
      enthalpyDPT = safe_recalloc(enthalpyDPT, len1, len2, o); /* dyn. programming table for dS and dH */
      entropyDPT = safe_recalloc(entropyDPT, len1, len2, o); /* enthalpyDPT is 3D array represented as 1D array */
      initMatrix();
      fillMatrix(a->maxLoop, o);
      SH = (double*) safe_malloc(2 * sizeof(double), o);
      /* calculate terminal basepairs */
      bestI = bestJ = 0;
      G1 = bestG = _INFINITY;
      if(a->type==1)
        for (i = 1; i <= len1; i++) {
           for (j = 1; j <= len2; j++) {
              RSH(i, j, SH);
              SH[0] = SH[0]+SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0 */
              SH[1] = SH[1]+SMALL_NON_ZERO;
              G1 = (EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) - TEMP_KELVIN*(EntropyDPT(i, j) + SH[0] + dplx_init_S);
              if(G1<bestG){
                 bestG = G1;
                 bestI = i;
                 bestJ = j;
              }
           }
        }
      int *ps1, *ps2;
      ps1 = (int*) safe_calloc(len1, sizeof(int), o);
      ps2 = (int*) safe_calloc(len2, sizeof(int), o);
      for (i = 0; i < len1; ++i)
        ps1[i] = 0;
      for (j = 0; j < len2; ++j)
        ps2[j] = 0;
      if(a->type == 2 || a->type == 3)        {
         /* THAL_END1 */
         bestI = bestJ = 0;
         bestI = len1;
         i = len1;
         G1 = bestG = _INFINITY;
         for (j = 1; j <= len2; ++j) {
            RSH(i, j, SH);
            SH[0] = SH[0]+SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0,
                                             that compiler could understand that SH is changed in this cycle */
            SH[1] = SH[1]+SMALL_NON_ZERO;
            G1 = (EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) - TEMP_KELVIN*(EntropyDPT(i, j) + SH[0] + dplx_init_S);
                if(G1<bestG){
                   bestG = G1;
                   bestJ = j;
            }
         }
      }
      if (!isFinite(bestG)) bestI = bestJ = 1;
      double dH, dS;
      RSH(bestI, bestJ, SH);
      dH = EnthalpyDPT(bestI, bestJ)+ SH[1] + dplx_init_H;
      dS = (EntropyDPT(bestI, bestJ) + SH[0] + dplx_init_S);
      /* tracebacking */
      for (i = 0; i < len1; ++i)
        ps1[i] = 0;
      for (j = 0; j < len2; ++j)
        ps2[j] = 0;
      if(isFinite(EnthalpyDPT(bestI, bestJ))){
         traceback(bestI, bestJ, RC, ps1, ps2, a->maxLoop, o);
         o->sec_struct=drawDimer(ps1, ps2, SHleft, dH, dS, mode, a->temp, o);
         o->align_end_1=bestI;
         o->align_end_2=bestJ;
      } else  {
         o->temp = 0.0;
         /* fputs("No secondary structure could be calculated\n",stderr); */
      }
      free(ps1);
      free(ps2);
      free(SH);
      free(oligo2_rev);
      free(enthalpyDPT);
      free(entropyDPT);
      free(numSeq1);
      free(numSeq2);
      free(oligo1);
      return;
   }
   return;
}
/*** END thal() ***/

/* Set default args */
void
set_thal_default_args(thal_args *a)
{
   memset(a, 0, sizeof(*a));
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.8; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = TEMP_KELVIN; /* Kelvin */
   a->dimer = 1; /* by default dimer structure is calculated */
}

/* Set default args for oligo */
void
set_thal_oligo_default_args(thal_args *a)
{
   memset(a, 0, sizeof(*a));
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.0; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = TEMP_KELVIN; /* Kelvin */
   a->dimer = 1; /* by default dimer structure is calculated */
}


static unsigned char
str2int(char c)
{
   switch (c) {
    case 'A': case '0':
      return 0;
    case 'C': case '1':
      return 1;
    case 'G': case '2':
      return 2;
    case 'T': case '3':
      return 3;
   }
   return 4;
}

/* memory stuff */

static double*
safe_recalloc(double* ptr, int m, int n, thal_results* o)
{
   return (double*) safe_realloc(ptr, m * n * sizeof(double), o);
}

static void*
safe_calloc(size_t m, size_t n, thal_results *o)
{
   void* ptr;
   if (!(ptr = calloc(m, n))) {
#ifdef DEBUG
      fputs("Error in calloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
   }
   return ptr;
}

static void*
safe_malloc(size_t n, thal_results *o)
{
   void* ptr;
   if (!(ptr = malloc(n))) {
#ifdef DEBUG
      fputs("Error in malloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
   }
   return ptr;
}

static void*
safe_realloc(void* ptr, size_t n, thal_results *o)
{
   ptr = realloc(ptr, n);
   if (ptr == NULL) {
#ifdef DEBUG
      fputs("Error in realloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
   }
   return ptr;
}

static int
max5(double a, double b, double c, double d, double e)
{
   if(a > b && a > c && a > d && a > e) return 1;
   else if(b > c && b > d && b > e) return 2;
   else if(c > d && c > e) return 3;
   else if(d > e) return 4;
   else return 5;
}

static void
push(struct tracer** stack, int i, int j, int mtrx, thal_results* o)
{
   struct tracer* new_top;
   new_top = (struct tracer*) safe_malloc(sizeof(struct tracer), o);
   new_top->i = i;
   new_top->j = j;
   new_top->mtrx = mtrx;
   new_top->next = *stack;
   *stack = new_top;
}

static void
reverse(unsigned char *s)
{
   int i,j;
   char c;
   for (i = 0, j = length_unsig_char(s)-1; i < j; i++, j--) {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}

#define INIT_BUF_SIZE 1024

static char*
readParamFile(const char* dirname, const char* fname, thal_results* o)
{
  FILE* file;
  char* ret = NULL;
  char* paramdir = NULL;
  paramdir = (char*) safe_malloc(strlen(dirname) + strlen(fname) + 2, o);
  strcpy(paramdir, dirname);
#ifdef OS_WIN
  if (paramdir[strlen(paramdir) - 1] != '\\') {
    strcat(paramdir, "\\\0");
  }
#else
  if (paramdir[strlen(paramdir) - 1] != '/') {
    strcat(paramdir, "/\0");
  }
#endif
  strcat(paramdir, fname);
  if (!(file = fopen(paramdir, "r"))) {
    sprintf(o->msg, "Unable to open file %s", paramdir);
    if (paramdir != NULL) {
      free(paramdir);
      paramdir = NULL;
    }
    longjmp(_jmp_buf, 1);
    return NULL;
  }
  if (paramdir != NULL) {
    free(paramdir);
    paramdir = NULL;
  }
  char c;
  int i = 0;
  size_t ssz = INIT_BUF_SIZE;
  size_t remaining_size;
  remaining_size = ssz;
  ret = (char*) safe_malloc(ssz, o);
  while (1) {
    if (feof(file)) {
      ret[i] = '\0';
      fclose(file);
      return ret;
    }
    c = fgetc(file);
    remaining_size -= sizeof(char);
    if (remaining_size <= 0) {
      if (ssz >= INT_MAX / 2) {
        strcpy(o->msg, "Out of memory");
   free(ret);
   longjmp(_jmp_buf, 1);
   return NULL;
      } else {
        ssz += INIT_BUF_SIZE;
   remaining_size += INIT_BUF_SIZE;
      }
      ret = (char *) safe_realloc(ret, ssz, o);
    }
    ret[i] = c;
    i++;
  }
}

int
thal_load_parameters(const char *path, thal_parameters *a, thal_results* o)
{
  thal_free_parameters(a);
  if (setjmp(_jmp_buf) != 0) {
    printf("longjump\n");
    return -1;
  }
  a->dangle_dh = readParamFile(path, "dangle.dh", o);
  a->dangle_ds = readParamFile(path, "dangle.ds", o);
  a->loops_dh = readParamFile(path, "loops.dh", o);
  a->loops_ds = readParamFile(path, "loops.ds", o);
  a->stack_dh = readParamFile(path, "stack.dh", o);
  a->stack_ds = readParamFile(path, "stack.ds", o);
  a->stackmm_dh = readParamFile(path, "stackmm.dh", o);
  a->stackmm_ds = readParamFile(path, "stackmm.ds", o);
  a->tetraloop_dh = readParamFile(path, "tetraloop.dh", o);
  a->tetraloop_ds = readParamFile(path, "tetraloop.ds", o);
  a->triloop_dh = readParamFile(path, "triloop.dh", o);
  a->triloop_ds = readParamFile(path, "triloop.ds", o);
  a->tstack_tm_inf_ds = readParamFile(path, "tstack_tm_inf.ds", o);
  a->tstack_dh = readParamFile(path, "tstack.dh", o);
  a->tstack2_dh = readParamFile(path, "tstack2.dh", o);
  a->tstack2_ds = readParamFile(path, "tstack2.ds", o);
  return 0;
}

static double
saltCorrectS (double mv, double dv, double dntp)
{
   if(dv<=0) dntp=dv;
   return 0.368*((log((mv+120*(sqrt(fmax(0.0, dv-dntp))))/1000)));
}

static char*
th_read_str_line(char **str, thal_results* o)
{
  if (*str == NULL) {
    return NULL;
  }
  char *ptr = *str;
  char *ini = *str;
  while(1) {
    if ((*ptr == '\n') || (*ptr == '\0')) {
      char *ret = NULL;
      if (!(ret = malloc(sizeof(char) * (ptr - ini + 1)))) {
#ifdef DEBUG
        fputs("Error in malloc()\n", stderr);
#endif
        THAL_OOM_ERROR;
      }
      /* copy line */
      strncpy(ret, ini, (ptr - ini + 1));
      ret[ptr - ini] = '\0';

      if (*ptr == '\0') { /* End of String */
        *str = NULL;
      } else {
        ptr++;
        if (*ptr == '\0') { /* End of String */
          *str = NULL;
        } else {
          *str = ptr;
        }
      }
      if (ptr == ini) {
        if (ret != NULL) {
          free(ret);
   }
        return NULL;
      } else {
        return ret;
      }
    }
    ptr++;
  }
}

/* These functions are needed as "inf" cannot be read on Windows directly */
static double
readDouble(char **str, thal_results* o)
{
  double result;
  char *line = th_read_str_line(str, o);
  /* skip any spaces at beginning of the line */
  while (isspace(*line)) line++;
  if (!strncmp(line, "inf", 3)) {
    free(line);
    return _INFINITY;
  }
  sscanf(line, "%lf", &result);
  if (line != NULL) {
    free(line);
  }
  return result;
}

/* Reads a line containing 4 doubles, which can be specified as "inf". */
static void
readLoop(char **str, double *v1, double *v2, double *v3, thal_results *o)
{
  char *line = th_read_str_line(str, o);
  char *p = line, *q;
  /* skip first number on the line */
  while (isspace(*p)) p++;
  while (isdigit(*p)) p++;
  while (isspace(*p)) p++;
  /* read second number */
  q = p;
  while (!isspace(*q)) q++;
  *q = '\0'; q++;
  if (!strcmp(p, "inf")) *v1 = _INFINITY;
  else sscanf(p, "%lf", v1);
  while (isspace(*q)) q++;
  /* read third number */
  p = q;
  while (!isspace(*p)) p++;
  *p = '\0'; p++;
  if (!strcmp(q, "inf")) *v2 = _INFINITY;
  else sscanf(q, "%lf", v2);
  while (isspace(*p)) p++;
  /* read last number */
  q = p;
  while (!isspace(*q) && (*q != '\0')) q++;
  *q = '\0';
  if (!strcmp(p, "inf")) *v3 = _INFINITY;
  else sscanf(p, "%lf", v3);
  if (line != NULL) {
    free(line);
  }
}

/* Reads a line containing a short string and a double, used for reading a triloop or tetraloop. */
static int
readTLoop(char **str, char *s, double *v, int triloop, thal_results *o)
{
  char *line = th_read_str_line(str, o);
  if (!line) return -1;
  char *p = line, *q;
  /* skip first spaces */
  while (isspace(*p)) p++;
  /* read the string */
  q = p;
  while (isalpha(*q)) q++;
  *q = '\0'; q++;
  if (triloop) {
    strncpy(s, p, 5);   /*triloop string has 5 characters*/
    s[5] = '\0';
  } else {
    strncpy(s, p, 6);   /*tetraloop string has 6 characters*/
    s[6] = '\0';
  }
  /* skip all spaces */
  while (isspace(*q)) q++;
  p = q;
  while (!isspace(*p) && (*p != '\0')) p++;
  *p = '\0';
  if (!strcmp(q, "inf")) *v = _INFINITY;
  else sscanf(q, "%lg", v);
  if (line != NULL) {
    free(line);
  }
  return 0;
}

static void
getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o)
{
   int i, j, ii, jj;
   char *pt_ds = tp->stack_ds;
   char *pt_dh = tp->stack_dh;
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
         for (j = 0; j < 5; ++j) {
            for (jj = 0; jj < 5; ++jj) {
               if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                  stackEntropies[i][ii][j][jj] = -1.0;
                  stackEnthalpies[i][ii][j][jj] = _INFINITY;
               } else {
                  stackEntropies[i][ii][j][jj] = readDouble(&pt_ds, o);
                  stackEnthalpies[i][ii][j][jj] = readDouble(&pt_dh, o);
                  if (!isFinite(stackEntropies[i][ii][j][jj]) || !isFinite(stackEnthalpies[i][ii][j][jj])) {
                     stackEntropies[i][ii][j][jj] = -1.0;
                     stackEnthalpies[i][ii][j][jj] = _INFINITY;
                  }
               }
            }
         }
      }
   }
}

static void
getStackint2(double stackint2Entropies[5][5][5][5], double stackint2Enthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o)
{
   int i, j, ii, jj;
   char *pt_ds = tp->stackmm_ds;
   char *pt_dh = tp->stackmm_dh;
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
         for (j = 0; j < 5; ++j) {
            for (jj = 0; jj < 5; ++jj) {
               if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                  stackint2Entropies[i][ii][j][jj] = -1.0;
                  stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
               } else {
                  stackint2Entropies[i][ii][j][jj] = readDouble(&pt_ds, o);
                  stackint2Enthalpies[i][ii][j][jj] = readDouble(&pt_dh, o);
                  if (!isFinite(stackint2Entropies[i][ii][j][jj]) || !isFinite(stackint2Enthalpies[i][ii][j][jj])) {
                     stackint2Entropies[i][ii][j][jj] = -1.0;
                     stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
                  }
               }
            }
         }
      }
   }
}

/*
static void
verifyStackTable(double stack[5][5][5][5], char* type)
{

   int i, j, ii, jj;
   for (i = 0; i < 4; ++i)
     for (j = 0; j < 4; ++j)
       for (ii = 0; ii < 4; ++ii)
         for (jj = 0; jj < 4; ++jj)
           if (stack[i][j][ii][jj] != stack[jj][ii][j][i])
#ifdef DEBUG
             fprintf(stderr, "Warning: symmetrical stacks _are_ _not_ equal: %c-%c/%c-%c stack %s is %g; %c-%c/%c-%c stack %s is %g\n",
#endif
                     BASES[i], BASES[j], BASES[ii], BASES[jj], type, stack[i][j][ii][jj], BASES[jj],
                     BASES[ii], BASES[j], BASES[i], type, stack[jj][ii][j][i]);
}
*/

static void
getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5],
          double dangleEnthalpies5[5][5][5], const thal_parameters *tp, thal_results* o)
{
   int i, j, k;
   char *pt_ds = tp->dangle_ds;
   char *pt_dh = tp->dangle_dh;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
          if (i == 4 || j == 4) {
             dangleEntropies3[i][k][j] = -1.0;
             dangleEnthalpies3[i][k][j] = _INFINITY;
          } else if (k == 4) {
             dangleEntropies3[i][k][j] = -1.0;
             dangleEnthalpies3[i][k][j] = _INFINITY;
          } else {
             dangleEntropies3[i][k][j] = readDouble(&pt_ds, o);
             dangleEnthalpies3[i][k][j] = readDouble(&pt_dh, o);
             if(!isFinite(dangleEntropies3[i][k][j]) || !isFinite(dangleEnthalpies3[i][k][j])) {
                dangleEntropies3[i][k][j] = -1.0;
                dangleEnthalpies3[i][k][j] = _INFINITY;
             }
          }
       }

   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
          if (i == 4 || j == 4) {
             dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
          } else if (k == 4) {
             dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
          } else {
             dangleEntropies5[i][j][k] = readDouble(&pt_ds, o);
             dangleEnthalpies5[i][j][k] = readDouble(&pt_dh, o);
             if(!isFinite(dangleEntropies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k])) {
                dangleEntropies5[i][j][k] = -1.0;
                dangleEnthalpies5[i][j][k] = _INFINITY;
             }
          }
       }
}

static void
getLoop(double hairpinLoopEntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropies[30],
        double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30],
        const thal_parameters *tp, thal_results* o)
{
   int k;
   char *pt_ds = tp->loops_ds;
   char *pt_dh = tp->loops_dh;
   for (k = 0; k < 30; ++k) {
      readLoop(&pt_ds, &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k], o);
      readLoop(&pt_dh, &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k], o);
   }
}

static void
getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o)
{
   int i1, j1, i2, j2;
   char *pt_ds = tp->tstack_tm_inf_ds;
   char *pt_dh = tp->tstack_dh;
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
         for (j2 = 0; j2 < 5; ++j2)
           if (i1 == 4 || j1 == 4) {
              tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
              tstackEntropies[i1][i2][j1][j2] = -1.0;
           } else if (i2 == 4 || j2 == 4) {
              tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
              tstackEnthalpies[i1][i2][j1][j2] = 0.0;
           } else {
              tstackEntropies[i1][i2][j1][j2] = readDouble(&pt_ds, o);
              tstackEnthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, o);
              if (!isFinite(tstackEntropies[i1][i2][j1][j2]) || !isFinite(tstackEnthalpies[i1][i2][j1][j2])) {
                 tstackEntropies[i1][i2][j1][j2] = -1.0;
                 tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
              }
           }
}

static void
getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o)
{

   int i1, j1, i2, j2;
   char *pt_ds = tp->tstack2_ds;
   char *pt_dh = tp->tstack2_dh;
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
         for (j2 = 0; j2 < 5; ++j2)
           if (i1 == 4 || j1 == 4)  {
              tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
              tstack2Entropies[i1][i2][j1][j2] = -1.0;
           } else if (i2 == 4 || j2 == 4) {
              tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
              tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
           } else {
              tstack2Entropies[i1][i2][j1][j2] = readDouble(&pt_ds, o);
              tstack2Enthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, o);
              if (!isFinite(tstack2Entropies[i1][i2][j1][j2]) || !isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
                 tstack2Entropies[i1][i2][j1][j2] = -1.0;
                 tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
              }
           }
}

static void
getTriloop(struct triloop** triloopEntropies, struct triloop** triloopEnthalpies, int* num, const thal_parameters *tp, thal_results* o)
{
   int i, size;
   double value;
   char *pt_ds = tp->triloop_ds;
   *num = 0;
   size = 16;
   if (*triloopEntropies != NULL) {
     free(*triloopEntropies);
     *triloopEntropies = NULL;
   }
   *triloopEntropies = (struct triloop*) safe_calloc(16, sizeof(struct triloop), o);
   while (readTLoop(&pt_ds, (*triloopEntropies)[*num].loop, &value, 1, o) != -1) {
      for (i = 0; i < 5; ++i)
        (*triloopEntropies)[*num].loop[i] = str2int((*triloopEntropies)[*num].loop[i]);
      (*triloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size)        {
         size *= 2;
         *triloopEntropies = (struct triloop*) safe_realloc(*triloopEntropies, size * sizeof(struct triloop), o);
      }
   }
   *triloopEntropies = (struct triloop*) safe_realloc(*triloopEntropies, *num * sizeof(struct triloop), o);

   char *pt_dh = tp->triloop_dh;
   *num = 0;
   size = 16;

   if (*triloopEnthalpies != NULL) {
     free(*triloopEnthalpies);
     *triloopEnthalpies = NULL;
   }
   *triloopEnthalpies = (struct triloop*) safe_calloc(16, sizeof(struct triloop), o);
   while (readTLoop(&pt_dh, (*triloopEnthalpies)[*num].loop, &value, 1, o) != -1) {
      for (i = 0; i < 5; ++i)
        (*triloopEnthalpies)[*num].loop[i] = str2int((*triloopEnthalpies)[*num].loop[i]);
      (*triloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *triloopEnthalpies = (struct triloop*) safe_realloc(*triloopEnthalpies, size * sizeof(struct triloop), o);
      }
   }
   *triloopEnthalpies = (struct triloop*) safe_realloc(*triloopEnthalpies, *num * sizeof(struct triloop), o);
}

static void
getTetraloop(struct tetraloop** tetraloopEntropies, struct tetraloop** tetraloopEnthalpies, int* num, const thal_parameters *tp, thal_results* o)
{
   int i, size;
   double value;
   char *pt_ds = tp->tetraloop_ds;
   *num = 0;
   size = 16;
   if (*tetraloopEntropies != NULL) {
     free(*tetraloopEntropies);
     *tetraloopEntropies = NULL;
   }
   *tetraloopEntropies = (struct tetraloop*) safe_calloc(16, sizeof(struct tetraloop), o);
   while (readTLoop(&pt_ds, (*tetraloopEntropies)[*num].loop, &value, 0, o) != -1) {
      for (i = 0; i < 6; ++i)
        (*tetraloopEntropies)[*num].loop[i] = str2int((*tetraloopEntropies)[*num].loop[i]);
      (*tetraloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *tetraloopEntropies = (struct tetraloop*) safe_realloc(*tetraloopEntropies, size * sizeof(struct tetraloop), o);
      }
   }
   *tetraloopEntropies = (struct tetraloop*) safe_realloc(*tetraloopEntropies, *num * sizeof(struct tetraloop), o);

   char *pt_dh = tp->tetraloop_dh;
   *num = 0;
   size = 16;
   if (*tetraloopEnthalpies != NULL) {
     free(*tetraloopEnthalpies);
     *tetraloopEnthalpies = NULL;
   }
   *tetraloopEnthalpies = (struct tetraloop*) safe_calloc(16, sizeof(struct tetraloop), o);
   while (readTLoop(&pt_dh, (*tetraloopEnthalpies)[*num].loop, &value, 0, o) != -1) {
      for (i = 0; i < 6; ++i)
        (*tetraloopEnthalpies)[*num].loop[i] = str2int((*tetraloopEnthalpies)[*num].loop[i]);
      (*tetraloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *tetraloopEnthalpies = (struct tetraloop*) safe_realloc(*tetraloopEnthalpies, size * sizeof(struct tetraloop), o);
      }
   }
   *tetraloopEnthalpies = (struct tetraloop*) safe_realloc(*tetraloopEnthalpies, *num * sizeof(struct tetraloop), o);
}

static void
tableStartATS(double atp_value, double atpS[5][5])
{

   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpS[i][j] = 0.00000000001;
   atpS[0][3] = atpS[3][0] = atp_value;
}


static void
tableStartATH(double atp_value, double atpH[5][5])
{

   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpH[i][j] = 0.0;

   atpH[0][3] = atpH[3][0] = atp_value;
}

static int
comp3loop(const void* loop1, const void* loop2)
{

     int i;
     const unsigned char* h1 = (const unsigned char*) loop1;
     const struct triloop *h2 = (const struct triloop*) loop2;

     for (i = 0; i < 5; ++i)
         if (h1[i] < h2->loop[i])
             return -1;
       else if (h1[i] > h2->loop[i])
           return 1;

     return 0;
}

static int
comp4loop(const void* loop1, const void* loop2)
{
   int i;
   const unsigned char* h1 = (const unsigned char*) loop1;
   const struct tetraloop *h2 = (const struct tetraloop*) loop2;

   for (i = 0; i < 6; ++i)
     if (h1[i] < h2->loop[i])
       return -1;
   else if (h1[i] > h2->loop[i])
     return 1;

   return 0;
}


static void
initMatrix()
{
   int i, j;
   for (i = 1; i <= len1; ++i) {
      for (j = 1; j <= len2; ++j) {
         if (bpIndx(numSeq1[i], numSeq2[j]) == 0)  {
            EnthalpyDPT(i, j) = _INFINITY;
            EntropyDPT(i, j) = -1.0;
         } else {
            EnthalpyDPT(i, j) = 0.0;
            EntropyDPT(i, j) = MinEntropy;
         }
      }
   }
}

static void
initMatrix2()
{
   int i, j;
   for (i = 1; i <= len1; ++i)
     for (j = i; j <= len2; ++j)
       if (j - i < MIN_HRPN_LOOP + 1 || (bpIndx(numSeq1[i], numSeq1[j]) == 0)) {
          EnthalpyDPT(i, j) = _INFINITY;
          EntropyDPT(i, j) = -1.0;
       } else {
          EnthalpyDPT(i, j) = 0.0;
          EntropyDPT(i, j) = MinEntropy;
       }

}

static void
fillMatrix(int maxLoop, thal_results *o)
{
   int d, i, j, ii, jj;
   double* SH;

   SH = (double*) safe_malloc(2 * sizeof(double), o);
   for (i = 1; i <= len1; ++i) {
      for (j = 1; j <= len2; ++j) {
         if(isFinite(EnthalpyDPT(i, j))) { /* if finite */
            SH[0] = -1.0;
            SH[1] = _INFINITY;
            LSH(i,j,SH);
            if(isFinite(SH[1])) {
               EntropyDPT(i,j) = SH[0];
               EnthalpyDPT(i,j) = SH[1];
            }
            if (i > 1 && j > 1) {
               maxTM(i, j); /* stack: sets EntropyDPT(i, j) and EnthalpyDPT(i, j) */
               for(d = 3; d <= maxLoop + 2; d++) { /* max=30, length over 30 is not allowed */
                  ii = i - 1;
                  jj = - ii - d + (j + i);
                  if (jj < 1) {
                     ii -= abs(jj-1);
                     jj = 1;
                  }
                  for (; ii > 0 && jj < j; --ii, ++jj) {
                     if (isFinite(EnthalpyDPT(ii, jj))) {
                        SH[0] = -1.0;
                        SH[1] = _INFINITY;
                        calc_bulge_internal(ii, jj, i, j, SH,0,maxLoop);
                        if(SH[0] < MinEntropyCutoff) {
                           /* to not give dH any value if dS is unreasonable */
                           SH[0] = MinEntropy;
                           SH[1] = 0.0;
                        }
                        if(isFinite(SH[1])) {
                           EnthalpyDPT(i, j) = SH[1];
                           EntropyDPT(i, j) = SH[0];
                        }
                     }
                  }
               }
            } /* if */
         }
      } /* for */
   } /* for */
   free(SH);
}

static void
fillMatrix2(int maxLoop, thal_results* o)
{
   int i, j;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double), o);
   for (j = 2; j <= len2; ++j)
     for (i = j - MIN_HRPN_LOOP - 1; i >= 1; --i) {
        if (isFinite(EnthalpyDPT(i, j))) {
           SH[0] = -1.0;
           SH[1] = _INFINITY;
           maxTM2(i,j); /* calculate stack */
           CBI(i, j, SH, 0,maxLoop); /* calculate Bulge and Internal loop and stack */
           SH[0] = -1.0;
           SH[1] = _INFINITY;
           calc_hairpin(i, j, SH, 0);
           if(isFinite(SH[1])) {
              if(SH[0] < MinEntropyCutoff){ /* to not give dH any value if dS is unreasonable */
                 SH[0] = MinEntropy;
                 SH[1] = 0.0;
              }
              EntropyDPT(i,j) = SH[0];
              EnthalpyDPT(i,j) = SH[1];
           }
        }
     }
   free(SH);
}


static void
maxTM(int i, int j)
{
   double T0, T1;
   double S0, S1;
   double H0, H1;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double), 0);
   T0 = T1 = -_INFINITY;
   S0 = EntropyDPT(i, j);
   H0 = EnthalpyDPT(i, j);
   RSH(i,j,SH);
   T0 = (H0 + dplx_init_H + SH[1]) /(S0 + dplx_init_S + SH[0] + RC); /* at current position */
   if(isFinite(EnthalpyDPT(i - 1, j - 1)) && isFinite(Hs(i - 1, j - 1, 1))) {
      S1 = (EntropyDPT(i - 1, j - 1) + Ss(i - 1, j - 1, 1));
      H1 = (EnthalpyDPT(i - 1, j - 1) + Hs(i - 1, j - 1, 1));
      T1 = (H1 + dplx_init_H + SH[1]) /(S1 + dplx_init_S + SH[0] + RC);
   } else {
      S1 = -1.0;
      H1 = _INFINITY;
      T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
   }

   if(S1 < MinEntropyCutoff) {
      /* to not give dH any value if dS is unreasonable */
      S1 = MinEntropy;
      H1 = 0.0;
   }
   if(S0 < MinEntropyCutoff) {
      /* to not give dH any value if dS is unreasonable */
      S0 = MinEntropy;
      H0 = 0.0;
   }
   if(T1 > T0) {
      EntropyDPT(i, j) = S1;
      EnthalpyDPT(i, j) = H1;
   } else if(T0 >= T1) {
      EntropyDPT(i, j) = S0;
      EnthalpyDPT(i, j) = H0;
   }
   free(SH);
}

static void
maxTM2(int i, int j)
{
   double T0, T1;
   double S0, S1;
   double H0, H1;
   T0 = T1 = -_INFINITY;
   S0 = EntropyDPT(i, j);
   H0 = EnthalpyDPT(i, j);
   T0 = (H0 + dplx_init_H) /(S0 + dplx_init_S + RC);
   if(isFinite(EnthalpyDPT(i, j))) {
      S1 = (EntropyDPT(i + 1, j - 1) + Ss(i, j, 2));
      H1 = (EnthalpyDPT(i + 1, j - 1) + Hs(i, j, 2));
   } else {
      S1 = -1.0;
      H1 = _INFINITY;
   }
   T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
   if(S1 < MinEntropyCutoff) {
      S1 = MinEntropy;
      H1 = 0.0;
   }
   if(S0 < MinEntropyCutoff) {
      S0 = MinEntropy;
      H0 = 0.0;
   }

   if(T1 > T0) {
      EntropyDPT(i, j) = S1;
      EnthalpyDPT(i, j) = H1;
   } else {
      EntropyDPT(i, j) = S0;
      EnthalpyDPT(i, j) = H0;
   }
}


static void
LSH(int i, int j, double* EntropyEnthalpy)
{
   double S1, H1, T1, G1;
   double S2, H2, T2, G2;
   S1 = S2 = -1.0;
   H1 = H2 = -_INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {
      EntropyDPT(i, j) = -1.0;
      EnthalpyDPT(i, j) = _INFINITY;
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(!isFinite(H1) || G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }
   /** If there is two dangling ends at the same end of duplex **/
   if((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1 ) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
        dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
        dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1<T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0) {
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   } else if ((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1<T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if (G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   } else if ((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2  && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0) {
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   G1 = H1 -TEMP_KELVIN*S1;
   G2 = H2 -TEMP_KELVIN*S2;
   if(isFinite(H1)) {
      if(T1 < T2) {
         EntropyEnthalpy[0] = S2;
         EntropyEnthalpy[1] = H2;
      } else {
         EntropyEnthalpy[0] = S1;
         EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}

static void
RSH(int i, int j, double* EntropyEnthalpy)
{
   double G1, G2;
   double S1, S2;
   double H1, H2;
   double T1, T2;
   S1 = S2 = -1.0;
   H1 = H2 = _INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(!isFinite(H1) || G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }

   if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
    G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }

      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }

   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2 >0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }

   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if (G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   G1 = H1 -TEMP_KELVIN*S1;
   G2 =  H2 -TEMP_KELVIN*S2;
   if(isFinite(H1)) {
      if(T1 < T2) {
         EntropyEnthalpy[0] = S2;
         EntropyEnthalpy[1] = H2;
      } else {
         EntropyEnthalpy[0] = S1;
         EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}

static double
Ss(int i, int j, int k)
{
   if(k==2) {
      if (i >= j)
        return -1.0;
      if (i == len1 || j == len2 + 1)
        return -1.0;

      if (i > len1)
        i -= len1;
      if (j > len2)
        j -= len2;
      return stackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
   } else {
      return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}


static double
Hs(int i, int j, int k)
{
   if(k==2) {
      if (i >= j)
        return _INFINITY;
      if (i == len1 || j == len2 + 1)
        return _INFINITY;

      if (i > len1)
        i -= len1;
      if (j > len2)
        j -= len2;
      if(isFinite(stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]])) {
         return stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
      } else {
         return _INFINITY;
      }
   } else {
      return stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}

static void
CBI(int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop)
{
   int d, ii, jj;
   for (d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop; --d)
     for (ii = i + 1; ii < j - d && ii <= len1; ++ii) {
        jj = d + ii;
        if(traceback==0) {
           EntropyEnthalpy[0] = -1.0;
           EntropyEnthalpy[1] = _INFINITY;
        }
        if (isFinite(EnthalpyDPT(ii, jj)) && isFinite(EnthalpyDPT(i, j))) {
           calc_bulge_internal2(i, j, ii, jj, EntropyEnthalpy, traceback,maxLoop);
           if(isFinite(EntropyEnthalpy[1])) {
              if(EntropyEnthalpy[0] < MinEntropyCutoff) {
                 EntropyEnthalpy[0] = MinEntropy;
                 EntropyEnthalpy[1] = 0.0;
              }
              if(traceback==0) {
                 EnthalpyDPT(i, j) = EntropyEnthalpy[1];
                 EntropyDPT(i, j) = EntropyEnthalpy[0];
              }
           }
        }
     }
   return;
}

static void
calc_hairpin(int i, int j, double* EntropyEnthalpy, int traceback)
{
   int loopSize = j - i - 1;
   double G1, G2;
   G1 = G2 = -_INFINITY;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double), 0);
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   if(loopSize < MIN_HRPN_LOOP) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   if (i <= len1 && len2 < j) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   } else if (i > len2) {
      i -= len1;
      j -= len2;
   }
   if(loopSize <= 30) {
      EntropyEnthalpy[1] = hairpinLoopEnthalpies[loopSize - 1];
      EntropyEnthalpy[0] = hairpinLoopEntropies[loopSize - 1];
   } else {
      EntropyEnthalpy[1] = hairpinLoopEnthalpies[29];
      EntropyEnthalpy[0] = hairpinLoopEntropies[29];
   }

   if (loopSize > 3) { /* for loops 4 bp and more in length, terminal mm are accounted */
      EntropyEnthalpy[1] += tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
      EntropyEnthalpy[0] += tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
   } else if(loopSize == 3){ /* for loops 3 bp in length at-penalty is considered */
      EntropyEnthalpy[1] += atPenaltyH(numSeq1[i], numSeq1[j]);
      EntropyEnthalpy[0] += atPenaltyS(numSeq1[i], numSeq1[j]);
   }

   if (loopSize == 3) {         /* closing AT-penalty (+), triloop bonus, hairpin of 3 (+) */
      struct triloop* loop;
      if (numTriloops) {
         if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEnthalpies, numTriloops, sizeof(struct triloop), comp3loop)))
           EntropyEnthalpy[1] += loop->value;
         if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEntropies, numTriloops, sizeof(struct triloop), comp3loop)))
           EntropyEnthalpy[0] += loop->value;
      }
   } else if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
      struct tetraloop* loop;
      if (numTetraloops) {
         if ((loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEnthalpies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
            EntropyEnthalpy[1] += loop->value;
         }
         if ((loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEntropies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
            EntropyEnthalpy[0] += loop->value;
         }
      }
   }
   if(!isFinite(EntropyEnthalpy[1])) {
      EntropyEnthalpy[1] = _INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   if(isPositive(EntropyEnthalpy[1]) && isPositive(EntropyEnthalpy[0]) && (!isPositive(EnthalpyDPT(i, j)) || !isPositive(EntropyDPT(i, j)))) { /* if both, S and H are positive */
      EntropyEnthalpy[1] = _INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   RSH(i,j,SH);
   G1 = EntropyEnthalpy[1]+SH[1] -TEMP_KELVIN*(EntropyEnthalpy[0]+SH[0]);
   G2 = EnthalpyDPT(i, j)+SH[1] -TEMP_KELVIN*(EntropyDPT(i, j)+SH[0]);
     if(G2 < G1 && traceback == 0) {
      EntropyEnthalpy[0] = EntropyDPT(i, j);
      EntropyEnthalpy[1] = EnthalpyDPT(i, j);
   }
   free(SH);
   return;
}


static void
calc_bulge_internal(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback, int maxLoop)
{
   int loopSize1, loopSize2, loopSize;
   double S,H,G1,G2;
   int N, N_loop;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double), 0);
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   S = -1.0;
   H = _INFINITY;
   loopSize1 = ii - i - 1;
   loopSize2 = jj - j - 1;
   if(ii < jj) {
      N = ((2 * i)/2);
      N_loop = N;
      if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
      if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
   } else {
      N = ((2 * j)/2);
      N_loop = 2 * jj;
      if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
      if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
      N_loop = (N_loop/2) - 1;
   }
#ifdef DEBUG
   if (ii <= i){
      fputs("Error in calc_bulge_internal(): ii is not greater than i\n", stderr);
   }
   if (jj <= j)
     fputs("Error in calc_bulge_internal(): jj is not greater than j\n", stderr);
#endif

#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      free(SH);
      return;
   }
#endif
#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      free(SH);
      return;
   }
#endif
   loopSize = loopSize1 + loopSize2-1;
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                              the intervening nn-pair must be added */

         if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
            H = bulgeLoopEnthalpies[loopSize] +
              stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
            S = bulgeLoopEntropies[loopSize] +
              stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
         }
         if(isPositive(H) || isPositive(S)){
            H = _INFINITY;
            S = -1.0;
         }
         H += EnthalpyDPT(i, j);
         S += EntropyDPT(i, j);
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
         RSH(ii,jj,SH);
         G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
         G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*((EntropyDPT(ii, jj)+SH[0]));
         if((G1< G2) || (traceback==1)) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

         H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
         H += EnthalpyDPT(i, j);

         S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
         S += EntropyDPT(i, j);
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
         if(isPositive(H) && isPositive(S)){
            H = _INFINITY;
            S = -1.0;
         }

     RSH(ii,jj,SH);
         G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
         G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
         if(G1< G2 || (traceback==1)){
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }

      }
   } else if (loopSize1 == 1 && loopSize2 == 1) {
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        stackint2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      S += EntropyDPT(i, j);

      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        stackint2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      H += EnthalpyDPT(i, j);
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
      if(isPositive(H) && isPositive(S)){
         H = _INFINITY;
         S = -1.0;
      }
     RSH(ii,jj,SH);
         G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
      G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
           if((G1< G2) || traceback==1) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      free(SH);
      return;
   } else { /* only internal loops */
      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        tstackEnthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]
        + (ILAH * abs(loopSize1 - loopSize2));
      H += EnthalpyDPT(i, j);

      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        tstackEntropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      S += EntropyDPT(i, j);
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
   if(isPositive(H) && isPositive(S)){
         H = _INFINITY;
         S = -1.0;
      }
     RSH(ii,jj,SH);
     G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
     G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
     if((G1< G2) || (traceback==1)){
             EntropyEnthalpy[0] = S;
             EntropyEnthalpy[1] = H;
      }
   }
   free(SH);
   return;
}

static void
calc_bulge_internal2(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback, int maxLoop)
{
   int loopSize1, loopSize2, loopSize;
   double T1, T2;
   double S,H;
   /* int N, N_loop; Triinu, please review */
   T1 = T2 = -_INFINITY;
   S = MinEntropy;
   H = 0.0;
   loopSize1 = ii - i - 1;
   loopSize2 = j - jj - 1;
   if (loopSize1 + loopSize2 > maxLoop) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   /* Triinu, please review the statements below. */
   /* if(i < (len1 -j)) { */
     /* N  = i; */
      /* N_loop = (i - 1); */
   /* } else { */
     /* N = len1-j;  */
      /* N_loop = len1 - j - 1; */
   /* } */
#ifdef DEBUG
   if (ii <= i)
     fputs("Error in calc_bulge_internal(): ii isn't greater than i\n", stderr);
   if (jj >= j)
     fputs("Error in calc_bulge_internal(): jj isn't less than j\n", stderr);
   if (ii >= jj)
     fputs("Error in calc_bulge_internal(): jj isn't greater than ii\n", stderr);

   if ((i <= len1 && len1 < ii) || (jj <= len2 && len2 < j))  {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
#endif

#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      return;
   }
#endif
#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      return;
   }
#endif

#ifdef DEBUG
   if (i > len1)
     i -= len1;
   if (ii > len1)
     ii -= len1;
   if (j > len2)
     j -= len2;
   if (jj > len2)
     jj -= len2;
#endif
   loopSize = loopSize1 + loopSize2 -1; /* for indx only */
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                              the intervening nn-pair must be added */
         if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
            H = bulgeLoopEnthalpies[loopSize] +
              stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
            S = bulgeLoopEntropies[loopSize] +
              stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
         }
         if(traceback!=1) {
            H += EnthalpyDPT(ii, jj); /* bulge koos otsaga, st bulge i,j-ni */
            S += EntropyDPT(ii, jj);
         }

         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }

         T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
         T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);

         if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }

      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

         H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
         if(traceback!=1)
           H += EnthalpyDPT(ii, jj);

         S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
         if(traceback!=1)
           S += EntropyDPT(ii, jj);
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }

         T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
         T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);

         if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      }
   } /* end of calculating bulges */
   else if (loopSize1 == 1 && loopSize2 == 1) {
      /* mismatch nearest neighbor parameters */

      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
        stackint2Entropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
        S += EntropyDPT(ii, jj);

      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
        stackint2Enthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
        H += EnthalpyDPT(ii, jj);
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }

      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);
      if((DBL_EQ(T1,T2) == 2) || traceback) {
         if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      }
      return;
   } else { /* only internal loops */

      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
        tstackEnthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]]
        + (ILAH * abs(loopSize1 - loopSize2));
      if(traceback!=1)
        H += EnthalpyDPT(ii, jj);

      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
        tstackEntropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      if(traceback!=1)
        S += EntropyDPT(ii, jj);
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }

      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);
      if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
         EntropyEnthalpy[0] = S;
         EntropyEnthalpy[1] = H;
      }
   }
   return;
}

static void
calc_terminal_bp(double temp) { /* compute exterior loop */
   int i;
   int max;
   SEND5(0) = SEND5(1) = -1.0;
   HEND5(0) = HEND5(1) = _INFINITY;
   for(i = 2; i<=(len1); i++) {
      SEND5(i) = MinEntropy;
      HEND5(i) = 0;
   }

   double T1, T2, T3, T4, T5;
   T1 = T2 = T3 = T4 = T5 = -_INFINITY;
   double G;
   /* adding terminal penalties to 3' end and to 5' end */
   for(i = 2; i <= len1; ++i) {
      max = 0;
      T1 = T2 = T3 = T4 = T5 = -_INFINITY;
      T1 = (HEND5(i - 1) + dplx_init_H) / (SEND5(i - 1) + dplx_init_S + RC);
      T2 = (END5_1(i,1) + dplx_init_H) / (END5_1(i,2) + dplx_init_S + RC);
      T3 = (END5_2(i,1) + dplx_init_H) / (END5_2(i,2) + dplx_init_S + RC);
      T4 = (END5_3(i,1) + dplx_init_H) / (END5_3(i,2) + dplx_init_S + RC);
      T5 = (END5_4(i,1) + dplx_init_H) / (END5_4(i,2) + dplx_init_S + RC);
      max = max5(T1,T2,T3,T4,T5);
      switch (max) {
       case 1:
         SEND5(i) = SEND5(i - 1);
         HEND5(i) = HEND5(i - 1);
         break;
       case 2:
         G = END5_1(i,1) - (temp * (END5_1(i,2)));
         if(G < G2) {
            SEND5(i) = END5_1(i,2);
            HEND5(i) = END5_1(i,1);
         } else {
            SEND5(i) = SEND5(i - 1);
            HEND5(i) = HEND5(i - 1);
         }
         break;
       case 3:
         G = END5_2(i,1) - (temp * (END5_2(i,2)));
         if(G < G2) {
            SEND5(i) = END5_2(i,2);
            HEND5(i) = END5_2(i,1);
         } else {
            SEND5(i) = SEND5(i - 1);
            HEND5(i) = HEND5(i - 1);
         }
         break;
       case 4:
         G = END5_3(i,1) - (temp * (END5_3(i,2)));
         if(G < G2) {
            SEND5(i) = END5_3(i,2);
            HEND5(i) = END5_3(i,1);
         } else {
            SEND5(i) = SEND5(i - 1);
            HEND5(i) = HEND5(i - 1);
         }
         break;
       case 5:
         G = END5_4(i,1) - (temp * (END5_4(i,2)));
         if(G < G2) {
            SEND5(i) = END5_4(i,2);
            HEND5(i) = END5_4(i,1);
         } else {
            SEND5(i) = SEND5(i - 1);
            HEND5(i) = HEND5(i - 1);
         }
         break;
       default:
#ifdef DEBUG
         printf ("WARNING: max5 returned character code %d ??\n", max);
#endif
         break;
      }
   }
}

static double
END5_1(int i,int hs)
{
   int k;
   double max_tm; /* energy min */
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   S_max = S = -1.0;
   T1 = T2 = -_INFINITY;
   max_tm = -_INFINITY;
   for(k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
         H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
         S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
         if(!isFinite(H) || H > 0 || S > 0) { /* H and S must be greater than 0 to avoid BS */
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
         H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
         S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
         if(S > MinEntropyCutoff) {
            H_max = H;
            S_max = S;
            max_tm = T1;
         }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double
END5_2(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
         H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
         S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
         H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
         S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
         if(S > MinEntropyCutoff) {
            H_max = H;
            S_max = S;
            max_tm = T1;
         }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double
END5_3(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
         H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
         S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
         H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
         S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
         if(S > MinEntropyCutoff) {
            H_max = H;
            S_max = S;
            max_tm = T1;
         }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double
END5_4(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for(k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
         H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
         S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
         H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
         S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
         if(S > MinEntropyCutoff) {
            H_max = H;
            S_max = S;
            max_tm = T1;
         }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}


static double
Sd5(int i, int j)
{
   return dangleEntropies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
}

static double
Hd5(int i, int j)
{
   return dangleEnthalpies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
}

static double
Sd3(int i, int j)
{
   return dangleEntropies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
}

static double
Hd3(int i, int j)
{
   return dangleEnthalpies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
}

static double
Ststack(int i, int j)
{
   return tstack2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
}

static double
Htstack(int i, int j)
{ /* e.g AG_TC 210 */
   return tstack2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
}

/* Return 1 if string is symmetrical, 0 otherwise. */
static int
symmetry_thermo(const unsigned char* seq)
{
   register char s;
   register char e;
   const unsigned char *seq_end=seq;
   int i = 0;
   int seq_len=length_unsig_char(seq);
   int mp = seq_len/2;
   if(seq_len%2==1) {
      return 0;
   }
   seq_end+=seq_len;
   seq_end--;
   while(i<mp) {
      i++;
      s=toupper(*seq);
      e=toupper(*seq_end);
      if ((s=='A' && e!='T')
          || (s=='T' && e!='A')
          || (e=='A' && s!='T')
          || (e=='T' && s!='A')) {
         return 0;
      }
      if ((s=='C' && e!='G')
          || (s=='G' && e!='C')
          || (e=='C' && s!='G')
          || (e=='G' && s!='C')) {
         return 0;
      }
      seq++;
      seq_end--;
   }
   return 1;
}

static int
length_unsig_char(const unsigned char * str)
{
   int i = 0;
   while(*(str++)) {
      i++;
      if(i == INT_MAX)
        return -1;
   }
   return i;
}

static void
tracebacku(int* bp, int maxLoop,thal_results* o) /* traceback for unimolecular structure */
{
   int i, j;
   i = j = 0;
   int ii, jj, k;
   struct tracer *top, *stack = NULL;
   double* SH1;
   double* SH2;
   double* EntropyEnthalpy;
   SH1 = (double*) safe_malloc(2 * sizeof(double), o);
   SH2 = (double*) safe_malloc(2 * sizeof(double), o);
   EntropyEnthalpy = (double*) safe_malloc(2 * sizeof(double), o);
   push(&stack,len1, 0, 1, o);
   while(stack) {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;
      if(top->mtrx==1) {
         while (equal(SEND5(i), SEND5(i - 1)) && equal(HEND5(i), HEND5(i - 1))) /* if previous structure is the same as this one */
           --i;
         if (i == 0)
           continue;
         if (equal(SEND5(i), END5_1(i,2)) && equal(HEND5(i), END5_1(i,1))) {
            for (k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k)
              if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i)) &&
                  equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i))) {
                 push(&stack, k + 1, i,0, o);
                 break;
              }
            else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i)) &&
                     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i))) {
               push(&stack, k + 1, i, 0, o);
               push(&stack, k, 0, 1, o);
               break;
            }
         }
         else if (equal(SEND5(i), END5_2(i,2)) && equal(HEND5(i), END5_2(i,1))) {
            for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
              if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i)) &&
                  equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i))) {
                 push(&stack, k + 2, i, 0, o);
                 break;
              }
            else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i)) &&
                     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i))) {
               push(&stack, k + 2, i, 0, o);
               push(&stack, k, 0, 1, o);
               break;
            }
         }
         else if (equal(SEND5(i), END5_3(i,2)) && equal(HEND5(i), END5_3(i,1))) {
            for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
              if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1))
                  && equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
                 push(&stack, k + 1, i - 1, 0, o);
                 break;
              }
            else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1)) &&
                     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
               push(&stack, k + 1, i - 1, 0, o); /* matrix 0  */
               push(&stack, k, 0, 1, o); /* matrix 3 */
               break;
            }
         }
         else if(equal(SEND5(i), END5_4(i,2)) && equal(HEND5(i), END5_4(i,1))) {
            for (k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k)
              if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
                  equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1))) {
                 push(&stack, k + 2, i - 1, 0, o);
                 break;
              }
            else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
                     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1)) ) {
               push(&stack, k + 2, i - 1, 0, o);
               push(&stack, k, 0, 1, o);
               break;
            }
         }
      }
      else if(top->mtrx==0) {
         bp[i - 1] = j;
         bp[j - 1] = i;
         SH1[0] = -1.0;
         SH1[1] = _INFINITY;
         calc_hairpin(i, j, SH1, 1); /* 1 means that we use this method in traceback */
         SH2[0] = -1.0;
         SH2[1] = _INFINITY;
         CBI(i,j,SH2,2,maxLoop);
         if (equal(EntropyDPT(i, j), Ss(i, j, 2) + EntropyDPT(i + 1, j - 1)) &&
             equal(EnthalpyDPT(i, j), Hs(i, j, 2) + EnthalpyDPT(i + 1, j - 1))) {
            push(&stack, i + 1, j - 1, 0, o);
         }
         else if (equal(EntropyDPT(i, j), SH1[0]) && equal(EnthalpyDPT(i,j), SH1[1]));
         else if (equal(EntropyDPT(i, j), SH2[0]) && equal(EnthalpyDPT(i, j), SH2[1])) {
            int d, done;
            for (done = 0, d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop && !done; --d)
              for (ii = i + 1; ii < j - d; ++ii) {
                 jj = d + ii;
                 EntropyEnthalpy[0] = -1.0;
                 EntropyEnthalpy[1] = _INFINITY;
                 calc_bulge_internal2(i, j, ii, jj,EntropyEnthalpy,1,maxLoop);
                 if (equal(EntropyDPT(i, j), EntropyEnthalpy[0] + EntropyDPT(ii, jj)) &&
                     equal(EnthalpyDPT(i, j), EntropyEnthalpy[1] + EnthalpyDPT(ii, jj))) {
                    push(&stack, ii, jj, 0, o);
                    ++done;
                    break;
                 }
              }
         } else {
         }
      }
      free(top);
   }
   free(SH1);
   free(SH2);
   free(EntropyEnthalpy);
}


static void
traceback(int i, int j, double RT, int* ps1, int* ps2, int maxLoop, thal_results* o)
{
   int d, ii, jj, done;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double), o);
   ps1[i - 1] = j;
   ps2[j - 1] = i;
   while(1) {
      SH[0] = -1.0;
      SH[1] = _INFINITY;
      LSH(i,j,SH);
      if(equal(EntropyDPT(i,j),SH[0]) && equal(EnthalpyDPT(i,j),SH[1])) {
         break;
      }
      done = 0;
      if (i > 1 && j > 1 && equal(EntropyDPT(i,j), Ss(i - 1, j - 1, 1) + EntropyDPT(i - 1, j - 1)) && equal(EnthalpyDPT(i,j), Hs(i - 1, j - 1, 1) + EnthalpyDPT(i - 1, j - 1))) {
         i = i - 1;
         j = j - 1;
         ps1[i - 1] = j;
         ps2[j - 1] = i;
         done = 1;
      }
      for (d = 3; !done && d <= maxLoop + 2; ++d) {
         ii = i - 1;
         jj = -ii - d + (j + i);
         if (jj < 1) {
            ii -= abs(jj-1);
            jj = 1;
         }
         for (; !done && ii > 0 && jj < j; --ii, ++jj) {
            SH[0] = -1.0;
            SH[1] = _INFINITY;
            calc_bulge_internal(ii, jj, i, j, SH,1,maxLoop);
            if (equal(EntropyDPT(i, j), SH[0]) && equal(EnthalpyDPT(i, j), SH[1])) {
               i = ii;
               j = jj;
               ps1[i - 1] = j;
               ps2[j - 1] = i;
               done = 1;
               break;
            }
         }
      }
   }
   free(SH);
}

char *
drawDimer(int* ps1, int* ps2, double temp, double H, double S, const thal_mode mode, double t37, thal_results *o)
{
   int  ret_space = 0;
   char *ret_ptr = NULL;
   int ret_nr, ret_pr_once;
   char ret_para[400];
   char* ret_str[4];
   int i, j, k, numSS1, numSS2, N;
   char* duplex[4];
   double G, t;
   t = G = 0;
   if (!isFinite(temp)){
      if((mode != THL_FAST) && (mode != THL_DEBUG_F) && (mode != THL_STRUCT)) {
         printf("No predicted secondary structures for given sequences\n");
      }
      o->temp = 0.0; /* lets use generalization here; this should rather be very negative value */
      strcpy(o->msg, "No predicted sec struc for given seq");
      return NULL;
   } else {
      N=0;
      for(i=0;i<len1;i++){
         if(ps1[i]>0) ++N;
      }
      for(i=0;i<len2;i++) {
         if(ps2[i]>0) ++N;
      }
      N = (N/2) -1;
      t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
      S = S + (N * saltCorrection);
      G = (H) - (t37 * S);
      o->temp = (double) t;
      o->H = (double) H;
      o->S = (double) S;
      o->dG = (double) G;
      
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
         o->temp = (double) t;
         /* maybe user does not need as precise as that */
         /* printf("Thermodynamical values:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\tN = %d, SaltC=%f, RC=%f\n",
                len1, (double) S, (double) H, (double) G, (double) t, (int) N, saltCorrection, RC); */
         if (mode != THL_STRUCT) {
           printf("Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
                  (double) S, (double) H, (double) G, (double) t);
         } else {
           sprintf(ret_para, "Tm: %.1f&deg;C  dG: %.0f cal/mol  dH: %.0f cal/mol  dS: %.0f cal/mol*K\\n",
                   (double) t, (double) G, (double) H, (double) S);
         }
      } else {
        /* do nothing, return later */
        // return NULL;
      }
   }

   duplex[0] = (char*) safe_malloc(len1 + len2 + 1, o);
   duplex[1] = (char*) safe_malloc(len1 + len2 + 1, o);
   duplex[2] = (char*) safe_malloc(len1 + len2 + 1, o);
   duplex[3] = (char*) safe_malloc(len1 + len2 + 1, o);
   duplex[0][0] = duplex[1][0] = duplex[2][0] = duplex[3][0] = 0;

   i = 0;
   numSS1 = 0;
   while (ps1[i++] == 0) ++numSS1;
   j = 0;
   numSS2 = 0;
   while (ps2[j++] == 0) ++numSS2;

   if (numSS1 >= numSS2){
      for (i = 0; i < numSS1; ++i) {
         strcatc(duplex[0], oligo1[i]);
         strcatc(duplex[1], ' ');
         strcatc(duplex[2], ' ');
      }
      for (j = 0; j < numSS1 - numSS2; ++j) strcatc(duplex[3], ' ');
      for (j = 0; j < numSS2; ++j) strcatc(duplex[3], oligo2[j]);
   } else {
      for (j = 0; j < numSS2; ++j) {
         strcatc(duplex[3], oligo2[j]);
         strcatc(duplex[1], ' ');
         strcatc(duplex[2], ' ');
      }
      for (i = 0; i < numSS2 - numSS1; ++i)
        strcatc(duplex[0], ' ');
      for (i = 0; i < numSS1; ++i)
        strcatc(duplex[0], oligo1[i]);
   }
   i = numSS1 + 1;
   j = numSS2 + 1;

   while (i <= len1) {
      while (i <= len1 && ps1[i - 1] != 0 && j <= len2 && ps2[j - 1] != 0) {
         strcatc(duplex[0], ' ');
         strcatc(duplex[1], oligo1[i - 1]);
         strcatc(duplex[2], oligo2[j - 1]);
         strcatc(duplex[3], ' ');
         ++i;
         ++j;
      }
      numSS1 = 0;
      while (i <= len1 && ps1[i - 1] == 0) {
         strcatc(duplex[0], oligo1[i - 1]);
         strcatc(duplex[1], ' ');
         ++numSS1;
         ++i;
      }
      numSS2 = 0;
      while (j <= len2 && ps2[j - 1] == 0) {
         strcatc(duplex[2], ' ');
         strcatc(duplex[3], oligo2[j - 1]);
         ++numSS2;
         ++j;
      }
      if (numSS1 < numSS2)
        for (k = 0; k < numSS2 - numSS1; ++k) {
           strcatc(duplex[0], '-');
           strcatc(duplex[1], ' ');
        }
      else if (numSS1 > numSS2)
        for (k = 0; k < numSS1 - numSS2; ++k) {
           strcatc(duplex[2], ' ');
           strcatc(duplex[3], '-');
        }
   }
   if (mode == THL_FAST) {
     strcpy(o->duplex0, duplex[0]);
     strcpy(o->duplex1, duplex[1]);
     strcpy(o->duplex2, duplex[2]);
     strcpy(o->duplex3, duplex[3]);
     return NULL;
   }
   if ((mode == THL_GENERAL) || (mode == THL_DEBUG)) {
     printf("SEQ\t");
     printf("%s\n", duplex[0]);
     printf("SEQ\t");
     printf("%s\n", duplex[1]);
     printf("STR\t");
     printf("%s\n", duplex[2]);
     printf("STR\t");
     printf("%s\n", duplex[3]);
   }
   if (mode == THL_STRUCT) {
     ret_str[3] = NULL;
     ret_str[0] = (char*) safe_malloc(len1 + len2 + 10, o);
     ret_str[1] = (char*) safe_malloc(len1 + len2 + 10, o);
     ret_str[2] = (char*) safe_malloc(len1 + len2 + 10, o);
     ret_str[0][0] = ret_str[1][0] = ret_str[2][0] = '\0';

     /* Join top primer */
     strcpy(ret_str[0], "   ");
     strcat(ret_str[0], duplex[0]);
     ret_nr = 0;
     while (duplex[1][ret_nr] != '\0') {
       if (duplex[1][ret_nr] == 'A' || duplex[1][ret_nr] == 'T' ||
           duplex[1][ret_nr] == 'C' || duplex[1][ret_nr] == 'G' ||
           duplex[1][ret_nr] == '-') {
         ret_str[0][ret_nr + 3] = duplex[1][ret_nr];
       }
       ret_nr++;
     }
     if (strlen(duplex[1]) > strlen(duplex[0])) {
       ret_str[0][strlen(duplex[1]) + 3] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[0]) - 1;
     while (ret_nr > 0 && (ret_str[0][ret_nr] == ' ' || ret_str[0][ret_nr] == '-')) {
       ret_str[0][ret_nr--] = '\0';
     }
     /* Write the 5' */
     ret_nr = 3;
     ret_pr_once = 1;
     while (ret_str[0][ret_nr] != '\0' && ret_pr_once == 1) {
       if (ret_str[0][ret_nr] == 'A' || ret_str[0][ret_nr] == 'T' ||
           ret_str[0][ret_nr] == 'C' || ret_str[0][ret_nr] == 'G' ||
           ret_str[0][ret_nr] == '-') {
         ret_str[0][ret_nr - 3] = '5';
    ret_str[0][ret_nr - 2] = '\'';
    ret_pr_once = 0;
       }
       ret_nr++;
     }

     /* Create the align tics */
     strcpy(ret_str[1], "     ");
     for (i = 0 ; i < strlen(duplex[1]) ; i++) {
       if (duplex[1][i] == 'A' || duplex[1][i] == 'T' ||
           duplex[1][i] == 'C' || duplex[1][i] == 'G' ) {
    ret_str[1][i + 3] = '|';
       } else {
         ret_str[1][i + 3] = ' ';
       }
       ret_str[1][i + 4] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[1]) - 1;
     while (ret_nr > 0 && ret_str[1][ret_nr] == ' ') {
       ret_str[1][ret_nr--] = '\0';
     }
     /* Join bottom primer */
     strcpy(ret_str[2], "   ");
     strcat(ret_str[2], duplex[2]);
     ret_nr = 0;
     while (duplex[3][ret_nr] != '\0') {
       if (duplex[3][ret_nr] == 'A' || duplex[3][ret_nr] == 'T' ||
           duplex[3][ret_nr] == 'C' || duplex[3][ret_nr] == 'G' ||
           duplex[3][ret_nr] == '-') {
         ret_str[2][ret_nr + 3] = duplex[3][ret_nr];
       }
       ret_nr++;
     }
     if (strlen(duplex[3]) > strlen(duplex[2])) {
       ret_str[2][strlen(duplex[3]) + 3] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[2]) - 1;
     while (ret_nr > 0 && (ret_str[2][ret_nr] == ' ' || ret_str[2][ret_nr] == '-')) {
       ret_str[2][ret_nr--] = '\0';
     }
     /* Write the 5' */
     ret_nr = 3;
     ret_pr_once = 1;
     while (ret_str[2][ret_nr] != '\0' && ret_pr_once == 1) {
       if (ret_str[2][ret_nr] == 'A' || ret_str[2][ret_nr] == 'T' ||
           ret_str[2][ret_nr] == 'C' || ret_str[2][ret_nr] == 'G' ||
           ret_str[2][ret_nr] == '-') {
         ret_str[2][ret_nr - 3] = '3';
         ret_str[2][ret_nr - 2] = '\'';
         ret_pr_once = 0;
       }
       ret_nr++;
     }

     save_append_string(&ret_str[3], &ret_space, o, ret_para);
     save_append_string(&ret_str[3], &ret_space, o, ret_str[0]);
     save_append_string(&ret_str[3], &ret_space, o, " 3\'\\n");
     save_append_string(&ret_str[3], &ret_space, o, ret_str[1]);
     save_append_string(&ret_str[3], &ret_space, o, "\\n");
     save_append_string(&ret_str[3], &ret_space, o, ret_str[2]);
     save_append_string(&ret_str[3], &ret_space, o, " 5\'\\n");


/*
     save_append_string(&ret_str, &ret_space, o, "SEQ ");
     save_append_string(&ret_str, &ret_space, o, duplex[0]);
     save_append_string(&ret_str, &ret_space, o, "\\nSEQ ");
     save_append_string(&ret_str, &ret_space, o, duplex[1]);
     save_append_string(&ret_str, &ret_space, o, "\\nSTR ");
     save_append_string(&ret_str, &ret_space, o, duplex[2]);
     save_append_string(&ret_str, &ret_space, o, "\\nSTR ");
     save_append_string(&ret_str, &ret_space, o, duplex[3]);
     save_append_string(&ret_str, &ret_space, o, "\\n");
*/
     ret_ptr = (char *) safe_malloc(strlen(ret_str[3]) + 1, o);
     strcpy(ret_ptr, ret_str[3]);
     if (ret_str[3]) {
       free(ret_str[3]);
     }
     free(ret_str[0]);
     free(ret_str[1]);
     free(ret_str[2]);
   }
   free(duplex[0]);
   free(duplex[1]);
   free(duplex[2]);
   free(duplex[3]);

   return ret_ptr;
}

char *
drawHairpin(int* bp, double mh, double ms, const thal_mode mode, double temp, thal_results *o)
{
   int  ret_space = 0;
   char *ret_ptr = NULL;
   int ret_last_l, ret_first_r, ret_center, ret_left_end, ret_right_start, ret_left_len, ret_right_len;
   int ret_add_sp_l, ret_add_sp_r;
   char ret_center_char;
   char ret_para[400];
   char* ret_str;
   /* Plain text */
   int i, N;
   N = 0;
   double mg, t;
   if (!isFinite(ms) || !isFinite(mh)) {
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
        if (mode != THL_STRUCT) {
          printf("0\tdS = %g\tdH = %g\tinf\tinf\n", (double) ms,(double) mh);
#ifdef DEBUG
          fputs("No temperature could be calculated\n",stderr);
#endif
        }
      } else {
         o->temp = 0.0; /* lets use generalization here */
         strcpy(o->msg, "No predicted sec struc for given seq\n");
      }
   } else {
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
         for (i = 1; i < len1; ++i) {
            if(bp[i-1] > 0) N++;
         }
      } else {
         for (i = 1; i < len1; ++i) {
            if(bp[i-1] > 0) N++;
         }
      }
      t = (mh / (ms + (((N/2)-1) * saltCorrection))) - ABSOLUTE_ZERO;
      mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));
      ms = ms + (((N/2)-1) * saltCorrection);
      o->temp = (double) t;
      o->H = (double) mh;
      o->S = (double) ms;
      o->dG = (double) mg;

      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
         o->temp = (double) t;
         if (mode != THL_STRUCT) {
           printf("Calculated thermodynamical parameters for dimer:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
                  len1, (double) ms, (double) mh, (double) mg, (double) t);
         } else {
           sprintf(ret_para, "Tm: %.1f&deg;C  dG: %.0f cal/mol  dH: %.0f cal/mol  dS: %.0f cal/mol*K\\n",
                   (double) t, (double) mg, (double) mh, (double) ms);
         }
      }
   }
   /* plain-text output */
   char asciiRow[len1+1];
   
   // asciiRow = (char*) safe_malloc(len1, o);
   for(i = 0; i <= len1; ++i) asciiRow[i] = '\0';
   for(i = 1; i < len1+1; ++i) {
      if(bp[i-1] == 0) {
         asciiRow[(i-1)] = '-';
      } else {
         if(bp[i-1] > (i-1)) {
            asciiRow[(bp[i-1]-1)]='\\';
         } else  {
            asciiRow[(bp[i-1]-1)]='/';
         }
      }
   }
   if (mode == THL_FAST) {
         strcpy(o->duplex0, asciiRow);
         strcpy(o->duplex1, oligo1);
         strcpy(o->duplex2, "");
         strcpy(o->duplex3, "");
         // fputs("return from drawHairpin\n",stderr);
         return NULL;
   }
   
   if ((mode == THL_GENERAL) || (mode == THL_DEBUG)) {
     printf("SEQ\t");
     for(i = 0; i < len1; ++i) printf("%c",asciiRow[i]);
     printf("\nSTR\t%s\n", oligo1);
   }
   if (mode == THL_STRUCT) {
     ret_str = NULL;

     save_append_string(&ret_str, &ret_space, o, ret_para);

     ret_last_l = -1;
     ret_first_r = -1;
     ret_center_char = '|';
     for(i = 0; i < len1; ++i) {
       if (asciiRow[i] == '/') {
         ret_last_l = i;
       }
       if ((ret_first_r == -1) && (asciiRow[i] == '\\')) {
         ret_first_r = i;
       }
     }
     ret_center = ret_first_r - ret_last_l;
     if (ret_center % 2 == 0) {
       /* ret_center is odd */
       ret_left_end = ret_last_l + (ret_first_r - ret_last_l) / 2 - 1;
       ret_center_char = (char) oligo1[ret_left_end + 1];
       ret_right_start = ret_left_end + 2;
     } else {
       /* ret_center is even */
       ret_left_end = ret_last_l + (ret_first_r - ret_last_l - 1) / 2;
       ret_right_start = ret_left_end + 1;
     }
     ret_left_len = ret_left_end + 1;
     ret_right_len = len1 - ret_right_start;
     ret_add_sp_l = 0;
     ret_add_sp_r = 0;
     if (ret_left_len > ret_right_len) {
       ret_add_sp_r = ret_left_len - ret_right_len + 1;
     }
     if (ret_right_len > ret_left_len) {
       ret_add_sp_l = ret_right_len - ret_left_len;
     }
     for (i = 0 ; i < ret_add_sp_l ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ');
     }
     save_append_string(&ret_str, &ret_space, o, "5' ");
     for (i = 0 ; i < ret_left_len ; i++) {
       save_append_char(&ret_str, &ret_space, o, (char) oligo1[i]);
     }
     save_append_string(&ret_str, &ret_space, o, "U+2510\\n   ");
     for (i = 0 ; i < ret_add_sp_l ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ');
     }
     for (i = 0 ; i < ret_left_len ; i++) {
       if (asciiRow[i] == '/') {
         save_append_char(&ret_str, &ret_space, o, '|');
       } else {
         save_append_char(&ret_str, &ret_space, o, ' ');
       }
     }
     if (ret_center_char == '|' ) {
       save_append_string(&ret_str, &ret_space, o, "U+2502");
     } else {
       save_append_char(&ret_str, &ret_space, o, ret_center_char);
     }
     save_append_string(&ret_str, &ret_space, o, "\\n");
     for (i = 0 ; i < ret_add_sp_r - 1 ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ');
     }
     save_append_string(&ret_str, &ret_space, o, "3' ");
     for (i = len1 ; i > ret_right_start - 1; i--) {
       save_append_char(&ret_str, &ret_space, o, (char) oligo1[i]);
     }
     save_append_string(&ret_str, &ret_space, o, "U+2518\\n");

/*
     save_append_string(&ret_str, &ret_space, o, "SEQ ");
     for(i = 0; i < len1; ++i) {
       save_append_char(&ret_str, &ret_space, o, asciiRow[i]);
     }
     save_append_string(&ret_str, &ret_space, o, "\\nSTR ");
     save_append_string(&ret_str, &ret_space, o, (const char*) oligo1);
     save_append_string(&ret_str, &ret_space, o, "\\n");
*/
     ret_ptr = (char *) safe_malloc(strlen(ret_str) + 1, o);
     strcpy(ret_ptr, ret_str);
     if (ret_str != NULL) {
       free(ret_str);
     }
   }
   free(asciiRow);
   return ret_ptr;
}

static void
save_append_string(char** ret, int *space, thal_results *o, const char *str) {
  int xlen, slen;
  if (str == NULL) {
    return;
  }
  if (*ret == NULL) {
    *ret = (char *) safe_malloc(sizeof(char)*500, o);
    *ret[0] = '\0';
    *space = 500;
  }
  xlen = strlen(*ret);
  slen = strlen(str);
  if (xlen + slen + 1 > *space) {
    *space += 4 * (slen + 1);
    *ret = (char *) safe_realloc(*ret, *space, o);
  }
  strcpy(*ret + xlen, str);
  return;
}

static void
save_append_char(char** ret, int *space, thal_results *o, const char str) {
  char fix[3];
  fix[0] = str;
  fix[1] = '\0';
  save_append_string(ret, space, o, fix);
}


static int
equal(double a, double b)
{
#ifdef INTEGER
   return a == b;
#endif

   if (!isfinite(a) || !isfinite(b))
     return 0;
   return fabs(a - b) < 1e-5;

   if (a == 0 && b == 0)
     return 1;
}

static void
strcatc(char* str, char c)
{
   str[strlen(str) + 1] = 0;
   str[strlen(str)] = c;
}

/* This file is created by thal_parameters_c_create.pl
   Do not edit this file, edit the script instead!
 */
static void * _thpr_safe_char_cp_malloc(const char *ct);

int set_default_thal_parameters(thal_parameters *a) {
  const char *dangle_dh = "0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n-500\n4700\n-4100\n-3800\n0\n"
            "0\n0\n0\n0\n0\n0\n0\n-5900\n-2600\n-3200\n-5200\n0\n0\n0\n0\n0\n0\n"
            "0\n0\n-2100\n-200\n-3900\n-4400\n0\n0\n0\n0\n0\n0\n0\n0\n-700\n4400\n-1600\n"
            "2900\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n"
            "0\n0\n0\n0\n0\n0\n0\n0\n-2900\n-4100\n-4200\n-200\n0\n0\n0\n0\n0\n"
            "0\n0\n0\n-3700\n-4000\n-3900\n-4900\n0\n0\n0\n0\n0\n0\n0\n0\n-6300\n-4400\n"
            "-5100\n-4000\n0\n0\n0\n0\n0\n0\n0\n0\n200\n600\n-1100\n-6900\n0\n0\n0\n"
            "0\n0\n0\n0\n0\n0\n0\n0\n0\n";

  const char *dangle_ds = "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-1.1\n14.2\n-13.1\n-12.6\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\n-16.5\n-7.4\n-10.4\n-15\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\n-3.9\n-0.1\n-11.2\n-13.1\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-0.8\n14.9\n-3.6\n"
            "10.4\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7.6\n-13\n-15\n-0.5\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\n-10\n-11.9\n-10.9\n-13.8\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-17.1\n-12.6\n"
            "-14\n-10.9\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n2.3\n3.3\n-1.6\n-20\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n";

  const char *loops_dh = "1\tinf\t0.0\tinf\n2\tinf\t0.0\tinf\n3\t0.0\t0.0\t0.0\n4\t0.0\t0.0\t0.0\n5\t0.0\t0.0\t0.0\n"
            "6\t0.0\t0.0\t0.0\n7\t0.0\t0.0\t0.0\n8\t0.0\t0.0\t0.0\n9\t0.0\t0.0\t0.0\n10\t0.0\t0.0\t0.0\n"
            "11\t0.0\t0.0\t0.0\n12\t0.0\t0.0\t0.0\n13\t0.0\t0.0\t0.0\n14\t0.0\t0.0\t0.0\n15\t0.0\t0.0\t0.0\n"
            "16\t0.0\t0.0\t0.0\n17\t0.0\t0.0\t0.0\n18\t0.0\t0.0\t0.0\n19\t0.0\t0.0\t0.0\n20\t0.0\t0.0\t0.0\n"
            "21\t0.0\t0.0\t0.0\n22\t0.0\t0.0\t0.0\n23\t0.0\t0.0\t0.0\n24\t0.0\t0.0\t0.0\n25\t0.0\t0.0\t0.0\n"
            "26\t0.0\t0.0\t0.0\n27\t0.0\t0.0\t0.0\n28\t0.0\t0.0\t0.0\n29\t0.0\t0.0\t0.0\n30\t0.0\t0.0\t0.0\n"
            "";

  const char *loops_ds = "1\t-1.0\t-12.89\t-1.0\n2\t-1.0\t-9.35\t-1.0\n3\t-10.31\t-9.99\t-11.28\n4\t-11.6\t-10.31\t-11.28\n5\t-12.89\t-10.64\t-10.64\n"
            "6\t-14.18\t-11.28\t-12.89\n7\t-14.83\t-11.92\t-13.54\n8\t-15.47\t-12.57\t-13.86\n9\t-15.79\t-13.21\t-14.5\n10\t-15.79\t-13.86\t-14.83\n"
            "11\t-16.26\t-14.32\t-15.29\n12\t-16.76\t-14.5\t-16.12\n13\t-17.15\t-14.89\t-16.5\n14\t-17.41\t-15.47\t-16.44\n15\t-17.74\t-15.81\t-16.77\n"
            "16\t-18.05\t-16.12\t-17.08\n17\t-18.34\t-16.41\t-17.38\n18\t-18.7\t-16.76\t-17.73\n19\t-18.96\t-17.02\t-17.99\n20\t-19.02\t-17.08\t-18.37\n"
            "21\t-19.25\t-17.32\t-18.61\n22\t-19.48\t-17.55\t-18.84\n23\t-19.7\t-17.76\t-19.05\n24\t-19.9\t-17.97\t-19.26\n25\t-20.31\t-18.05\t-19.66\n"
            "26\t-20.5\t-18.24\t-19.85\n27\t-20.68\t-18.42\t-20.04\n28\t-20.86\t-18.6\t-20.21\n29\t-21.03\t-18.77\t-20.38\n30\t-21.28\t-19.02\t-20.31\n"
            "";

  const char *stack_dh = "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7900\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-8400\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7800\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7200\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\n-8500\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\n-8000\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\n-10600\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\n-7800\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-8200\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-9800\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-8000\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-8400\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7200\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\n-8200\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\n-8500\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\n-7900\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\n";

  const char *stack_ds = "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.2\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.4\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-21.0\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-20.4\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.7\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\n-19.9\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\n-27.2\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\n-21.0\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.2\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-24.4\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-19.9\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.4\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-21.3\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\ninf\ninf\n-22.2\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\ninf\ninf\n-22.7\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\ninf\n-22.2\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n"
            "inf\n";

  const char *stackmm_dh = "inf\ninf\ninf\n4700\ninf\ninf\ninf\n7600\ninf\ninf\ninf\n3000\n1200\n2300\n-600\ninf\ninf\n"
            "inf\n-2900\ninf\ninf\ninf\n-700\ninf\ninf\ninf\n500\ninf\n5300\n-10\ninf\n700\ninf\n-900\n"
            "inf\ninf\ninf\n600\ninf\ninf\ninf\n-4000\ninf\ninf\n-700\ninf\n-3100\n1000\n1200\ninf\ninf\n"
            "inf\n5300\ninf\ninf\ninf\n-700\ninf\ninf\ninf\ninf\n-1200\n-2500\n-2700\ninf\ninf\ninf\n3400\n"
            "inf\ninf\ninf\n6100\n-900\n1900\n-700\ninf\ninf\ninf\ninf\n1000\ninf\ninf\n5200\ninf\ninf\n"
            "inf\n3600\ninf\n600\n-1500\ninf\n-800\ninf\ninf\n5200\ninf\ninf\n1900\ninf\ninf\ninf\n-1500\n"
            "inf\ninf\n-4000\ninf\n-4900\n-4100\ninf\n-1500\ninf\ninf\n2300\ninf\ninf\ninf\n-10\ninf\ninf\n"
            "inf\ninf\n-1500\n-2800\n-5000\n-1200\ninf\ninf\ninf\ninf\ninf\ninf\n700\n-2900\n5200\n-600\ninf\n"
            "inf\ninf\ninf\n1600\ninf\ninf\ninf\n-1300\ninf\ninf\n-600\ninf\n-700\n3600\ninf\n2300\ninf\n"
            "inf\n-6000\ninf\ninf\ninf\n-4400\ninf\ninf\n-700\ninf\ninf\n500\ninf\n-6000\n3300\ninf\n-4900\n"
            "inf\ninf\ninf\n-2800\ninf\n5800\n-600\ninf\ninf\ninf\ninf\n5200\n-4400\n-2200\n-3100\ninf\ninf\n"
            "inf\n-2500\ninf\n4100\ninf\ninf\n3400\n700\ninf\ninf\ninf\ninf\n1200\ninf\ninf\ninf\n-100\n"
            "inf\ninf\ninf\n200\n7600\n6100\ninf\n1200\ninf\ninf\n2300\ninf\ninf\ninf\n3300\ninf\ninf\n"
            "inf\n-2200\ninf\n3000\ninf\n1600\n-100\ninf\n-800\ninf\ninf\ninf\n-4100\ninf\n-1400\ninf\n-5000\n"
            "inf\ninf\ninf\n1000\n-1300\n200\n700\ninf\ninf\ninf\n1000\ninf\n5800\ninf\n-2700\ninf\ninf\n"
            "inf\n";

  const char *stackmm_ds = "inf\ninf\ninf\n12.9\ninf\ninf\ninf\n20.2\ninf\ninf\ninf\n7.4\n1.7\n4.6\n-2.3\ninf\ninf\n"
            "inf\n-9.8\ninf\ninf\ninf\n-3.8\ninf\ninf\ninf\n3.2\ninf\n14.6\n-4.4\ninf\n0.2\ninf\n-4.2\n"
            "inf\ninf\ninf\n-0.6\ninf\ninf\ninf\n-13.2\ninf\ninf\n-2.3\ninf\n-9.5\n0.9\n1.7\ninf\ninf\n"
            "inf\n14.6\ninf\ninf\ninf\n-2.3\ninf\ninf\ninf\ninf\n-6.2\n-8.3\n-10.8\ninf\ninf\ninf\n8.0\n"
            "inf\ninf\ninf\n16.4\n-4.2\n3.7\n-2.3\ninf\ninf\ninf\ninf\n0.7\ninf\ninf\n14.2\ninf\ninf\n"
            "inf\n8.9\ninf\n-0.6\n-7.2\ninf\n-4.5\ninf\ninf\n13.5\ninf\ninf\n3.7\ninf\ninf\ninf\n-7.2\n"
            "inf\ninf\n-13.2\ninf\n-15.3\n-11.7\ninf\n-6.1\ninf\ninf\n4.6\ninf\ninf\ninf\n-4.4\ninf\ninf\n"
            "inf\ninf\n-6.1\n-8.0\n-15.8\n-6.2\ninf\ninf\ninf\ninf\ninf\ninf\n0.7\n-9.8\n14.2\n-1.0\ninf\n"
            "inf\ninf\ninf\n3.6\ninf\ninf\ninf\n-5.3\ninf\ninf\n-1.0\ninf\n-3.8\n8.9\ninf\n5.4\ninf\n"
            "inf\n-15.8\ninf\ninf\ninf\n-12.3\ninf\ninf\n-2.3\ninf\ninf\n3.2\ninf\n-15.8\n10.4\ninf\n-15.3\n"
            "inf\ninf\ninf\n-8.0\ninf\n16.3\n-2.3\ninf\ninf\ninf\ninf\n13.5\n-12.3\n-8.4\n-9.5\ninf\ninf\n"
            "inf\n-8.3\ninf\n9.5\ninf\n12.9\n8.0\n0.7\ninf\ninf\ninf\ninf\n0.7\ninf\ninf\ninf\n-1.7\n"
            "inf\ninf\ninf\n-1.5\n20.2\n16.4\ninf\n0.7\ninf\ninf\n5.4\ninf\ninf\ninf\n10.4\ninf\ninf\n"
            "inf\n-8.4\ninf\n7.4\ninf\n3.6\n-1.7\ninf\n-4.5\ninf\ninf\ninf\n-11.7\ninf\n-6.2\ninf\n-15.8\n"
            "inf\ninf\ninf\n0.7\n-5.3\n-1.5\n0.2\ninf\ninf\ninf\n0.9\ninf\n16.3\ninf\n-10.8\ninf\ninf\n"
            "inf\n";

  const char *tetraloop_dh = "AAAAAT\t500\nAAAACT\t700\nAAACAT\t1000\nACTTGT\t0\nAGAAAT\t-1100\n"
            "AGAGAT\t-1100\nAGATAT\t-1500\nAGCAAT\t-1600\nAGCGAT\t-1100\nAGCTTT\t200\n"
            "AGGAAT\t-1100\nAGGGAT\t-1100\nAGGGGT\t500\nAGTAAT\t-1600\nAGTGAT\t-1100\n"
            "AGTTCT\t800\nATTCGT\t-200\nATTTGT\t0\nATTTTT\t-500\nCAAAAG\t500\n"
            "CAAACG\t700\nCAACAG\t1000\nCAACCG\t0\nCCTTGG\t0\nCGAAAG\t-1100\n"
            "CGAGAG\t-1100\nCGATAG\t-1500\nCGCAAG\t-1600\nCGCGAG\t-1100\nCGCTTG\t200\n"
            "CGGAAG\t-1100\nCGGGAG\t-1000\nCGGGGG\t500\nCGTAAG\t-1600\nCGTGAG\t-1100\n"
            "CGTTCG\t800\nCTTCGG\t-200\nCTTTGG\t0\nCTTTTG\t-500\nGAAAAC\t500\n"
            "GAAACC\t700\nGAACAC\t1000\nGCTTGC\t0\nGGAAAC\t-1100\nGGAGAC\t-1100\n"
            "GGATAC\t-1600\nGGCAAC\t-1600\nGGCGAC\t-1100\nGGCTTC\t200\nGGGAAC\t-1100\n"
            "GGGGAC\t-1100\nGGGGGC\t500\nGGTAAC\t-1600\nGGTGAC\t-1100\nGGTTCC\t800\n"
            "GTTCGC\t-200\nGTTTGC\t0\nGTTTTC\t-500\nTAAAAA\t500\nTAAACA\t700\n"
            "TAACAA\t1000\nTCTTGA\t0\nTGAAAA\t-1100\nTGAGAA\t-1100\nTGATAA\t-1600\n"
            "TGCAAA\t-1600\nTGCGAA\t-1100\nTGCTTA\t200\nTGGAAA\t-1100\nTGGGAA\t-1100\n"
            "TGGGGA\t500\nTGTAAA\t-1600\nTGTGAA\t-1100\nTGTTCA\t800\nTTTCGA\t-200\n"
            "TTTTGA\t0\nTTTTTA\t-500\n";

  const char *tetraloop_ds = "AAAAAT\t-650\nAAAACT\t1610\nAAACAT\t1610\nACTTGT\t4190\nAGAAAT\t1610\n"
            "AGAGAT\t1610\nAGATAT\t1610\nAGCAAT\t1610\nAGCGAT\t1610\nAGCTTT\t1610\n"
            "AGGAAT\t1610\nAGGGAT\t1610\nAGGGGT\t640\nAGTAAT\t1610\nAGTGAT\t1610\n"
            "AGTTCT\t1610\nATTCGT\t1610\nATTTGT\t1610\nATTTTT\t1610\nCAAAAG\t-1290\n"
            "CAAACG\t0\nCAACAG\t0\nCAACCG\t0\nCCTTGG\t2570\nCGAAAG\t0\n"
            "CGAGAG\t0\nCGATAG\t0\nCGCAAG\t0\nCGCGAG\t0\nCGCTTG\t0\n"
            "CGGAAG\t0\nCGGGAG\t0\nCGGGGG\t-970\nCGTAAG\t0\nCGTGAG\t0\n"
            "CGTTCG\t0\nCTTCGG\t0\nCTTTGG\t0\nCTTTTG\t0\nGAAAAC\t-3230\n"
            "GAAACC\t0\nGAACAC\t0\nGCTTGC\t2570\nGGAAAC\t0\nGGAGAC\t0\n"
            "GGATAC\t0\nGGCAAC\t0\nGGCGAC\t0\nGGCTTC\t0\nGGGAAC\t0\n"
            "GGGGAC\t0\nGGGGGC\t-970\nGGTAAC\t0\nGGTGAC\t0\nGGTTCC\t0\n"
            "GTTCGC\t0\nGTTTGC\t0\nGTTTTC\t0\nTAAAAA\t320\nTAAACA\t1610\n"
            "TAACAA\t1610\nTCTTGA\t4190\nTGAAAA\t1610\nTGAGAA\t1610\nTGATAA\t1610\n"
            "TGCAAA\t1610\nTGCGAA\t1610\nTGCTTA\t1610\nTGGAAA\t1610\nTGGGAA\t1610\n"
            "TGGGGA\t640\nTGTAAA\t1610\nTGTGAA\t1610\nTGTTCA\t1610\nTTTCGA\t1610\n"
            "TTTTGA\t1610\nTTTTTA\t1610\n";

  const char *triloop_dh = "AGAAT\t-1500\nAGCAT\t-1500\nAGGAT\t-1500\nAGTAT\t-1500\nCGAAG\t-2000\n"
            "CGCAG\t-2000\nCGGAG\t-2000\nCGTAG\t-2000\nGGAAC\t-2000\nGGCAC\t-2000\n"
            "GGGAC\t-2000\nGGTAC\t-2000\nTGAAA\t-1500\nTGCAA\t-1500\nTGGAA\t-1500\n"
            "TGTAA\t-1500\n";

  const char *triloop_ds = "AGAAT\t0\nAGCAT\t0\nAGGAT\t0\nAGTAT\t0\nCGAAG\t0\n"
            "CGCAG\t0\nCGGAG\t0\nCGTAG\t0\nGGAAC\t0\nGGCAC\t0\n"
            "GGGAC\t0\nGGTAC\t0\nTGAAA\t0\nTGCAA\t0\nTGGAA\t0\n"
            "TGTAA\t0\n";

  const char *tstack_dh = "0\n0\n0\n-2500\n0\n0\n0\n-2700\n0\n0\n0\n-2400\n-3100\n-1600\n-1900\n0\n0\n"
            "0\n-8000\n0\n0\n0\n-3200\n0\n0\n0\n-4600\n0\n-1800\n-100\n0\n-900\n0\n-4300\n"
            "0\n0\n0\n-2700\n0\n0\n0\n-6000\n0\n0\n-2500\n0\n-1100\n-3200\n-3100\n0\n0\n"
            "0\n-1800\n0\n0\n0\n-2500\n0\n0\n0\n0\n-2300\n-3500\n-2400\n0\n0\n0\n-2300\n"
            "0\n0\n0\n-700\n-4300\n-2600\n-3900\n0\n0\n0\n0\n-700\n0\n0\n-5000\n0\n0\n"
            "0\n-3900\n0\n-2700\n-2100\n0\n-3200\n0\n0\n-3000\n0\n0\n-2600\n0\n0\n0\n-2100\n"
            "0\n0\n-6000\n0\n-3800\n-3800\n0\n-3900\n0\n0\n-1600\n0\n0\n0\n-100\n0\n0\n"
            "0\n0\n-3900\n-6600\n-6100\n-2300\n0\n0\n0\n0\n0\n0\n-2000\n-8000\n-5000\n-4300\n0\n"
            "0\n0\n0\n-1100\n0\n0\n0\n-3600\n0\n0\n-4300\n0\n-3200\n-3900\n0\n-4900\n0\n"
            "0\n-700\n0\n0\n0\n-5900\n0\n0\n-3900\n0\n0\n-4600\n0\n-700\n-5700\n0\n-3800\n"
            "0\n0\n0\n-6600\n0\n0\n-1900\n0\n0\n0\n0\n-3000\n-5900\n-7400\n-1100\n0\n0\n"
            "0\n-3500\n0\n0\n0\n-2500\n-2300\n-2000\n-7200\n0\n0\n0\n-2500\n0\n0\n0\n-3900\n"
            "0\n0\n0\n-3200\n-2700\n-700\n0\n-2500\n0\n0\n-4900\n0\n0\n0\n-5700\n0\n0\n"
            "0\n-7400\n0\n-2400\n0\n-1100\n-3900\n0\n-3200\n0\n0\n0\n-3800\n0\n0\n0\n-6100\n"
            "0\n0\n0\n-700\n-3600\n-3200\n-900\n0\n0\n0\n-3200\n0\n0\n0\n-2400\n0\n0\n"
            "0\n";

  const char *tstack2_dh = "0\n0\n0\n-2500\n0\n0\n0\n-2700\n0\n0\n0\n-2400\n-3100\n-1600\n-1900\n-5000\n0\n"
            "0\n-8000\n0\n0\n0\n-3200\n0\n0\n0\n-4600\n0\n-1800\n-100\n-6000\n-900\n0\n-4300\n"
            "0\n0\n0\n-2700\n0\n0\n0\n-6000\n0\n0\n-2500\n-6000\n-1100\n-3200\n-3100\n0\n0\n"
            "0\n-1800\n0\n0\n0\n-2500\n0\n0\n0\n-5000\n-2300\n-3500\n-2400\n0\n0\n0\n-2300\n"
            "0\n0\n0\n-700\n-4300\n-2600\n-3900\n-6000\n0\n0\n0\n-700\n0\n0\n-5000\n0\n0\n"
            "0\n-3900\n0\n-2700\n-2100\n-7000\n-3200\n0\n0\n-3000\n0\n0\n-2600\n0\n0\n0\n-2100\n"
            "0\n0\n-6000\n-7000\n-3800\n-3800\n0\n-3900\n0\n0\n-1600\n0\n0\n0\n-100\n0\n0\n"
            "0\n-6000\n-3900\n-6600\n-6100\n-2300\n0\n0\n0\n0\n0\n0\n-2000\n-8000\n-5000\n-4300\n-6000\n"
            "0\n0\n0\n-1100\n0\n0\n0\n-3600\n0\n0\n-4300\n0\n-3200\n-3900\n-7000\n-4900\n0\n"
            "0\n-700\n0\n0\n0\n-5900\n0\n0\n-3900\n0\n0\n-4600\n-7000\n-700\n-5700\n0\n-3800\n"
            "0\n0\n0\n-6600\n0\n0\n-1900\n0\n0\n0\n-6000\n-3000\n-5900\n-7400\n-1100\n0\n0\n"
            "0\n-3500\n0\n0\n0\n-2500\n-2300\n-2000\n-5000\n0\n0\n0\n-2500\n0\n0\n0\n-3900\n"
            "0\n0\n0\n-3200\n-2700\n-700\n-6000\n-2500\n0\n0\n-4900\n0\n0\n0\n-5700\n0\n0\n"
            "0\n-7400\n0\n-2400\n-6000\n-1100\n-3900\n0\n-3200\n0\n0\n0\n-3800\n0\n0\n0\n-6100\n"
            "0\n0\n-5000\n-700\n-3600\n-3200\n-900\n0\n0\n0\n-3200\n0\n0\n0\n-2400\n0\n0\n"
            "0\n";

  const char *tstack2_ds = "inf\ninf\ninf\n-6.3\ninf\ninf\ninf\n-7.0\ninf\ninf\ninf\n-5.8\n-7.8\n-4.0\n-4.4\n-13.5\ninf\n"
            "inf\n-22.5\ninf\ninf\ninf\n-7.1\ninf\ninf\ninf\n-11.4\ninf\n-3.8\n-0.5\n-16.1\n-1.7\ninf\n-10.7\n"
            "inf\ninf\ninf\n-6.0\ninf\ninf\ninf\n-15.5\ninf\ninf\n-5.9\n-16.1\n-2.1\n-8.7\n-7.8\ninf\ninf\n"
            "inf\n-3.8\ninf\ninf\ninf\n-5.9\ninf\ninf\ninf\n-13.6\n-6.3\n-9.4\n-6.5\ninf\ninf\ninf\n-5.9\n"
            "inf\ninf\ninf\n-1.3\n-10.7\n-5.9\n-9.6\n-16.1\ninf\ninf\ninf\n-1.2\ninf\ninf\n-13.8\ninf\ninf\n"
            "inf\n-10.6\ninf\n-6.0\n-5.1\n-19.3\n-8.0\ninf\ninf\n-7.8\ninf\ninf\n-5.9\ninf\ninf\ninf\n-5.1\n"
            "inf\ninf\n-15.5\n-19.3\n-9.5\n-9.0\ninf\n-10.6\ninf\ninf\n-4.0\ninf\ninf\ninf\n-0.5\ninf\ninf\n"
            "inf\n-16.1\n-10.6\n-18.7\n-16.9\n-6.3\ninf\ninf\ninf\ninf\ninf\ninf\n-4.7\n-22.5\n-13.8\n-11.1\n-16.1\n"
            "inf\ninf\ninf\n-2.7\ninf\ninf\ninf\n-9.8\ninf\ninf\n-11.1\ninf\n-7.1\n-10.6\n-19.3\n-13.5\ninf\n"
            "inf\n-19.2\ninf\ninf\ninf\n-16.1\ninf\ninf\n-9.6\ninf\ninf\n-11.4\n-19.3\n-19.2\n-15.9\ninf\n-9.5\n"
            "inf\ninf\ninf\n-18.7\ninf\ninf\n-4.4\ninf\ninf\ninf\n-16.1\n-7.8\n-16.1\n-21.2\n-2.1\ninf\ninf\n"
            "inf\n-9.4\ninf\ninf\ninf\n-6.3\n-5.9\n-4.7\n-14.2\ninf\ninf\ninf\n-6.3\ninf\ninf\ninf\n-10.5\n"
            "inf\ninf\ninf\n-8.9\n-7.0\n-1.3\n-16.1\n-6.3\ninf\ninf\n-13.5\ninf\ninf\ninf\n-15.9\ninf\ninf\n"
            "inf\n-21.2\ninf\n-5.8\n-16.1\n-2.7\n-10.5\ninf\n-8.0\ninf\ninf\ninf\n-9.0\ninf\ninf\ninf\n-16.9\n"
            "inf\ninf\n-13.5\n-1.2\n-9.8\n-8.9\n-1.7\ninf\ninf\ninf\n-8.7\ninf\ninf\ninf\n-6.5\ninf\ninf\n"
            "inf\n";

  const char *tstack_tm_inf_ds = "inf\ninf\ninf\n-6.3\ninf\ninf\ninf\n-7.0\ninf\ninf\ninf\n-5.8\n-7.8\n-4.0\n-4.4\ninf\ninf\n"
            "inf\n-22.5\ninf\ninf\ninf\n-7.1\ninf\ninf\ninf\n-11.4\ninf\n-3.8\n-0.5\ninf\n-1.7\ninf\n-10.7\n"
            "inf\ninf\ninf\n-6.0\ninf\ninf\ninf\n-15.5\ninf\ninf\n-5.9\ninf\n-2.1\n-8.7\n-7.8\ninf\ninf\n"
            "inf\n-3.8\ninf\ninf\ninf\n-5.9\ninf\ninf\ninf\ninf\n-6.3\n-9.4\n-6.5\ninf\ninf\ninf\n-5.9\n"
            "inf\ninf\ninf\n-1.3\n-10.7\n-5.9\n-9.6\ninf\ninf\ninf\ninf\n-1.2\ninf\ninf\n-13.8\ninf\ninf\n"
            "inf\n-10.6\ninf\n-6.0\n-5.1\ninf\n-8.0\ninf\ninf\n-7.8\ninf\ninf\n-5.9\ninf\ninf\ninf\n-5.1\n"
            "inf\ninf\n-15.5\ninf\n-9.5\n-9.0\ninf\n-10.6\ninf\ninf\n-4.0\ninf\ninf\ninf\n-0.5\ninf\ninf\n"
            "inf\ninf\n-10.6\n-18.7\n-16.9\n-6.3\ninf\ninf\ninf\ninf\ninf\ninf\n-4.7\n-22.5\n-13.8\n-11.1\ninf\n"
            "inf\ninf\ninf\n-2.7\ninf\ninf\ninf\n-9.8\ninf\ninf\n-11.1\ninf\n-7.1\n-10.6\ninf\n-13.5\ninf\n"
            "inf\n-19.2\ninf\ninf\ninf\n-16.1\ninf\ninf\n-9.6\ninf\ninf\n-11.4\ninf\n-19.2\n-15.9\ninf\n-9.5\n"
            "inf\ninf\ninf\n-18.7\ninf\ninf\n-4.4\ninf\ninf\ninf\ninf\n-7.8\n-16.1\n-21.2\n-2.1\ninf\ninf\n"
            "inf\n-9.4\ninf\ninf\ninf\n-6.3\n-5.9\n-4.7\ninf\ninf\ninf\ninf\n-6.3\ninf\ninf\ninf\n-10.5\n"
            "inf\ninf\ninf\n-8.9\n-7.0\n-1.3\ninf\n-6.3\ninf\ninf\n-13.5\ninf\ninf\ninf\n-15.9\ninf\ninf\n"
            "inf\n-21.2\ninf\n-5.8\ninf\n-2.7\n-10.5\ninf\n-8.0\ninf\ninf\ninf\n-9.0\ninf\ninf\ninf\n-16.9\n"
            "inf\ninf\ninf\n-1.2\n-9.8\n-8.9\n-1.7\ninf\ninf\ninf\n-8.7\ninf\ninf\ninf\n-6.5\ninf\ninf\n"
            "inf\n";

  thal_free_parameters(a);

  a->dangle_dh = _thpr_safe_char_cp_malloc(dangle_dh);
  a->dangle_ds = _thpr_safe_char_cp_malloc(dangle_ds);
  a->loops_dh = _thpr_safe_char_cp_malloc(loops_dh);
  a->loops_ds = _thpr_safe_char_cp_malloc(loops_ds);
  a->stack_dh = _thpr_safe_char_cp_malloc(stack_dh);
  a->stack_ds = _thpr_safe_char_cp_malloc(stack_ds);
  a->stackmm_dh = _thpr_safe_char_cp_malloc(stackmm_dh);
  a->stackmm_ds = _thpr_safe_char_cp_malloc(stackmm_ds);
  a->tetraloop_dh = _thpr_safe_char_cp_malloc(tetraloop_dh);
  a->tetraloop_ds = _thpr_safe_char_cp_malloc(tetraloop_ds);
  a->triloop_dh = _thpr_safe_char_cp_malloc(triloop_dh);
  a->triloop_ds = _thpr_safe_char_cp_malloc(triloop_ds);
  a->tstack_tm_inf_ds = _thpr_safe_char_cp_malloc(tstack_tm_inf_ds);
  a->tstack_dh = _thpr_safe_char_cp_malloc(tstack_dh);
  a->tstack2_dh = _thpr_safe_char_cp_malloc(tstack2_dh);
  a->tstack2_ds = _thpr_safe_char_cp_malloc(tstack2_ds);

  return 0;
};

static void * _thpr_safe_char_cp_malloc(const char *ct) {
  void *r = malloc((strlen(ct) + 1) * sizeof(char));
  if (NULL == r) {
    fprintf(stderr, "out of memory in thal_parameters\n");
    exit(-2);
  }
  strcpy(r, ct);
  return r;
}

/* Beginning of main 
'mv', 'dv', 'maxLoop', 'dntp', 'dna_conc']
*/
void align(char* oligo1, char* oligo2, int type, double mv, double dv, int maxLoop, double dntp, double dna_conc, double celsius)
{
   const char* usage;
   int tmp_ret = 0;
   thal_results o;
   thal_args a;
   set_thal_default_args(&a);

   a.type = type;
   a.mv = mv;
   a.dv = dv;
   a.maxLoop = maxLoop;
   a.dntp = dntp;
   a.dna_conc = dna_conc;
   a.temp = (double) celsius + ABSOLUTE_ZERO;
   
   thal_mode mode = THL_FAST; /* by default print only melting temperature,
                                    do not draw structure or print any additional parameters */
   int thal_debug = 0;
   int thal_only = 0;
   
   /* read default thermodynamic parameters */
   thal_parameters thermodynamic_parameters;
   thal_set_null_parameters(&thermodynamic_parameters);
   set_default_thal_parameters(&thermodynamic_parameters);
   get_thermodynamic_values(&thermodynamic_parameters, &o);

   /* alignment */
   thal(oligo1,oligo2,&a,mode,&o);

   /* encountered error during thermodynamical calc */
   if (o.temp == THAL_ERROR_SCORE) {
      tmp_ret = fprintf(stderr, "Error: %s\n", o.msg);
      exit(-1);
   }

   /* cleanup */
   destroy_thal_structures();
   thal_free_parameters(&thermodynamic_parameters);
   
   Inline_Stack_Vars;
 
  Inline_Stack_Reset;
  Inline_Stack_Push(sv_2mortal(newSVnv(o.temp)));
  Inline_Stack_Push(sv_2mortal(newSVnv(o.H)));
  Inline_Stack_Push(sv_2mortal(newSVnv(o.S)));
  Inline_Stack_Push(sv_2mortal(newSVnv(o.dG)));

  Inline_Stack_Push(sv_2mortal(newSVpv(o.duplex0, strlen(o.duplex0))));
  Inline_Stack_Push(sv_2mortal(newSVpv(o.duplex1, strlen(o.duplex1))));
  Inline_Stack_Push(sv_2mortal(newSVpv(o.duplex2, strlen(o.duplex2))));
  Inline_Stack_Push(sv_2mortal(newSVpv(o.duplex3, strlen(o.duplex3))));

  Inline_Stack_Done;
  
  /*  return o.temp; */
}


