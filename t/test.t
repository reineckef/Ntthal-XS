use Test::More tests => 44;
use Test::Warn;

BEGIN {
  use_ok('Ntthal::XS') || print "Bail out!\n";
}
diag("Testing Ntthal::XS $Ntthal::XS::VERSION, Perl $], $^X");
require Ntthal::XS;
my $Ntthal = new_ok('Ntthal::XS');

# set values back to defaults
my %Defaults = (
  only     => 'hash',
  type     => 'ANY',
  mv       => 50,
  dv       => 0,
  maxLoop  => 8,
  dntp     => 0.8,
  dna_conc => 50,
  temp     => 37
);
foreach my $arg (qw<mv dv maxLoop dntp dna_conc temp>) {
  my $value = sprintf "%.2f", rand(20);
  $Ntthal->$arg($value);
  is( $Ntthal->$arg, $value, "Failed to set $arg to $value" );
  $Ntthal->$arg( $Defaults{$arg} );    # return to default
}

# helpers
diag("Testing helper functions");
ok( Ntthal::XS::cleanseq("123 TTCAGTAggcATCAtcACGxxTA \n") eq 'TTCAGTAGGCATCATCACGNNTA', "clean sequence");
ok( Ntthal::XS::revcomp("TTCAGTAGGCATCATCACGTA") eq 'TACGTGATGATGCCTACTGAA',             "revcomp sequence");
ok( $Ntthal->args('command') eq "-a ANY -mv 50 -dv 0 -maxloop 8 -n 0.8 -d 50 -t 37",     "cmd argument list");

# run tm convenience method
my $tm = $Ntthal->tm(qw<TTCAGTAGGCATCATCACGTA>);
diag("Testing tm method");
ok( $tm > 50, "object t > 50");
ok( $tm < 51, "object t < 51");

# run "full" ntthal
my $result = $Ntthal->ntthal(qw<TTCAGTAGGCATCATCACGTA TACGTTATGATGCCTTACTGAA>);
diag("Testing hash return");
ok( $result->{t} > 33.84, "object t > 33.84");
ok( $result->{t} < 33.9,  "object t < 33.9");

# check temperature ranges
$Ntthal->only('tm');    # switch mode
$tm = $Ntthal->ntthal(qw<TTCAGTAGGCATCATCACGTA TACGTTATGATGCCTTACTGAA>);
diag("Testing only switch");
ok( $tm > 33.84, "raw tm $tm > 33.84");
ok( $tm < 33.9,  "raw tm $tm < 33.9");
diag("Testing salt effects");

# more salt
$Ntthal->mv(100);
my $tm_mv = $Ntthal->ntthal(qw<TTCAGTAGGCATCATCACGTA TACGTTATGATGCCTTACTGAA>);
ok( $tm_mv > $tm, "more mv should increase tm");

# even more salt
$Ntthal->dv(2);
my $tm_dv = $Ntthal->ntthal(qw<TTCAGTAGGCATCATCACGTA TACGTTATGATGCCTTACTGAA>);
ok( $tm_dv > $tm_mv, "more dv should increase tm");

# roll back
map { $Ntthal->$_( $Defaults{$_} ) } keys %Defaults;

# test modes
#
#  ntthal -a ANY -mv 50 -dv 0 -maxloop 8 -n 0.8 -d 50 -t 37 -s1 AATCAGTAGGCATCATCACGTA -s2 CCGACGTGATGATGCCTACTGAA
#  Calculated thermodynamical parameters for dimer:	dS = -425.344	dH = -150200	dG = -18279.6	t = 52.3063
#  SEQ	AA                   A--
#  SEQ	  TCAGTAGGCATCATCACGT
#  STR	  AGTCATCCGTAGTAGTGCA
#  STR	 A                   GCC
diag("Testing warning for unsupported mode");
$Ntthal->type('INVALID');
warning_like { $Ntthal->tm('AATCAGTAGGCATCATCACGTA') } qr/undefined/;

diag("Testing alignment mode ANY");
$Ntthal->type('ANY');
my $any = $Ntthal->ntthal(qw<AATCAGTAGGCATCATCACGTA CCGACGTGATGATGCCTACTGAA>);
ok( $any->{t} > 52,      "mode ANY tm $any->{t} is in expected range");
ok( $any->{t} < 53,      "mode ANY tm is in expected range");
ok( $any->{dH} < -150000, "mode ANY reports reasonable H");
ok( $any->{dS} < -425,    "mode ANY reports  reasonable S");
ok( $any->{dG} < -18270, "mode ANY reports stable dG");
#
#
#  ntthal -a END1 -mv 50 -dv 0 -maxloop 8 -n 0.8 -d 50 -t 37 -s1 AATCAGTAGGCATCATCACGTA -s2 CCGACGTGATGATGCCTACTGAA
#  Calculated thermodynamical parameters for dimer:	dS = -326.046	dH = -110100	dG = -8976.7	t = 30.8185
#  SEQ	AA               ACGT ------
#  SEQ	  TCAGTAGGCATCATC    A
#  STR	  AGTCATCCGTAGTAG    T
#  STR	 A               ---- GCAGCC
diag("Testing alignment mode END1");
$Ntthal->type('END1');
my $end1 = $Ntthal->ntthal(qw<AATCAGTAGGCATCATCACGTA CCGACGTGATGATGCCTACTGAA>);
ok( $end1->{t} < $any->{t}, "mode END1 is less stable");
ok( $end1->{t} > 30,        "mode END1 tm is in expected range");
ok( $end1->{t} < 31,        "mode END1 tm is in expected range");
ok( $end1->{dH} < -110000,   "mode END1 reports reasonable H");
ok( $end1->{dS} < -326,      "mode END1 reports  reasonable S");
ok( $end1->{dG} < -8976,    "mode END1 reports stable dG");

#  ntthal -a END2 -mv 50 -dv 0 -maxloop 8 -n 0.8 -d 50 -t 37 -s1 AATCAGTAGGCATCATCACGTA -s2 CCGACGTGATGATGCCTACTGAA
#  Calculated thermodynamical parameters for dimer:	dS = -437.134	dH = -149700	dG = -14123	t = 43.1426
#  SEQ	CCG                  A --
#  SEQ	   ACGTGATGATGCCTACTG A
#  STR	   TGCACTACTACGGATGAC T
#  STR	  A                  - AA
diag("Testing alignment mode END2");
$Ntthal->type('END2');
my $end2 = $Ntthal->ntthal(qw<AATCAGTAGGCATCATCACGTA CCGACGTGATGATGCCTACTGAA>);
ok( $end2->{t} < $any->{t},  "mode END2 is less stable than ANY");
ok( $end2->{t} > $end1->{t}, "mode END2 is more stable than END1");
ok( $end2->{t} > 43,         "mode END2 tm is in expected range");
ok( $end2->{t} < 44,         "mode END2 tm is in expected range");
ok( $end2->{dH} < -149000,    "mode END2 reports reasonable H");
ok( $end2->{dS} < -437,       "mode END2 reports reasonable S");
ok( $end2->{dG} < -14100,    "mode END2 reports stable dG");

# TTCAGTAGGCATCATCACGTAATGCCTACTGAA
#  ntthal -a HAIRPIN -mv 50 -dv 0 -maxloop 8 -n 0.8 -d 50 -t 37 -s1 TTCAGTAGGCATCATCACGTAATGCCTACTGAA
#  Calculated thermodynamical parameters for dimer:	33	dS = -257.124	dH = -86900	dG = -7152.9	t = 64.8188
#  SEQ	-///////////---------\\\\\\\\\\\-
#  STR	TTCAGTAGGCATCATCACGTAATGCCTACTGAA
diag("Testing alignment mode HAIRPIN");
$Ntthal->type('HAIRPIN');
my $hp = $Ntthal->ntthal(qw<TTCAGTAGGCATCATCACGTAATGCCTACTGAA>);
ok( $hp->{t} > 64,       "mode HAIRPIN tm is in expected range");
ok( $hp->{t} < 65,       "mode HAIRPIN tm is in expected range");
ok( $end2->{dH} < -86800, "mode HAIRPIN reports reasonable H");
ok( $end2->{dS} < -257,   "mode HAIRPIN reports  reasonable S");
ok( $end2->{dG} < -7150, "mode HAIRPIN reports stable dG");
#
diag("Testing counter");
ok( $Ntthal->counter == 10 );
