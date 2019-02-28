use lib 'inc';
use strict;
use warnings;
use ExtUtils::MakeMaker;
use Inline::Module;
#
WriteMakefile(
  NAME               => 'Ntthal::XS',
  AUTHOR             => q{Frank Reinecke <frank.reinecke@qiagen.com>},
  VERSION_FROM       => 'lib/Ntthal/XS.pm',
  ABSTRACT_FROM      => 'lib/Ntthal/XS.pm',
  LIBS               => [],
  EXE_FILES          => ['bin/ntthal_xs.pl'],
  LICENSE            => 'GPL',
  PL_FILES           => {},
  MIN_PERL_VERSION   => '5.006',
  CONFIGURE_REQUIRES => {
    'ExtUtils::MakeMaker' => '0',
  },
  BUILD_REQUIRES => {
    'Test::More' => '0',
    'Test::Warn' => '0',
  },
  PREREQ_PM => {
    'Data::Dump::Color'  => '0',
    'Mouse'              => '0',
    'Method::Signatures' => '0',
    'Inline::C'          => '0',
    'Inline::Module'     => '0',
    'Scalar::Util'       => '0',
  },
  dist      => {COMPRESS => 'gzip -9f',     SUFFIX => 'gz',},
  clean     => {FILES    => ['Ntthal-XS-*', '_Inline']},
  postamble => {
    inline => {
      module => 'Ntthal::XS',
      stub   => 'Ntthal::XS::Inline',
      ilsm   => 'Inline::C',
    },
  },
);
