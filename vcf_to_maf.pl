#!/usr/bin/perl

use strict;
use warnings;
use Text::CSV;

sub parse_vcf_header {
  my $header = @_;
  if ($header =~ /^##\w+=<.*>$/) {
    print("Good header: $header");
  } else {
    print("Bad header: $header");
  }
}

sub parse_vcf_entry {
  print("ok");
}

my ($vcf_input, $maf_output) = @ARGV;

my $csv = Text::CSV->new ({ sep_char => "\t" });

open( my $vcf, "<:encoding(UTF-8)", $vcf_input ) or die "Can't open $vcf_input: $!";

while ( my $line = <$vcf> ) {
  chomp $line;
  next unless $line;
  die $csv->error_diag
    unless $csv->parse( $line )
    ;

  if ($line =~ "^##") {
      parse_vcf_header($line);
  } else {
      parse_vcf_entry($line);
  }
}
