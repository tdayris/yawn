#!/usr/bin/perl

use strict;
use warnings;
my ($vcf_path) = @ARGV;

open(my $vcf, "<:encoding(UTF-8)", $vcf_path) or die "Can't open $vcf_path: $!";

sub RemoveFieldList {
  my $cline = shift @_;

  foreach my $field (@_) {
    $cline = $cline =~ s/$field=[^;\t]+;?//gr;
  }

  $cline = $cline =~ s/;;/;/gr;
  $cline = $cline =~ s/;\t/\t/gr;

  print "$cline\n";
}

while ( my $line = <$vcf> ) {
  chomp $line;

  if ( $line =~ /^#/ ) {
    print "$line\n";
  } else {
    RemoveFieldList(
      $line,
      "ANN", "mir_ACC", "GeneNameRefSeq", "Statut",
      "IG_hom-het", "1000G_AF", "EVS_AF", "H3_AF",
      "FREQ_HOM", "FREQ_HTZ", "C_used", "A_used",
      "G_used", "T_used", "ref_upstream", "ref/indel",
      "ref_downstream", "max_gtype", "Qmax_gtype",
      "alt_reads", "indel_reads", "other_reads",
      "repeat_unit", "ref_repeat_count", "EVS",
      "indel_repeat_count", "EVS_CA","EVS_MAF", "ID",
      "Qdepth"
    );
  }
}
