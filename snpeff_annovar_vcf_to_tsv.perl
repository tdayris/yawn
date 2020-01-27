#!/usr/bin/perl

use strict;
use warnings;
use Text::CSV;

my ($annovar_txt) = @ARGV;

my $csv = Text::CSV->new ({ sep_char => "\t" });


print("Chr\tStart\tEnd\tGene_ID\tInheritance\tGnomAD_NFE\tExac_NFE\tGrantham_Score\tREF\tALT\tRegion\tGranham_Annotation\tCosmic70\tSift\tPolyphen\tMutationTaster\tCADD\tAA_Change\tFather_pval\tMother_pval\tChild_pval\tFather_freq\tMother_freq\tChild_freq\n");

open(my $annovar, "<:encoding(UTF-8)", $annovar_txt) or die "Can't open $annovar_txt: $!";


while ( my $line = <$annovar> ) {
    chomp $line;
    next unless $line;
    die $csv->error_diag
        unless $csv->parse( $line )
        ;

    next unless defined( my $chr = ( $csv->fields )[0] );
    next unless defined( my $start = ( $csv->fields )[1] );
    next unless defined( my $end = ( $csv->fields )[2] );

    next unless defined( my $ref = ( $csv->fields )[3] );
    next unless defined( my $alt = ( $csv->fields )[4] );

    next unless defined( my $cosmic70 = ( $csv->fields )[10] );
    next unless defined( my $sift = ( $csv->fields )[13] );
    next unless defined( my $polyphen = ( $csv->fields )[16] );
    next unless defined( my $mutation_taster = ( $csv->fields )[25] );
    next unless defined( my $cadd = ( $csv->fields )[48] );

    next unless defined( my $gene_func = ( $csv->fields )[77] );
    next unless defined( my $gene_name = ( $csv->fields )[129] );
    next unless defined( my $gene_id = ( $csv->fields )[78] );
    next unless defined( my $aa_change = ( $csv->fields )[81] );

    next unless defined( my $exac_nfe = ( $csv->fields )[87] );
    next unless defined( my $gnomad_nfe = ( $csv->fields )[96] );

    next unless defined( my $other = ( $csv->fields )[143] );
    next unless defined( my $father = ( $csv->fields )[145] );
    next unless defined( my $mother = ( $csv->fields )[146] );
    next unless defined( my $child = ( $csv->fields )[147] );

    my $not_chrm = ($chr =~ /^chr[1-9XY]+$/);
    my $exonic = ($gene_func =~ "exonic" and $gene_func !~ "^ncRNA");
    my $gnomad_threshold = ($gnomad_nfe =~ /^\.$/ or $gnomad_nfe <= 0.001);
    my $exac_threshold = ($exac_nfe =~ /^\.$/ or $exac_nfe <= 0.001);

    my ($inheritance) = $other =~ /Mutation_Inheritance=([^;]+);/;
    my $is_inherited = ($inheritance !~ /^Nobody_has_mutation_under_confidence_interval$/);

    my ($gran_score) = $other =~ /Grantham_score=([^;]+);/;
    my ($gran_annot) = $other =~ /Grantham_annotation=([^;]+);/;

    my @father_info = split /:/, $father;
    my @mother_info = split /:/, $mother;
    my @child_info = split /:/, $child;

    my $father_no_doubt = $father_info[7] !~ /E-/;
    my $mother_no_doubt = $mother_info[7] !~ /E-/;
    my $de_novo = $inheritance =~ /De_Novo/;

    my $father_freq = (split /%/, $father_info[6]) <= 40;
    my $mother_freq = (split /%/, $mother_info[6]) <= 40;
    my $child_freq = (split /%/, $child_info[6]) >= 40 & (split /%/, $child_info[6]) <= 60;

    my $cadd_threshold = ($cadd !~ /^\.$/ and $cadd >= 20);

    # Add allelic freq

    my $localization = $not_chrm & $exonic;
    my $population = $gnomad_threshold & $exac_threshold;
    my $strict = $father_no_doubt & $mother_no_doubt & $cadd_threshold;
    my $transmission = $de_novo | ($father_freq & $mother_freq & $child_freq)

    if ($localization & $population & $is_inherited & $strict & $transmission) {
      print(
        $chr, "\t", $start, "\t", $end, "\t", $gene_name, "\t", $inheritance, "\t",
        $gnomad_nfe, "\t", $exac_nfe, "\t", $gran_score, "\t", $ref, "\t", $alt, "\t",
        $gene_func, "\t", $cosmic70, "\t", $sift, "\t", $polyphen, "\t",
        $mutation_taster, "\t", $gran_annot, "\t", $cadd, "\t", $aa_change, "\t",
        $father_info[7], "\t", $mother_info[7], "\t", $child_info[7], "\t",
        $father_info[6], "\t", $mother_info[6], "\t", $child_info[6], "\n"
      );
    }
}
close $annovar_txt;
