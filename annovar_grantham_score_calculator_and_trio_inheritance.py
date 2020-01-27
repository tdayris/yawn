#!/usr/bin/python3.7
# conding: utf-8

"""
This script iterates over a VCF file and computes a Grantham score based on the
eponym matrix.

INFO: Grantham works for all VCF, including indels.
WARNING: Optional inheritance especially designed for Varscan TRIO.
"""

import argparse           # Parse command line arguments
import logging            # Traces and loggings
import logging.handlers   # Logging behaviour
import os                 # Dealing with OS related function
import os.path as op      # Dealing with paths
import pandas as pd       # Parsing large quantification tables
import sys                # System related functions

from pathlib import Path                   # Easily deal with paths
from snakemake.utils import makedirs       # Build output directories
from typing import Any, Dict               # Type hints


logger = logging.getLogger(
    os.path.splitext(os.path.basename(sys.argv[0]))[0]
)


# Building custom class for help formatter
class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    """
    This class is used only to allow line breaks in the documentation,
    without breaking the classic argument formatting.
    """
    pass


# Handling logging options
# No tests for this function
def setup_logging(args: argparse.ArgumentParser) -> None:
    """
    Configure logging behaviour
    """
    root = logging.getLogger("")
    root.setLevel(logging.WARNING)
    logger.setLevel(args.debug and logging.DEBUG or logging.INFO)
    if not args.quiet:
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter(
            "%(levelname)s [%(name)s]: %(message)s"
        ))
        root.addHandler(ch)


def parse_args(args: Any = sys.argv[1:]) -> argparse.ArgumentParser:
    """
    This function parses command lines and returns and ArgumentParser object
    """
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter
    )

    # Parsing positional argument
    main_parser.add_argument(
        "vcf",
        help="Path to a VCF file",
        type=str
    )

    main_parser.add_argument(
        "grantham",
        help="Path to a TSV file containing the Grantham matrix",
        type=str
    )

    main_parser.add_argument(
        "-o", "--output",
        help="Path to output file",
        type=str,
        default="grantham_annotated.vcf",
        metavar="PATH"
    )

    main_parser.add_argument(
        "-f", "--field-separator",
        help="For optional trio annotation, the format field separator",
        type=str,
        default=":",
        metavar="CHAR"
    )

    main_parser.add_argument(
        "-t", "--threshold",
        help="For optional trio annotation, the PVAL threshold",
        type=float,
        default=0.05,
        metavar="FLOAT"
    )

    main_parser.add_argument(
        "-a", "--trio-annotation",
        help="Trio optional inheritance details added into annotation column",
        action="store_true",
        default=False
    )

    # Logging options
    log = main_parser.add_mutually_exclusive_group()
    log.add_argument(
        "-d", "--debug",
        help="Set logging in debug mode",
        default=False,
        action='store_true'
    )

    log.add_argument(
        "-q", "--quiet",
        help="Turn off logging behaviour",
        default=False,
        action='store_true'
    )

    return main_parser.parse_args(args)


def read_grantham(path: Path = "grantham.tsv") -> pd.DataFrame:
    """
    This function returns the given matrix as a DataFrame
    """
    logger.debug("Loading Grantham matrix")
    data = pd.read_csv(
        path,
        sep="\t",
        header=0,
        index_col=0
    )
    logger.debug(data.head())
    return data


def read_vcf(path: Path) -> pd.DataFrame:
    """
    This function reads a VCF file and returns a pandas DataFrame
    """
    logger.debug("Loading VCF file")
    data = pd.read_csv(
        path,
        sep="\t",
        header=None,
        index_col=None,
        comment="#"
    )
    data.columns = [
        "CHROM", "POS", "ID",
        "REF", "ALT", "QUALT",
        "FILTER", "INFO", "FORMAT",
        "Father", "Mother", "Child"
    ]
    logger.debug(data.head())
    return data


def get_format(format: str, field_separator: str = ":") -> Dict[str, int]:
    """
    This function reads the format column and returns a hash of the field
    and its corresponding index.
    """
    return {
        val: idx
        for idx, val in enumerate(format.split(field_separator))
    }


def get_pval_from_format(index: str,
                         patient: str,
                         field_separator: str = ":") -> float:
    """
    This function returns the requested pvalue form a format filed
    """
    try:
        return float(patient.split(field_separator)[index["PVAL"]])
    except ValueError:
        return float(1)


def get_inheritance(format: str,
                    father: str,
                    mother: str,
                    child: str,
                    field_separator: str = ":",
                    threshold: float = 0.05) -> str:
    """
    This function returns a transmition status given a pvalue threshold
    """
    idx = get_format(format, field_separator)
    father_pval = get_pval_from_format(idx, father, field_separator)
    mother_pval = get_pval_from_format(idx, mother, field_separator)
    child_pval = get_pval_from_format(idx, child, field_separator)

    father_mut = father_pval <= threshold
    mother_mut = mother_pval <= threshold
    child_mut = child_pval <= threshold

    if father_mut is True:
        # Father has mutation
        if mother_mut is True:
            # Both father and mother have mutation
            if child_mut is True:
                # Whold trio has mutation
                return "Whole_Trio"
            else:
                # Only parents, but child did not inheritate
                return "Parents_Only"
        else:
            # Only father has mutation
            if child_mut is True:
                # Like father, like child
                return "Inherited_from_Father"
            else:
                # Father did not transmitt
                return "Father_did_not_transmitt_mutation"
    else:
        # Father has no mutation
        if mother_mut is True:
            # Only mother has mutation
            if child_mut is True:
                # Like Mother, like child
                return "Inherited_from_Mother"
            else:
                # Mother did not transmitt
                return "Mother_did_not_transmitt_mutation"
        else:
            # No parent have mutation
            if child_mut is True:
                # Only child has mutation
                return "De_Novo"
            else:
                # No mutation at all
                return "Nobody_has_mutation_under_confidence_interval"


def annotate_trio_inheritance(data: pd.DataFrame,
                              field_separator: str = ":",
                              threshold: float = 0.05) -> pd.DataFrame:
    """
    This function adds an annotation field in order to clarify the inheritance
    process within mutation acquisition in trios
    """
    data["INFO"] = [
        info + ";Mutation_Inheritance=" + get_inheritance(
            format, father, mother, child, field_separator, threshold
        )
        for info, format, father, mother, child in zip(
            data["INFO"], data["FORMAT"], data["Father"],
            data["Mother"], data["Child"]
        )
    ]
    logger.debug(data.head())
    return data


def compute_grantham(ref: str,
                     alt: str,
                     grantham: pd.DataFrame) -> str:
    """
    Compute grantham score
    """
    score = 0
    annot = "No_Change"

    # Easy case shortcut
    if ref == alt:
        return "Grantham_score={};Grantham_annotation={}".format(score, annot)

    # Case of deletions
    elif ref == "." or alt == ".":
        annot = "Radical"
        return "Grantham_score={};Grantham_annotation={}".format(score, annot)

    # Indel support
    for calt in alt:
        for cref in ref:
            try:
                score += grantham[cref][calt]
            except KeyError:
                pass

    # Annotation
    if score > 150:
        annot = "Radical"
    elif score > 100:
        annot = "Moderately_Radical"
    elif score > 50:
        annot = "Moderately_Conservative"
    else:
        annot = "Conservative"

    return "Grantham_score={};Grantham_annotation={}".format(score, annot)


def add_grantham(data: pd.DataFrame, grantham: pd.DataFrame) -> pd.DataFrame:
    """
    This function adds both grantham annotation and score in the given VCF file
    """
    data["INFO"] = [
        info + ";" + compute_grantham(ref, alt, grantham)
        for info, ref, alt in zip(data["INFO"], data["REF"], data["ALT"])
    ]
    return data


def write_vcf_with_comment(vcf: Path,
                           output: Path,
                           data: pd.DataFrame,
                           trio_annotation: bool = False) -> None:
    """
    This function writes the given dataframe whilst conserving the comments
    from original vcf. Some new commands will be added as well.
    """
    logger.debug("Writing results to file")
    mut_in = (
        '##INFO=<ID=Mutation_Inheritance,Number=1,Type=String,'
        'Description="Variant status in trio (Whole_Trio, Parents_Only,'
        ' Inherited_from_Father, Father_did_not_transmitt_mutation,'
        ' Inherited_from_Mother, Mother_did_not_transmitt_mutation,'
        ' De_Novo, Nobody_has_mutation_under_confidence_interval)>"\n'
    )
    grant_score = (
        '##INFO=<ID=Grantham_score,Number=1,Type=Integer,'
        'Description="Grantham score computed from eponym matrix">\n'
    )
    grant_annotation = (
        '##INFO=<ID=Grantham_annotation,Number=1,Type=String,'
        'Description="Grantham annotation obtained from eponym matrix">\n'
    )
    with open(vcf, "r") as vcf_file:
        with open(output, "a") as output_file:
            kept = ""
            for line in vcf_file:

                if line.startswith("#"):
                    if kept != "":
                        output_file.write(kept)
                    kept = line
                else:
                    # The last line contains column header.
                    # We have to add out modifications to the comments
                    if trio_annotation is True:
                        output_file.write(mut_in)
                        output_file.write(grant_score)
                        output_file.write(grant_annotation)
                        output_file.write(kept)
                        data.to_csv(output_file, sep="\t",
                                    header=False, index=False)
                    break


def main(args: argparse.ArgumentParser) -> None:
    """
    This function perfoms the whole annotation process accoding to command line
    options
    """
    data = read_vcf(args.vcf)
    grantham = read_grantham(args.grantham)

    if args.trio_annotation:
        logger.debug("Adding Trio annotations")
        data = annotate_trio_inheritance(
            data, args.field_separator, args.threshold
        )

    logger.debug("Computing Grantham")
    data = add_grantham(data, grantham)

    logger.debug(data.head())
    write_vcf_with_comment(args.vcf, args.output, data, args.trio_annotation)


# Running programm if not imported
if __name__ == '__main__':
    # Parsing command line
    args = parse_args()
    setup_logging(args)

    try:
        logger.debug("Annotating VCF")
        main(args)
    except Exception as e:
        logger.exception("%s", e)
        sys.exit(1)
    sys.exit(0)
