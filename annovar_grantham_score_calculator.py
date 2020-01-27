#!/usr/bin/python3.7
# conding: utf-8

"""
This script iterates over a VCF file and computes a Grantham score based on the
eponym matrix.

INFO: Grantham works for all VCF, including indels.
WARNING: Optional inheritance especially designed for Varscan TRIO.
"""

import argparse
import logging            # Traces and loggings
import logging.handlers   # Logging behaviour
import os                 # Dealing with OS related function
import os.path as op      # Dealing with paths
import pandas as pd       # Parsing large quantification tables
import sys                # System related functions

from pathlib import Path              # Easily deal with paths
from snakemake.utils import makedirs  # Build output directories
from typing import List, Any               # Type hints


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
def setup_logging() -> None:
    """
    Configure logging behaviour
    """
    root = logging.getLogger("")
    root.setLevel(logging.WARNING)
    logger.setLevel(logging.DEBUG)
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
        "--field-separator",
        help="For optional trio annotation, the format field separator",
        type=str,
        default=":"
    )

    main_parser.add_argument(
        "--threshold",
        help="For optional trio annotation, the PVAL threshold",
        type=float,
        default=0.05
    )

    main_parser.add_argument(
        "--trio-annotation",
        help="Trio optional inheritance details added into annotation column",
        action="store_true"
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
    return pd.read_csv(
        path,
        sep="\t",
        header=0,
        index_col=0
    )


def read_vcf(path: Path) -> pd.DataFrame:
    """
    This function reads a VCF file and returns a pandas DataFrame
    """
    data = pd.read_csv(
        path,
        sep="\t",
        header=None,
        index_col=None,
        commant="#"
    )
    data.columns = [
        "CHROM", "POS", "ID",
        "REF", "ALT", "QUALT",
        "FILTER", "INFO", "FORMAT",
        "Father", "Mother", "Child"
    ]
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
    return float(patient.split(field_separator)[index["PVAL"]])


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

    if father_mut:
        # Father has mutation
        if mother_mut:
            # Both father and mother have mutation
            if child_mut:
                # Whold trio has mutation
                return "Whole_Trio"
            else:
                # Only parents, but child did not inheritate
                return "Parents_Only"
        else:
            # Only father has mutation
            if child_mut:
                # Like father, like child
                return "Inherited_from_Father"
            else:
                # Father did not transmitt
                return "Father_did_not_transmitt_mutation"
    else:
        # Father has no mutation
        if mother_mut:
            # Only mother has mutation
            if child_mut:
                # Like Mother, like child
                return "Inherited_from_Mother"
            else:
                # Mother did not transmitt
                return "Mother_did_not_transmitt_mutation"
        else:
            # No parent have mutation
            if child_mut:
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
    return data


def compute_grantham(ref: str,
                     alt: str,
                     grantham: pd.DataFrame) -> List[Union[int, str]]:
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
            score += grantham[cref][catl]

    # Annotation
    if score > 150:
        annot = "Radical"
    elif score > 100:
        annot = "Moderately_Radical"
    elif score > 50:
        annot = "Moderately_Conservative"
    annot = "Conservative"

    return "Grantham_score={};Grantham_annotation={}".format(score, annot)


def add_grantham(data: pd.DataFrame, grantham: pd.DataFrame) -> pd.DataFrame:
    """
    This function adds both grantham annotation and score in the given VCF file
    """
    data["INFO"] = [
        info + compute_grantham(ref, alt, grantham)
        for info, ref, alt in zip(data["INFO"], data["REF"], data["ALT"])
    ]
    return data


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
    data.to_csv(args.output, sep="\t")


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
