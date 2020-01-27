#!/usr/bin/python3
# Author: Dayris Thibault

"""
This script extracts a tsv from a GTF. This tsv contains gene_id,
transcript_id, gene_name.
"""

import argparse


def process_gtf(gtf_path: str, positions: bool = False, sep: str = "\t"):
    """
    Build a Transcript to Gene table.

    TSV text file with 3 columns:
    - Gene ID
    - Transcript ID
    - Gene Name

    If positions is True, then two columns are added:
    - Chromosome
    - Start
    - Stop

    No column name given.
    """
    with open(gtf_path, "r") as gtf_file:
        for line in gtf_file:
            if line.startswith("#"):
                # The line is a comment if it starts wirh #
                continue

            # If the line is not a transcript, then it will no contain
            # the required information
            chomp = line[:-1].split("\t")
            if chomp[2] != "transcript":
                continue

            start = chomp[3]
            stop = chomp[4]
            chrom = chomp[0]
            strand = chomp[6]

            # We are interested in the last column
            chomp = {
                attr.split('"')[0].strip(): attr.split('"')[1].strip()
                for attr in chomp[8].split(";")
                if attr != '' and '"' in attr
            }

            # Some genes have an ID but no name ...
            try:
                result = [
                    chomp["gene_id"],
                    chomp["transcript_id"],
                    chomp["gene_name"]
                ]
            except KeyError:
                result = [
                    chomp["gene_id"],
                    chomp["transcript_id"],
                    chomp["gene_id"]
                ]

            if positions:
                print(sep.join(result + [chrom, start, stop, strand]))
            else:
                print(sep.join(result))


def main():
    main_parser = argparse.ArgumentParser(
        description="""Build a Transcript to Gene table
TSV text file with three columns:
    1- Gene ID
    2- Transcript ID
    3- Gene Name

    If positions is True, then three columns are added:
    4- Chromosome
    5- Start
    6- Stop
        """,
        epilog="Designed for GTF files only",
        formatter_class=argparse.RawTextHelpFormatter
    )

    main_parser.add_argument(
        "gtf_path",
        help="Path to a GTF file",
        metavar="PATH",
        type=str
    )

    main_parser.add_argument(
        "-c", "--colnames",
        help="Add column names (default: False)",
        action="store_true"
    )

    main_parser.add_argument(
        "-p", "--positions",
        help="Add start and stop for each transcript",
        action="store_true"
    )

    main_parser.add_argument(
        "--csv",
        help="Print output as CSV file instead of TSV",
        action="store_true"
    )

    args = main_parser.parse_args()

    sep = "\t"

    if args.colnames:
        cols = (
            ["Gene_ID", "Transcript_ID", "Gene_Name",
             "Chromosome", "Start", "Stop", "Strand"]
            if args.positions else
            ["Gene_ID", "Transcript_ID", "Gene_Name"]
        )
        sep = ("," if args.csv else "\t")
        print(sep.join(cols))

    process_gtf(args.gtf_path, args.positions, sep)


if __name__ == '__main__':
    main()
