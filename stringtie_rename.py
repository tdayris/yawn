#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script is here to handle StringTie unpleasant habit
to rename known transcripts
"""

import argparse

from collections import OrderedDict
from pathlib import Path
from typing import List, Generator

__version__ = "1.0.0"
__author__ = ["Thibault Dayris"]

Pathtype = "Path"


def gtf_generator(file_path: Pathtype, merged=False) \
        -> Generator[str, str, str]:
    """Return a generator of modified gtf lines"""
    with file_path.open() as gtf:
        for line in gtf:
            if line.startswith("#"):
                yield line[:-1]
                continue

            if "ref_gene_id" not in line:
                yield line[:-1]
                continue

            chomp = line[:-1].split("\t")

            attr = OrderedDict()
            for i in chomp[-1].split(";"):
                try:
                    key, value = i.split('"')[:2]
                    attr[key.strip()] = value.strip()
                except ValueError:
                    continue

            attr["gene_id"] = attr["ref_gene_id"]
            del attr["ref_gene_id"]
            if not merged:
                attr["gene_name"] = attr["ref_gene_name"]
                del attr["ref_gene_name"]
                attr["transcript_id"] = attr["reference_id"]
                del attr["reference_id"]

            str_attr = " ".join([
                "%s \"%s\";" % (key, value)
                for key, value in attr.items()
            ])

            yield "\t".join(chomp[0:8] + [str_attr])


def main(*gtf: List[str], merged=False) -> None:
    """Main function"""
    for gtf_file in gtf:
        print("Working on %s" % gtf_file)
        gtf_path = Path(gtf_file)
        assert gtf_path.exists(), "%s not found" % str(gtf_path)

        out_file = Path(gtf_path.parent, "%s_renamed.gtf" % gtf_path.stem)
        assert not out_file.exists(), "%s already exists" % str(out_file)

        out_file.write_text("\n".join(list(gtf_generator(gtf_path, merged))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Remove the known transcripts renaming done by StringTie",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "gtf",
        nargs="+",
        type=str,
        help="Path to one or more GTF files (NO GFF)"
    )

    parser.add_argument(
        "--merged",
        help="The gtf comes from stringtie merge, not raw stringtie",
        action="store_true"
    )

    args = parser.parse_args()

    main(*args.gtf, merged=args.merged)
