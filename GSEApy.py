#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Perform a GSEA analysis

Unlincense terms of use:
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"""

import argparse
import os.path
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

from typing import List


def main(cls_file: str,
         gmt_like_file: str,
         gene_set: str='GO_Biological_Process_2017b',
         permutation_type: str='phenotype',
         method: str="ratio_of_classes",
         output_dir: str="gsea_report",
         format: str="pdf",
         permutation_num: int=1000,
         threads: int=1) -> None:
    """
    Perform GSEA processing
    """
    phenoA, phenoB, class_vector = gp.parser.gsea_cls_parser(cls_file)

    gene_exp = pd.read_table(gmt_like_file, header=0, index_col=0, comment="#")

    gs_res = gp.gsea(
        data=gene_exp,
        cls=class_vector,
        permutation_type=permutation_type,
        permutation_num=permutation_num,
        outdir=output_dir,
        gene_sets=gene_set,
        method=method,
        processes=threads,
        format=format
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "cls_file",
        help="Path to a cls file",
        type=str
    )

    parser.add_argument(
        "gmt_like_file",
        help="Path to a gmt-like file with second line discarded",
        type=str
    )

    parser.add_argument(
        "-g", "--gene_set",
        help="Gene set name",
        type=str,
        choices=gp.get_library_name(),
        default='KEGG_2017'
    )

    parser.add_argument(
        "-p", "--permutation_type",
        help="Type of permutation used within GSEA",
        type=str,
        choices=['gene_set', 'phenotype'],
        default='phenotype'
    )

    parser.add_argument(
        "-o", "--output_dir",
        help="Path to output_dir",
        type=str,
        default='gsea_report'
    )

    parser.add_argument(
        "-f", "--format",
        help="Output format for all figures",
        type=str,
        choices=['pdf', 'png', 'jpeg', 'ps', 'eps', 'svg'],
        default="pdf"
    )

    parser.add_argument(
        "-n", "--permutation_num",
        help="Number of permutation used within GSEA",
        type=int,
        default=1000
    )

    parser.add_argument(
        "-m", "--method",
        help="Methods to calculate correlations of ranking metrics",
        type=str,
        choices=['signal_to_noise', 't_test', 'ratio_of_classes',
                 'diff_of_classes', 'log2_ratio_of_classes'],
        default="ratio_of_classes"
    )

    parser.add_argument(
        "-t", "--threads",
        help="Maximum number of threads used",
        type=int,
        default=1
    )

    args = parser.parse_args()
    main(**vars(args))
