#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script is here to build expression tables from Salmon quant.sf to Suppa
"""

import pandas as pd

from argparse import ArgumentParser, RawTextHelpFormatter
from pathlib import Path
from typing import List

__version__ = "1.0.0"
__author__ = ["Thibault Dayris"]

PathType = "Path"
DataFrameType = "DataFrame"


def read_salmon(salmon_path: PathType) -> DataFrameType:
    """Return a filtered DataFrame from a Salmon path"""
    data = pd.read_csv(
        str(salmon_path),
        sep="\t",
        index_col=0,
        header=0,
        dtype={
            0: str, 1: pd.np.int64, 2: pd.np.float64,
            3: pd.np.float64, 4: pd.np.float64
        }
    )

    data = data[["TPM"]]
    data.columns = [salmon_path.parts[-2]]

    return data


def main(*salmon_results: List[str]) -> None:
    """Main function"""
    out_file = Path("merged_Salmon_tables.tsv")
    assert not out_file.exists(), "%s already exists" % str(out_file)

    global_table = None
    for salmon_file in salmon_results:
        print("Working on %s" % salmon_file)
        salmon_path = Path(salmon_file)
        assert salmon_path.exists(), "%s not found" % str(salmon_path)

        salmon_data = read_salmon(salmon_path)

        try:
            global_table = pd.merge(
                global_table,
                salmon_data,
                left_index=True,
                right_index=True)
        except ValueError:
            global_table, salmon_data = salmon_data, global_table

    global_table.to_csv(str(out_file), sep="\t", index_label=False)


if __name__ == '__main__':
    parser = ArgumentParser(
        description="Build a unique, merged quantification table from Salmon"
                    " quant.sf files. Especially designed for Suppa.",
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "salmon_result",
        nargs="+",
        type=str,
        help="Path to a quant.sf file from Salmon 0.8.0 or higher"
    )

    args = parser.parse_args()

    main(*args.salmon_result)
