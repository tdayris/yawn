#!/usr/bin/python3
# conding: utf-8

import argparse
import pandas as pd


def salmon_reader(path: "Path to a tsv file") -> pd.DataFrame:
    """
    The salmon specific tsv file reader
    """
    return pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        header=0,
        dtype={
            0: str,
            1: pd.np.float,
            2: pd.np.float,
            3: pd.np.float,
            4: pd.np.float
        },
        na_values=""
    )


def sleuth_reader(path: "Path to a tsv file") -> pd.DataFrame:
    """
    The sleuth specific tsv file reader
    """
    return pd.read_csv(
        path,
        sep="\t",
        index_col=[12, 11, 0],
        header=0,
        dtype={
            0: str,
            1: pd.np.float,
            2: pd.np.float,
            3: pd.np.float,
            4: pd.np.float,
            5: pd.np.float,
            6: pd.np.float,
            7: pd.np.float,
            8: pd.np.float,
            9: pd.np.float,
            10: pd.np.float,
            11: str,
            12: str
        },
        na_values=""
    )


def merges(*salmon):
    merged = None
    for s in salmon:
        data = salmon_reader(s)



def main(salmon, sleuth):
    merges(*salmon)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Perform multiple subsets and merges operations"
                    " between Salmon and Sleuth files. 'Cause these"
                    " merges and subsets are like ALWAYS asked and "
                    "I'm tired to do them manually.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--salmon",
        nargs="+",
        type=str,
        help="Path to multiple salmon files"
    )

    parser.add_argument(
        "--sleuth",
        nargs="+",
        type=str,
        help="Path to multiple sleuth files"
    )

    args = parser.parse_args()

    main(args.salmon, args.sleuth)
