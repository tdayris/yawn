#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

import begins
import logging
import traceback

from Bio import SeqIO


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(fasta_path=str, output_path=str)
def main(*fasta_paths: "Space separated list of paths to fasta file",
         frame: "The frame shift to apply to each sequence" = 0,
         reverse: "Weather to revers complement or not" = False
         output_path: "Path to output file" = "AA.fasta"):
    """
    This script translates DNA sequence from fasta file into all 6 ORF in AA
    """
    for fasta_file in fasta_paths:
        for sequence in SeqIO.parse(fasta_file, "fasta"):
            logging.debug(sequence)
            print(sequence.id)
            if reverse is True:
                print(sequence.seq.reverse_complement()[frame:].translate())
            else:
                print(sequence.seq[frame:].translate())
