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

import begin
import logging
import os.path
import traceback

import pandas as pd


def parse_se(chomp: "Splitted identifier"):
    # BED starts are 0 based
    first_block_start = int(chomp[2].split("-")[0])
    second_block_start = int(chomp[3].split("-")[0])
    # BED ends are 1 based
    first_block_end = int(chomp[2].split("-")[1]) + 1
    second_block_end = int(chomp[3].split("-")[1]) + 1

    # We need to compute sizes and intervals for the columns 8 to 12.
    # Thanks to the One based system, the length of each section is:
    first_block_size = abs(
        first_block_end - first_block_start
    )

    second_block_size = abs(
        second_block_end - second_block_start
    )

    # We need to find out max and min for columns 2 and 3
    maximum_position = max(
        first_block_end, first_block_start,
        second_block_start, second_block_end
    )
    minimum_position = min(
        first_block_end, first_block_start,
        second_block_start, second_block_end
    )

    return {
        "SecondExonStart": second_block_start,
        "SecondExonStop": second_block_end,
        "FirstExonStart": first_block_start,
        "FirstExonStop": first_block_end,
        "RGB": "100,90,200",
        "Blocks": 2,
        "Sizes": "%i,%i" % (
            first_block_size,
            second_block_size
        ),
        "Starts": "%i,%i" % (
            0,
            second_block_start - minimum_position
        ),
        "MinPos": minimum_position,
        "MaxPos": maximum_position,
        "thickStart": minimum_position,
        "thickEnd": maximum_position,
        "RGB": "255,200,0"
    }


@begin.convert(isa3=begin.utils.tobool, strand=begin.utils.tobool)
def parse_a35(chomp: "splitted identifier",
              isa3: "True if A3, else A5",
              strand: "True if positive, else False"):
    # Find common value
    # Define variable A3/5 site
    # ---==EEEE--- or ---EEE==---
    # With = -> Common junction
    #      E -> Variable donnor/acceptor
    #      - -> Sides

    # BED starts are 0 based
    first_block_start = int(chomp[2].split("-")[0])
    second_block_start = int(chomp[3].split("-")[0])
    # BED ends are 1 based
    first_block_end = int(chomp[2].split("-")[1]) + 1
    second_block_end = int(chomp[3].split("-")[1]) + 1

    # We define the order of the values
    anchors = [first_block_start, second_block_start,
               first_block_end, second_block_end]

    # We need to compute sizes and intervals for the columns 8 to 12.
    # Thanks to the One based system, the length of each section is:
    first_block_size = abs(
        first_block_end - first_block_start
    )

    second_block_size = abs(
        second_block_end - second_block_start
    )

    # We need to find out max and min for columns 2 and 3
    maximum_position = max(
        first_block_end, first_block_start,
        second_block_start, second_block_end
    )
    minimum_position = min(
        first_block_end, first_block_start,
        second_block_start, second_block_end
    )

    if first_block_size > second_block_size:
        # Then we are in the --===EEE-- conformation
        # Forst block is:    --FFFFFF--
        # Second block is:   --SSS-----
        # Either in that sense or the reverse one.
        sizes = "%i,%i" % (
            first_block_size - second_block_size,
            second_block_size
        )
        thick_start = first_block_start
        thick_end = second_block_start
        if strand:
            thick_end, thick_start = thick_start, thick_end
        starts = "%i,%i" % (
            0,
            first_block_size - second_block_size
        )
    else:
        # Then we are in the --===EEE-- conformation
        # Forst block is:    --FFF-----
        # Second block is:   --SSSSSS--
        # Either in that sense or the reverse one.
        sizes = "%i,%i" % (
            second_block_size - first_block_size,
            first_block_size
        )
        thick_start = second_block_start
        thick_end = first_block_start
        if strand:
            thick_end, thick_start = thick_start, thick_end
        starts = "%i,%i" % (
            0,
            second_block_size - first_block_size
        )

    result = {
        "SecondExonStart": second_block_start,
        "SecondExonStop": second_block_end,

        "FirstExonStart": first_block_start,
        "FirstExonStop": first_block_end,
        # Bed ends are one based.
        # A comma-separated list of the block sizes.
        # The number of items in this list should correspond to blockCount
        "Sizes": sizes,
        # A comma-separated list of block starts. All of the blockStart
        # positions should be calculated relative to chromStart.
        # The number of items in this list should correspond to blockCount

        # Yet, this is not possible here: overlapping exons will be
        # conseidered in this BED12 format.
        "Starts": starts,
        "MaxPos": maximum_position,
        "MinPos": minimum_position,
        "thickStart": thick_start,
        "thickEnd": thick_end
    }

    # Colors are arbitrary defined
    if isa3:
        return {"RGB": "255,120,0", **result}
    return {"RGB": "255,160,0", **result}


def parse_ri(chomp: "splitted identifier"):
    # BED starts are 0 based
    intron_block_start = int(chomp[3].split("-")[0])
    # BED ends are 1 based
    intron_block_end = int(chomp[3].split("-")[1]) + 1

    # Alternative acceptor/donnor position
    first_exon_start = int(chomp[2])
    second_exon_start = int(chomp[4])

    # We need to compute sizes and intervals for the columns 8 to 12.
    intron_block_size = abs(
        intron_block_end - intron_block_start
    )
    first_exon_size = abs(
        intron_block_start - first_exon_start
    )
    second_exon_size = abs(
        second_exon_start - intron_block_end
    )

    # We need to find out max and min for columns 2 and 3
    maximum_position = max(
        intron_block_start, intron_block_end,
        first_exon_start, second_exon_start
    )
    minimum_position = min(
        intron_block_start, intron_block_end,
        first_exon_start, second_exon_start
    )

    return {
        "FirstStart": first_exon_start,
        "IntronStart": intron_block_start,
        "IntronStop": intron_block_end,
        "FirstStop": second_exon_start,
        "RGB": "100,90,200",
        "Blocks": 3,
        "Sizes": "%i,%i,%i" % (
            first_exon_size,
            intron_block_size,
            second_exon_size
        ),
        "Starts": "%i,%i,%i" % (
            0,
            intron_block_start - minimum_position,
            intron_block_end - minimum_position
        ),
        "MinPos": minimum_position,
        "MaxPos": maximum_position,
        "thickStart": intron_block_start,
        "thickEnd": intron_block_end
    }


def parse_mx(chomp: "Splitted identifier"):
    # BED starts are 0 based
    first_block_start = int(chomp[2].split("-")[0])
    second_block_start = int(chomp[3].split("-")[0])
    third_block_start = int(chomp[4].split("-")[0])
    fourth_block_start = int(chomp[5].split("-")[0])
    # BED ends are 1 based
    first_block_end = int(chomp[3].split("-")[1]) + 1
    second_block_end = int(chomp[5].split("-")[1]) + 1
    third_block_end = int(chomp[4].split("-")[1]) + 1
    fourth_block_end = int(chomp[5].split("-")[1]) + 1

    # We need to compute sizes and intervals for the columns 8 to 12.
    first_block_size = abs(first_block_end - first_block_start)
    second_block_size = abs(second_block_end - second_block_start)
    third_block_size = abs(third_block_end - third_block_start)
    fourth_block_size = abs(fourth_block_end - fourth_block_start)

    # We need to find out max and min for columns 2 and 3
    maximum_position = max(
        first_block_end, first_block_start, second_block_end,
        second_block_start, third_block_end, third_block_start,
        fourth_block_end, fourth_block_start
    )
    minimum_position = min(
        first_block_end, first_block_start, second_block_end,
        second_block_start, third_block_end, third_block_start,
        fourth_block_end, fourth_block_start
    )

    return {
        "FirstExonStart": first_block_start,
        "FirstExonStop": first_block_end,
        "SecondExonStart": second_block_start,
        "SecondExonStop": second_block_end,
        "ThirdExonStart": third_block_start,
        "ThirdExonStop": third_block_end,
        "FourthExonStart": fourth_block_start,
        "FourthExonStop": fourth_block_end,
        # A comma-separated list of the block sizes.
        # The number of items in this list should correspond to blockCount
        "Sizes": "%i,%i,%i,%i" % (
            first_block_size, second_block_size,
            third_block_size, fourth_block_size
        ),
        # A comma-separated list of block starts. All of the blockStart
        # positions should be calculated relative to chromStart.
        # The number of items in this list should correspond to blockCount
        "Starts": "%i,%i,%i,%i" % (
            0,
            maximum_position - second_block_start,
            maximum_position - third_block_start,
            maximum_position - fourth_block_start,
            # first_block_start, second_block_start,
            # third_block_start, fourth_block_start
        ),
        "MaxPos": maximum_position,
        "MinPos": minimum_position,
        "RGB": "255,220,0",
        "Blocks": 4,
        "thickStart": minimum_position,
        "thickEnd": maximum_position
    }


@begin.convert(isal=begin.utils.tobool)
def parse_afal(chomp: "Splitted identifier",
               isal: "True if event is AL",
               strand: "True if positive stranding"):
    try:
        # BED starts are 0 based
        first_block_start = int(chomp[2].split("-")[0])
        second_block_start = int(chomp[4].split("-")[0])
        # BED ends are 1 based
        first_block_end = int(chomp[2].split("-")[1]) + 1
        second_block_end = int(chomp[4].split("-")[1]) + 1

        # Alternative acceptor/donnor position
        first_ad = int(chomp[3])
        second_ad = int(chomp[5])
    except IndexError:
        # BED starts are 0 based
        first_block_start = int(chomp[3].split("-")[0])
        second_block_start = int(chomp[5].split("-")[0])
        # BED ends are 1 based
        first_block_end = int(chomp[3].split("-")[1]) + 1
        second_block_end = int(chomp[5].split("-")[1]) + 1

        # Alternative acceptor/donnor position
        first_ad = int(chomp[2])
        second_ad = int(chomp[4])

    # We need to compute sizes and intervals for the columns 8 to 12.
    first_block_size = abs(
        first_block_end - first_block_start
    )
    second_block_size = abs(
        second_block_end - second_block_start
    )
    sorted_block_size = list(sorted([first_block_size, second_block_size]))

    # We need to find out max and min for columns 2 and 3
    maximum_position = max(
        first_block_start, first_block_end,
        second_block_start, second_block_end,
        first_ad, second_ad
    )
    minimum_position = min(
        first_block_start, first_block_end,
        second_block_start, second_block_end,
        first_ad, second_ad
    )

    total_size = maximum_position - minimum_position

    result = {
        "FirstExonStart": first_block_start,
        "FirstExonStop": first_block_end,
        "FirstStart": first_ad,
        "SecondExonStart": second_block_start,
        "SecondExonStop": second_block_end,
        "SecondStart": second_ad,
        # A comma-separated list of the block sizes.
        # The number of items in this list should correspond
        # to blockCount
        "Sizes": "%i,%i" % (
            sorted_block_size[0],
            total_size - sorted_block_size[0]
            ),
        # A comma-separated list of block starts. All of the blockStart
        # positions should be calculated relative to chromStart.
        # The number of items in this list should correspond
        # to blockCount
        "Starts": "%i,%i" % (
            0,
            sorted_block_size[0]
        ),
        "MaxPos": maximum_position,
        "MinPos": minimum_position,
        "thickStart": min(first_ad, second_ad),
        "thickEnd": max(second_ad, first_ad)
    }

    if isal:
        result.update(**{"RGB": "80,160,150"})
    else:
        result.update(**{"RGB": "80,160,220"})
    return result


@begin.convert(event=str)
def parse_event(event: "Event identifier without targer id"):
    # Event ID are both : and - separated.
    chomp = event.split(":")

    # These fields are common (except for Blocks, that will be updated)
    result = {
            "EventType": chomp[0],
            "Chromosome": chomp[1],
            "Strand": chomp[-1],
            "Blocks": 2
    }

    if event.startswith("RI"):
        result.update(**parse_ri(chomp))
    elif event.startswith("SE"):
        result.update(**parse_se(chomp))
    elif event.startswith("MX"):
        result.update(**parse_mx(chomp))
    # Many events shall be processed the same way
    elif event.startswith(("A5", "A3")):
        result.update(**parse_a35(chomp, chomp[0] == "A3", chomp[-1] == "+"))
    elif event.startswith(("AF", "AL")):
        result.update(**parse_afal(chomp, chomp[0] == "AL", chomp[-1] == "+"))
    else:
        raise ValueError(event)

    return result


@begin.convert(dpsi=str, nb_cols=int)
def dpsi_event_to_bed(dpsi: "Path to a dpsi file",
                      nb_cols: "Number of columns"=12,
                      drop_dpsi: "Threashold to drop dpsi"=None):
    data = pd.read_csv(dpsi, sep="\t", header=None, index_col=0)
    data.fillna("NaN", inplace=True)
    data.columns = ["DPSI", "P-Val"]
    logging.debug("Original data's shape: %s" % str(data.shape))

    if drop_dpsi is not None:
        data = data[data["DPSI"].astype(float).abs() >= float(drop_dpsi)]
        logging.debug("Removing low dpsi: %s" % str(data.shape))

    data.reset_index(inplace=True)
    data.columns = ["index"] + data.columns.tolist()[1:]
    bed = pd.DataFrame(data["index"].str.split(";", expand=True))
    bed.columns = ["Gene_ID", "Event"]
    bed = pd.merge(
        bed, pd.DataFrame([parse_event(i) for i in bed["Event"]]),
        left_index=True, right_index=True
    )
    bed = pd.merge(bed, data, left_index=True, right_index=True)
    if nb_cols == 3:
        bed = bed[[
            "Chromosome", "MinPos", "MaxPos"
        ]]
    elif nb_cols == 6:
        bed = bed[[
            "Chromosome", "MinPos", "MaxPos",
            "index", "DPSI", "Strand"
        ]]
    elif nb_cols == 9:
        bed = bed[[
            "Chromosome", "MinPos", "MaxPos",
            "index", "DPSI", "Strand",
            "thickStart", "thickEnd", "RGB"
        ]]
    else:
        bed = bed[[
            "Chromosome", "MinPos", "MaxPos",
            "index", "DPSI", "Strand",
            "thickStart", "thickEnd", "RGB",
            "Blocks", "Sizes", "Starts"
        ]]

    logging.debug("Final BED size: %s" % str(bed.shape))
    bed.to_csv("%s.bed" % dpsi, sep="\t", index=False, header=False)


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(dpsi_files=str, dpsi_event=begin.utils.tobool,
               bed3=begin.utils.tobool, bed9=begin.utils.tobool,
               bed12=begin.utils.tobool)
def main(*dpsi_files: "Space separated list of dpsi files",
         dpsi_event: "The dpi file is an event dpsi, to a Transcript one"=True,
         drop_dpsi: "Drop DPSI if value is lower than absolute given value "
                    "(None = No drop out)"=None,
         bed3: "Output a BED3 file instead of BED6"=False,
         bed9: "Output a BED9 file instead of BED6"=False,
         bed12: "Output a BED12 file instead of BED6"=False):
    """
    This file takes a variable number of dpsi files from Suppa and produces
    corresponding bed files in order to be opened within IGV or used within
    downstream analysis.

    WARNING: Bed 12 is not supported
    """
    cols = 6
    if bed3:
        cols = 3
    elif bed12:
        cols = 12

    if dpsi_event:
        for dpsi in dpsi_files:
            if os.path.exists("%s.bed" % dpsi):
                logging.error("%s.bed already exists. Skipping." % dpsi)
                continue
            dpsi_event_to_bed(dpsi, cols, drop_dpsi)
    else:
        raise NotImplementedError("DPSI Event not implemented.")
