#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script merges multiple ioe files into one tsv with the number of events
per transcript.

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
import traceback

from pandas import DataFrame
from typing import Generator, List


def ioe_generator(*ioe_paths: List[str]) -> Generator[str, str, str]:
    """Yield lines of ioe files one by one, file after file"""
    for path in ioe_paths:
        with open(path, "r") as ioe:
            for line in ioe:
                logging.debug(line[:-1])
                yield line[:-1]


def parse_line(line: str) -> List[str]:
    """Return the parsed line as a list of values"""
    chomp = line.split("\t")
    try:
        gene = chomp[1]
        event = chomp[2]
        transcripts = chomp[3].split(",")
        ev_type = event.split(";")[1].split(":")[0]
        for transcript in transcripts:
            logging.debug([gene, event, transcript])
            yield [gene, event, transcript, ev_type]
    except IndexError:
        return None


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(output=str)
def main(*ioe_paths: "Path to multiple ioe files",
         output: "Path to output file prefix"= "Event_Per_Transcript"):
    """
    Convert multiple ioe files from Suppa to one merged tsv with:
    Col1: Gene ID
    Col2: Transcript ID
    Col3: Event ID
    Col4: Event type
    """
    # result = []
    # for line in ioe_generator(*ioe_paths):
    #     for correspondancy in parse_line(line):
    #         result.append(correspondancy)
    #
    # result = DataFrame(result[1:])

    result = DataFrame(
        [correspondancy
         for line in ioe_generator(*ioe_paths)
         for correspondancy in parse_line(line)][1:]
    )
    result.columns = ["Gene_ID", "Event_ID", "Transcript_ID", "Event_Type"]
    result.drop_duplicates(inplace=True)
    result.set_index(["Transcript_ID", "Event_ID"], inplace=True)
    logging.debug(result.head())
    result.to_csv(output, sep="\t")
