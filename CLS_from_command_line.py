#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script builds a CLS file from command line
"""

import begin
import logging
import os.path
import pandas as pd
import traceback

from typing import List


@begin.start
@begin.convert(output_file=str)
@begin.logging
@begin.tracebacks
def main(classes: "Comma separated list of class names",
         nb_cls: "Comma separated list of number of samples per class",
         output_file: "Path to output file"="class.cls"):
    """
    This script builds a CLS file from command line
    """
    if os.path.exists(output_file):
        raise FileExistsError("%s already exists" % output_file)

    # Split comma separated input
    classes = classes.split(",")
    nb_cls = nb_cls.split(",")

    # Classes mus have unique names
    if any(classes.count(x) != 1 for x in classes):
        raise ValueError("Duplicate class name!")

    # Save file
    with open(output_file, "w") as outfile:
        outfile.write("%i %i 1\n" % (sum(map(int, nb_cls)), len(classes)))
        outfile.write("# %s\n" % " ".join(classes))
        outfile.write(
            "%s\n" % " ".join(
                [" ".join([cls]*int(nb)) for cls, nb in zip(classes, nb_cls)]
            )
        )
