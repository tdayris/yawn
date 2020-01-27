#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script builds a CLS file from command line

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
