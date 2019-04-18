#!/bin/bash
# author: Thibault Dayris

# This script intends to get the sequences length of fasta files.

# Unlincense terms of use:
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org/>


# Initialization
FASTA_INPUT=""
OUTPUT_FILE=""
SEQUENCES=()

# Command line parsing
while [[ $# -gt 0 ]]; do
  key="${1}"
  case "${key}" in
    -f|--fasta)
      # Path to input fasta file
      FASTA_INPUT="${2}"
      shift; shift
      ;;
    -o|--output)
      # Path to output file
      OUTPUT_FILE="${2}"
      shift; shift
      ;;
    *)
      # Sequences names that are to be measured
      SEQUENCES+=("${1}")
      shift
      ;;
  esac
done

# Checking input parameters
if [ ! -f "${FASTA_INPUT}" ]; then
  echo "ERROR at line ${LINENO}: ${FASTA_INPUT} not found"
  exit 1
fi

if [ -f "${OUTPUT_FILE}" ]; then
  echo "ERROR at line ${LINENO}: ${OUTPUT_FILE} already exists"
  exit 1
fi

if [ ! "$(which samtools)" ]; then
  echo "ERROR at line ${LINENO}: Could not find samtools"
  exit 1
fi

echo "Everything's ok, I'll start to count your sequences' lengths."

# Executing length estimation
for seq in "${SEQUENCES[@]}"; do
  echo -e "${seq}\t$(samtools faidx ${FASTA_INPUT} ${seq} | sed '1d' | wc -c)";
done > "${OUTPUT_FILE}"
