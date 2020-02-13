#!/bin/bash


function message() {
  local status=${1}
  local message="${2:-1}"

  if [ ${status} = INFO ]; then
    >&2 echo -e "\033[1;36m@INFO:\033[0m ${message}"
  elif [ ${status} = CMD ]; then
    >&2 echo -e "\033[1;32m@CMD:\033[0m ${message}"
  elif [ ${status} = ERROR ]; then
    >&2 echo -e "\033[41m@ERROR:\033[0m ${message}"
  elif [ ${status} = DOC ]; then
    >&2 echo -e "\033[0;33m@DOC:\033[0m ${message}"
  else
    error_handling ${LINENO} 1 "Unknown message type"
  fi
}

function error_handling() {
  # This function will take error messages and exit the program

  # Geathering input parameter (message, third parameter is optionnal)
  echo -ne "\n"
  local parent_lineno="$1"
  local code="$2"
  local message="${3:-1}"

  # Checking the presence or absence of message
  if [[ -n "$message" ]] ; then
    # Case message is present
    message ERROR "Error on or near line ${parent_lineno}:\n ${message}"
    message ERROR "Exiting with status ${code}"
  else
    # Case message is not present
    message ERROR "Error on or near line ${parent_lineno}"
    message ERROR "Exiting with status ${code}"
  fi

  # Exiting with given error code
  exit "${code}"
}

function help_message() {
  message DOC """
  Hi, thanks for using me as your script to prepare automatic snapshots within
  IGV. Please make sure IGV is available in your path.

  Warning: a very long bed file will take ages to be snapshoted by IGV !

  The expected input includes a bed-like file:
  Column 1: CHR
  Column 2: Start
  Column 3: End
  Column 4: Gene name or ID

  Good? So let's begin!

  I understan very fiew things, and here they are:
  -b|--bed4     Path to a bed-like file
  -m|--bam      Path to a bam file
  -g|--genome   Genome ID (hg38, hg19, mm10, ...)
  -h|--help     Print this help then exit
  """
}

function check_file() {
  FILE="${1}"
  if [ ! -f ${FILE} ]; then
    error_handling ${LINENO} 3 "Could not find file: ${FILE}"
  fi
}

function check_dir() {
  DIR="${1}"
  if [ ! -d ${DIR} ]; then
    mkdir --parents ${DIR} || error_handling ${LINENO} 3 "Could not create ${DIR}"
  fi
}

BED4=""
BAM=""
GENOME="hg38"
OUTPUT_DIR=${PWD}

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -b|--bed4) BED4="${2}"; shift 2;;
    -m|--bam) BAM="${2}"; shift 2;;
    -g|--genome) GENOME="${2}"; shift 2;;
    -d|--output-dir) OUTPUT_DIR="${2}"; shift 2;;
    -h|--help) help_message; exit 1;;
    -*|*) error_handling ${LINENO} 1 "Unknown ergument ${1}"; exit 1;;
  esac
done

help_message
message INFO "Checking input files"
check_file ${BED4}
check_file ${BAM}
check_dir ${OUTPUT_DIR}

message INFO "Generating IGV script"
echo "new"
echo "genome ${GENOME}"
echo "load ${BAM}"
echo "snapshotDirectory ${OUTPUT_DIR}"
awk 'BEGIN{FS="\t"} {print "goto "$1":"$2"-"$3"\nsort base\nexpand\nsnapshot "$1":"$2"-"$3"_"$4".png" }' BED4

message INFO "Success!"
