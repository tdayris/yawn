#!/bin/bash

# This function will take error messages and exit the program
function error_handling() {
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

# This function only changes echo headers
# for user's sake.
function message() {
  # Define local variables
  local status=${1}         # Either INFO, CMD, ERROR or DOC
  local message="${2:-1}"   # Your message

  # Classic switch based on status
  if [ ${status} = INFO ]; then
    echo -e "\033[1;36m@INFO:\033[0m ${message}"
  elif [ ${status} = CMD ]; then
    echo -e "\033[1;32m@CMD:\033[0m ${message}"
  elif [ ${status} = ERROR ]; then
    >&2 echo -e "\033[41m@ERROR:\033[0m ${message}"
  elif [ ${status} = DOC ]; then
    echo -e "\033[0;33m@DOC:\033[0m ${message}"
  else
    error_handling ${LINENO} 1 "Unknown message type"
  fi
}

function read_check(){
  local name="${1}"
  local seq="${2}"
  local recall="${3}"
  local qual="${4}"

  # 1. Read name exists
  if [ "${name}" == "" ]; then
    error_handling ${LINENO} 1 "No name was given for the read..."
  fi

  # 2. Sequence exists
  if [ "${seq}" == "" ]; then
    error_handling ${LINENO} 1 "No sequence was given for the read ${name}"
  fi

  # 3. Recall exists
  if [ "${recall}" == "" ]; then
    error_handling ${LINENO} 1 "No name recall was given for the read ${name}"
  fi

  # 4. Quality exists
  if [ "${qual}" == "" ]; then
    error_handling ${LINENO} 1 "No quality score was given for the read ${name}"
  fi

  # 5. Name starts with @
  if [[ "${name}" != @* ]]; then
    error_handling ${LINENO} 2 "Name does not start with '@' in read ${name}"
  fi

  # 6. Recall starts with +
  if [[ "${recall}" != +* ]]; then
    error_handling ${LINENO} 2 "Name recall does not start with '+' in read ${name}: ${recall}"
  fi

  # 7. Seq and qual have the same length
  if [[ ${#seq} -ne ${#qual} ]]; then
    error_handling ${LINENO} 3 "Quality and Sequence do not have mathing length in read: ${name}:\n${seq}\n${qual}"
  fi

  # 8. Qual must be in: !"#$%&'\''()*+,-.\/0123456789:;<=>?@ABCDEFGHIJKL
  if [[ "${qual}" =~ ^![\!\"\#\$%\&\'\(\)\*\+,-\./0123456789:\;\<\=\>\?@ABCDEFGHIJKL]*$ ]]; then
    error_handling ${LINENO} 3 "Quality score contains unexpected value in ${name}: ${qual}"
  fi

  # 9. Seq must be in: ATGCNatgcn
  if [[ "${qual}" =~ ^![ATGCNatgcn]*$ ]]; then
    error_handling ${LINENO} 3 "Read sequence score contains unexpected value in ${name}: '${seq}'"
  fi
}

function help_message() {
  message DOC "Hi, thanks for using me as your script for controlling"
  message DOC "your fastq files format. I'm very proud to be your script"
  message DOC "today, and I hope you'll enjoy working with me."
  echo ""
  message DOC "Every time you'll see a line starting with '@', it will be because I speak."
  message DOC "In fact, I always start my speech with '\033[0;33m@DOC:\033[0m' when i't about my functions,"
  message DOC "'\033[1;36m@INFO:\033[0m' when it's about my intentions, "
  message DOC "with '\033[1;32m@CMD:\033[0m' I give you all the command lines I use,"
  message DOC "and with '\033[41m@ERROR:\033[0m', I tell you when things go wrong."
  echo ""
  message DOC "I understand very fiew things, and here they are:"
  # message DOC "-p | --process     The number of process per batch."
  message DOC "-h | --help        Print this help message, then exit."
  echo ""
  message DOC "A typical command line would be:"
  message DOC "bash fastq_format_checker.sh my_fq_file.fq"
  message DOC "bash fastq_format_checker.sh <(gunzip -c my_fq_file.fq.gz)"
  echo ""
  message DOC "If you dig in my code, you'll see that I have everything needed"
  message DOC "to be multithreaded. However, my read check function it so damn"
  message DOC "fast, that spawning a thread to make the job slows down my"
  message DOC "final runtime. So ... Single threaded I am."
  exit 0
}

# Global variables
# Maximum number of process per batch
PROCESS=1
# Positional arguments
FQFILE=""

while [ "$#" -gt 0 ]; do
  # Parse command line
  case "${1}" in
    # -p|--process) PROCESS="${2}"; shift 2;;
    -h|--help) help_message; exit 1;;
    -*) error_handling ${LINENO} 1 "Unknown arguments ${1}"; exit 1;;
    *) FQFILE="${1}"; shift;;
  esac
done

message INFO "Number of process: ${PROCESS}"
message INFO "Path to fastq file: ${FQFILE}"

# Uncomment the following to activate multi threads
# However, comment the next section!
# i=0
# (
# while IFS=$'\t' read NAME SEQ RECALL QUAL; do
#   ((i=i%PROCESS)); ((i++==0)) && wait
#   message INFO "Checking: ${NAME}"
#   read_check "${NAME}" "${SEQ}" "${RECALL}" "${QUAL}" &
# done < <(cat ${FQFILE} | paste - - - -); wait
# )

# Comment the following to deactivate monothreaded activity
# However, uncomment the previous section!
while IFS=$'\t' read NAME SEQ RECALL QUAL; do
  message INFO "Checking: ${NAME}"
  read_check "${NAME}" "${SEQ}" "${RECALL}" "${QUAL}"
done < <(cat ${FQFILE} | paste - - - -)
