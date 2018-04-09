#!/bin/bash

# This script takes accession numbers as input, then searches for data, and
# copy them on ${PWD}.

MAPPYACTS="/pandas/Data_Patient/MAPPYACTS"
MOSCATO=("/pandas/Data_Patient/MOSCATOPED" /data_bioinfo/MP/{MOSCATO01,MOSCATO02})
MATCHR="/data_bioinfo/MP/MATCHR"
data_pid=()  # The identifier of the data
fastq=false  # Boolean, weather to retrieve fastq only or the whole directory

errorhandling() {
  # This function will take error messages and exit the program

  # Geathering input parameter (message, third parameter is optionnal)
  echo -ne "\n"
  local parent_lineno="$1"
  local code="$2"
  local message="${3:-1}"

  # Checking the presence or absence of message
  if [[ -n "$message" ]] ; then
    # Case message is present
    echo -e "Error on or near line ${parent_lineno}:\n ${message}"
    echo "Exiting with status ${code}"
  else
    # Case message is not present
    echo "Error on or near line ${parent_lineno}"
    echo "Exiting with status ${code}"
  fi

  # Exiting with given error code
  exit "${code}"
}

for option in "$@"; do
  if [[ "${option}" =~ MAP[0-9]*[_[a-zA-Z]*]? ]]; then
    echo "Looking for this MAPPYACTS data: ${option}"
    data_pid+=($(find "${MAPPYACTS}" -maxdepth 1 -type d -name "${option}" || errorhandling ${LINENO} 2 "Could not find ${option}"))
  elif [[ "${option}" =~ MR[0-9]*[_[a-zA-Z]*]? ]]; then
    echo "Looking for this MATCHR data: ${option}"
    data_pid+=($(find "${MATCHR}" -maxdepth 1 -type d -name "${option}" || errorhandling ${LINENO} 2 "Could not find ${option}"))
  elif [[ "${option}" =~ M[0-9]*[_[a-zA-Z]*]? ]]; then
    echo "Looking for this MOSCATO data: ${option}"
    data_pid+=($(find "${MOSCATO[@]}" -maxdepth 1 -type d -name "${option}" || errorhandling ${LINENO} 2 "Could not find ${option}"))
  elif [[ "${option}" =~ (-|--)?(F|f)(q|astq|ASTQ) ]]; then
    echo "Retrieving fastq files only"
    fastq=true
  else
    errorhandling ${LINENO} 1 "Unknown option ${option}"
  fi
done

for data in "${data_pid[@]}"; do
  echo "Retrieving ${data} , and saving it to ${PWD}"
  if [ ${fastq} == true ]; then
    while read file; do
      if ! [[ -f "${PWD}/$(basename ${file})" ]]; then
        cp -v "${file}" "${PWD}" || errorhandling ${LINENO} 2 "Could not copy fastq data: ${data}"
      else
        echo "${file} seems to already be copied. Skipping."
      fi
    done < $(find ${data} -type f -iname "*fastq*")
  else
    if ! [[ -d "${PWD}/$(basename ${data})" ]]; then
      cp -rv "${data}" ${PWD} || errorhandling ${LINENO} 2 "Could not copy data: ${data}"
    else
      echo "${data} seems to already be copied. Skipping."
    fi
  fi
done
