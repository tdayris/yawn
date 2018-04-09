#!/bin/bash

# This script takes accession numbers as input, then searches for data, and
# copy them on ${PWD}.

data_paths=(/pandas/Data_Patient/{MAPPYACTS,MOSCATOPED} /data_bioinfo/MP/{MOSCATO01,MOSCATO02,MATCHR})
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
  case "${option}" in
    MAP[0-9]*|MR[0-9*]|M[0-9]*)
      echo -e "Will be looking for ${option}"
      data_pid+=("${option}")
      ;;
    --fastq|fastq|--fq)
      echo "Will only retrieve fastq data"
      fastq=true
      ;;
    -*|--|*)
      echo "Unknown option"
      errorhandling ${LINENO} 1 "Unknown option"
  esac
done

for data in "${data_pid[@]}"; do
  if ${fastq} == true; then
    while read path; do
      find ${path} -type f -iname "*${data}*fastq*" -exec bash -c echo {} ${PWD} \;
    done < $(find "${data_paths[@]}" -maxdepth 1 -type d -name "${data}")
  else
    find "${data_paths[@]}" -maxdepth 1 -type d -name "${data}" -exec bash -c echo {} ${PWD} \;
  fi
done
