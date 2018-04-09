#!/bin/bash

# This script takes accession numbers as input, then searches for data, and
# copy them on ${PWD}.
# CMD example: bash /data/dayris/Scripts/yawn/get_data.sh M737 MAP213 -fq

# List of error codes
# 1: Unknown command line option
# 2: directory creation error
# 3: Copy error
# 4: Move error
# 5: Write error

# These are the paths where data are searched into
MAPPYACTS="/pandas/Data_Patient/MAPPYACTS"
MOSCATO=("/pandas/Data_Patient/MOSCATOPED" /data_bioinfo/MP/{MOSCATO01,MOSCATO02})
MATCHR="/data_bioinfo/MP/MATCHR"

# instantiation of empty vars
data_pid=()  # The identifier of the data
fastq=false  # Boolean, weather to retrieve fastq only or the whole directory

# Path to DeVa related files
deva_path="/data/dayris/B17_LOVE/DEVA"

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

# Build the list of available data
for option in "$@"; do
  # Testing file project name input
  if [[ "${option}" =~ MAP[0-9]*[_[a-zA-Z]*]? ]]; then
    echo "Looking for this MAPPYACTS data: ${option}"
    data_pid+=($(find "${MAPPYACTS}" -maxdepth 1 -type d -name "${option}" || errorhandling ${LINENO} 2 "Could not find ${option}"))
  elif [[ "${option}" =~ MR[0-9]*[_[a-zA-Z]*]? ]]; then
    echo "Looking for this MATCHR data: ${option}"
    data_pid+=($(find "${MATCHR}" -maxdepth 1 -type d -name "${option}" || errorhandling ${LINENO} 2 "Could not find ${option}"))
  elif [[ "${option}" =~ M[0-9]*[_[a-zA-Z]*]? ]]; then
    echo "Looking for this MOSCATO data: ${option}"
    data_pid+=($(find "${MOSCATO[@]}" -maxdepth 1 -type d -name "${option}" || errorhandling ${LINENO} 2 "Could not find ${option}"))
  # Testing fastq option
  elif [[ "${option}" =~ (-|--)?(F|f)(q|astq|ASTQ)? ]]; then
    echo "Retrieving fastq files only"
    fastq=true
  else
    errorhandling ${LINENO} 1 "Unknown option ${option}"
  fi
done

# Copy data to working directory
for data in "${data_pid[@]}"; do
  echo "Retrieving ${data} , and saving it to ${PWD}"
  if [ ${fastq} == true ]; then
    # Copying files one after each others
    while read file; do
      if ! [[ -f "${PWD}/$(basename ${file})" ]]; then
        cp -v "${file}" "${PWD}" || errorhandling ${LINENO} 3 "Could not copy fastq data: ${data}"
      else
        echo "${file} seems to already be copied. Skipping."
      fi
    done < $(find ${data} -type f -iname "*fastq*")
  else
    # Copying whole directories if requested by user
    if ! [[ -d "${PWD}/$(basename ${data})" ]]; then
      cp -rv "${data}" ${PWD} || errorhandling ${LINENO} 3 "Could not copy data: ${data}"
    else
      echo "${data} seems to already be copied. Skipping."
    fi
  fi
done

# Sort fastq data into a tree understandable by DeVa
# Renaming if needed
if [ ${fastq} == true ]; then
  for file in {.,}*; do
    # Checking weather it's a regular file or not
    if [[ "${file}" =~ M(AP|R)?[0-9]*(-|_)(T|N)[0-9]?(_|\.|-)(ARN|ADN)?(_|\.)?R(1|2)(_|\.)fastq\.gz ]];
    then
      # Project name
      project=$(echo "${file}" | sed -n "s/^\(MR[0-9]*\|MAP[0-9]*\|M[0-9]*\).*$/\1/p")
      if ! [ -d "${project}" ]; then
        # identified as new project
        echo "Creation of the samples directories: ${project} , then copying DEVA"
        mkdir ${project}/{,ADN,ARN} || errorhandling ${LINENO} 2 "Could not create directories"
        cp -r "${deva_path}" "${project}/DEVA" || errorhandling ${LINENO} 3 "Could not copy DeVa"
      fi
      if [[ "${file}" =~ ADN ]]; then
        # identified as DNA
        echo "DNA file identified: ${file}"
        mv "${file}" "${project}/ADN/" || errorhandling ${LINENO} 4 "Could not move ${file}"
        file_path="${project}/ADN/$(basename "${file}")"
        if [[ "${file_path}" =~ ^.*_fastq.gz$ ]]; then
          mv "${file_path}" "${file_path%%_fastq.gz}".fastq.gz || errorhandling ${LINENO} 4 "Could not move ${file_path}"
        fi
      else
        # If regulat file and not DNA, then RNA
        echo "RNA file identified: ${file}"
        mv "${file}" "${project}/ARN/" || errorhandling ${LINENO} 4 "Could not move ${file}"
        file_path="${project}/ADN/$(basename "${file}")"
        if [[ "${file_path}" =~ ^.*_fastq.gz$ ]]; then
          mv "${file_path}" "${file_path%%_fastq.gz}".fastq.gz || errorhandling ${LINENO} 4 "Could not move ${file_path}"
        fi
      fi
    else
      if [ -f ${file} ]; then
        # Not regular file
        echo "Unknown file will be registered for deletion: ${file}"
        echo "rm ${file}" >> cleaning.sh || errorhandling ${LINENO} 5 "Could not write into cleaning.sh"
      elif [ -d ${file} ]; then
        echo "Skipping directory ${file}"
      fi
    fi
  done
fi

# Gunzipping files since DeVa don't like compressed data
while read file; do
  echo "Gunzipping ${file}"
  echo gunzip "${file}" | qsub -N "$(basename "${file}")" -V -d "${PWD}" -j oe -M thibault.dayris@gustaveroussy.fr -m be
done < $(find "${PWD}" -type f -name "*.fastq.gz")

echo "Come back after all unzipping processed are over to run DeVa."
