#!/bin/bash

# This script takes accession numbers as input, then searches for data, and
# copy them on ${PWD}.
# CMD example: bash get_data.sh M737 MAP213 -fq -unzip

# List of error codes
# 1: Unknown command line option
# 2: directory creation error
# 3: Copy error
# 4: Move error
# 5: Write error
# 6: DEVA_create_rawdata error

# These are the paths where data are searched into
MAPPYACTS="/pandas/Data_Patient/MAPPYACTS"
MOSCATO=("/pandas/Data_Patient/MOSCATOPED" /data_bioinfo/MP/{MOSCATO01,MOSCATO02})
MATCHR="/data_bioinfo/MP/MATCHR"

# instantiation of empty vars
data_pid=()  # The identifier of the data
projects_list=()  # The name of the projects
fastq=f  # Boolean, weather to retrieve fastq only or the whole directory
raw=f  # Perform a gunzip routine
create_raw_data=f  # Perform a raw data creation for DeVa

# Path to DeVa related files
deva_path="/data/dayris/B17_LOVE/DEVA"
deva_run_script="run_all_deva.sh"

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
  if [[ "${option}" =~ (MAP|NK_)[0-9]*[_?[a-zA-Z]*]? ]]; then
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
    fastq=t
  elif [[ "${option}" =~ (-|--)?g?unzip ]]; then
    raw=t
  elif [[ "${option}" =~ (-|--)?c(reate)? ]]; then
    create_raw_data=t
  else
    errorhandling ${LINENO} 1 "Unknown option ${option}"
  fi
done

# Copy data to working directory
for data in "${data_pid[@]}"; do
  echo "Retrieving ${data} , and saving it to ${PWD}"
  if [ ${fastq} == t ]; then
    # Copying files one after each others
    while read file; do
      if ! [[ -f "${PWD}/$(basename ${file})" ]]; then
        cp -v "${file}" "${PWD}" || errorhandling ${LINENO} 3 "Could not copy fastq data: ${data}"
      else
        echo "${file} seems to already be copied. Skipping."
      fi
    done < <(find ${data} -type f -iname "*fastq*")
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
if [ ${fastq} == t ]; then
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
        projects_list+=("${project}")
      fi
      if [[ "${file}" =~ (ADN|(-|_)N) ]]; then
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
if [ "${raw}" == t ]; then
  while read file; do
    echo "Gunzipping ${file}"
    echo gunzip "${file}" | qsub -N "$(basename "${file}")" -V -d "${PWD}" -j oe -M thibault.dayris@gustaveroussy.fr -m be
  done < <(find "${PWD}" -type f -name "*.fastq.gz")

  echo "The script will now wait until there are no gzippe files anymore."
else
  echo "Files are prepared."



  while [[ $(find "${PWD}" -type f -name "*.fastq.gz") -gt 0 ]]; do
    echo "Waiting 5 minutes bicause there are compressed fastq files"
    sleep 300
  done
  echo "All files seems gunzipped."
fi

# Trying to create rawdata for deva
if [[ "${create_raw_data}" == t ]]; then
  for project in "${projects_list[@]}"; do
    echo "${project} is being prepared for DeVa."
    arn=(${project}/ARN/*)
    adn=(${project}/ADN/*)
    if [ ${#arn[@]} -gt 0 ]; then
      echo "DNA present in ${project}"
      (
        cd "${project}/DEVA" || exit
        echo "${PWD}/../ARN" | bash DEVA_create_rawdata.sh RNA || errorhandling ${LINENO} 6 "Could not create RNA rawdata for ${project}"
        echo "${PWD}/rawData_RNA.tsv created"
        cd -
      )
      echo "cd ${project}/DEVA" >> "${deva_run_script}"
      echo "bash run_deva.sh rawData_RNA.tsv nocorrect RNA Fusion" >> "${deva_run_script}"
      echo "bash run_deva.sh rawData_RNA.tsv nocorrect RNA VariantCalling" >> "${deva_run_script}"
      echo "cd -" >> "${deva_run_script}"
    fi
    if [ ${#adn[@]} -gt 0 ]; then
      echo "DNA present in ${project}"
      (
        cd "${project}/DEVA" || exit
        echo "${PWD}/../ARN" | bash DEVA_create_rawdata.sh WES || errorhandling ${LINENO} 6 "Could not create WES rawdata for ${project}"
        echo "${PWD}/rawData_WES.tsv created"
        cd -
      )
      echo "cd ${project}/DEVA" >> "${deva_run_script}"
      echo "bash run_deva.sh rawData_RNA.tsv nocorrect WES" >> "${deva_run_script}"
      echo "cd -" >> "${deva_run_script}"
    fi
  done
fi
