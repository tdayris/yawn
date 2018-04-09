#!/bin/bash

# This script takes a directory in input and sorts files in it.

for file in {.,}*; do
  echo ${file}
  if [[ "${file}" =~ M(AP|R)?[0-9]*(-|_)(T|N)[0-9]?(_|\.|-)(ARN|ADN)?(_|\.)?R(1|2)(_|\.)fastq\.gz ]];
  then
    # echo "Fastq file deteted: ${file}"
    project=$(echo "${file}" | sed -n "s/^\(MR[0-9]*\|MAP[0-9]*\|M[0-9]*\).*$/\1/p")
    if ! [ -d "${project}" ]; then
      echo "Creation of the samples directories: ${project}"
      mkdir ${project}/{,ADN,ARN}
    fi
    if [[ "${file}" =~ ADN ]]; then
      echo "DNA file identified: ${file}"
      mv "${file}" "${project}/ADN/"
    else
      echo "RNA file identified: ${file}"
      mv "${file}" "${project}/ARN/"
    fi
  else
    echo "Unknown file will be deleted: ${file}"
    if [ -f ${file} ]; then
      echo "rm ${file}" >> cleaning.sh
    elif [ -d ${file} ]; then
      echo "Skipping directory ${file}"
    fi
  fi
done
