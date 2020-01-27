#!/bin/bash

work_dir=$(readlink -e "${1}")
command_file="$(date +%s)_Run_Deva_Clinics.sh"

cd "${work_dir}" || exit 1
echo "Now working on ${work_dir}"

while read i; do
  echo "Renaming: ${i}"
  mv "${i}" "${i%%_fastq.gz}".fastq.gz
done < <(find "${PWD}" -type f -name "*_fastq.gz")

while read i; do
  echo "Gunzipping: ${i}"
  echo gunzip "${i}" | qsub -N "$(basename "${i}")" -V -d "${PWD}" -j oe -M thibault.dayris@gustaveroussy.fr -m be
done < <(find "${PWD}" -type f -name "*.fastq.gz")


cd -
echo "Done, back on ${PWD}"

echo "Command will be: bash /data/dayris/DEVA/DEVA_run_clinics.sh ${work_dir}"
echo "bash ${PWD}/DEVA/DEVA_run_clinics.sh ${work_dir}" >> ${command_file}
echo "Commannd line saved at the end of ${command_file}"
