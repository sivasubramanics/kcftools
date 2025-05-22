#!/bin/bash
# This script runs the kcftools program with the specified parameters.
# the inputs are the tab delimited file (data) with below columns
#   1. sample name
#   2. fq1
#   3. fq2
# reference genome fasta file
# kmer size
# window length
# output directory

set -u
set -e
set -o pipefail

# start clock
start=$(date +%s)

# print log
print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] LOG: $1"
}

# run command
run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

# convert seconds to hh:mm:ss
seconds_to_hhmmss() {
  local seconds=$1
  local hours=$(printf "%02d" $((seconds / 3600)))
  local minutes=$(printf "%02d" $(((seconds / 60) % 60)))
  local seconds=$(printf "%02d" $((seconds % 60)))
  echo "$hours:$minutes:$seconds"
}

# usage function
usage(){
    echo "Usage: $0 -d <data_file> -r <ref_genome> -k <kmer_size> -w <window_length> -o <output_dir> -t <n_threads> -m <mem>"
    echo "Options:"
    echo "  -d  tab delimited file with sample_name, fq1, fq2"
    echo "  -r  reference genome fasta file"
    echo "  -n  reference genome name"
    echo "  -k  kmer size for analysis [31]"
    echo "  -w  window length for analysis [5000]"
    echo "  -o  directory to save output files [output]"
    echo "  -t  number of threads to use [4]"
    echo "  -m  memory to use for each thread (in GB) [12]"
}

# set default values
kmer_size=31
window_length=5000
output_dir="k${kmer_size}"
n_threads=4
mem=12

# parse command line arguments
while getopts ":d:r:n:k:w:o:t:m:h" opt; do
  case ${opt} in
    d) data_file="${OPTARG}"
    ;;
    r) ref_genome="${OPTARG}"
    ;;
    n) ref_name="${OPTARG}"
    ;;
    k) kmer_size="${OPTARG}"
    ;;
    w) window_length="${OPTARG}"
    ;;
    o) output_dir="${OPTARG}"
    ;;
    t) n_threads="${OPTARG}"
    ;;
    m) mem="${OPTARG}"
    ;;
    h) usage
       exit 0
    ;;
    \?) echo "Invalid option -${OPTARG}" >&2
        usage
        exit 1
    ;;
  esac
done

# check if required arguments are provided
if [ -z "${data_file:-}" ] || [ -z "${ref_genome:-}" ] || [ -z "${ref_name:-}" ]; then
    echo "Error: data_file and ref_genome are required arguments."
    usage
    exit 1
fi

# creat kmc database for each sample from the data file
print_log "Creating kmc database for each sample from the data file"
mkdir -p ${output_dir}
while IFS=$'\t' read -r sample fq1 fq2; do
    print_log "FASTQ2KMC sample: $sample"
    if [ -f "${output_dir}/${sample}.k${kmer_size}.kmc.done" ] && [ -f "${output_dir}/${sample}.k${kmer_size}.kmc_pre" ] && [ -f "${output_dir}/${sample}.k${kmer_size}.kmc_suf" ]; then
        print_log "Sample $sample already processed, skipping..."
        continue
    fi
    # create fastq list file, with each line containing fq file name
    echo "${fq1}" > ${output_dir}/$sample.fq.list
    echo "${fq2}" >> ${output_dir}/$sample.fq.list
    mkdir -p ${output_dir}/${sample}_kmc_tmp
    run_cmd "{ time kmc \
        -k${kmer_size} \
        -m${mem} \
        -t${n_threads} \
        -c0 \
        -p9 \
        -fq @${output_dir}/${sample}.fq.list \
        ${output_dir}/${sample}.k${kmer_size} \
        ${output_dir}/${sample}_kmc_tmp ; } &> ${output_dir}/$sample.k${kmer_size}.kmc.log 2>&1"
    touch ${output_dir}/${sample}.k${kmer_size}.kmc.done
done < "${data_file}"


# screen for kcftools variations using getVariations plugin
print_log "Screening for kcftools variations using getVariations plugin"
while IFS=$'\t' read -r sample fq1 fq2; do
    print_log "Screening for kcftools variations using getVariations plugin for sample: $sample"
    if [ -f "${output_dir}/${ref_name}.${sample}.k${kmer_size}.w${window_length}.kcf.done" ] && [ -f "${output_dir}/${ref_name}.${sample}.k${kmer_size}.w${window_length}.kcf" ]; then
        print_log "Sample $sample already processed, skipping..."
        continue
    fi
    run_cmd "{ time kcftools getVariations \
        -k ${output_dir}/${sample}.k${kmer_size} \
        -r ${ref_genome} \
        -s ${sample} \
        -o ${output_dir}/${ref_name}.${sample}.k${kmer_size}.w${window_length}.kcf \
        -w ${window_length} \
        -f window \
        -m \
        -t ${n_threads} ; } &> ${output_dir}/${ref_name}.${sample}.k${kmer_size}.w${window_length}.kcf.log 2>&1"
    touch ${output_dir}/${ref_name}.${sample}.k${kmer_size}.w${window_length}.kcf.done
done < "${data_file}"

# cohort the kcf files
print_log "Cohorting the kcf files"
# remove list file if exists
rm -f ${output_dir}/kcf_list.txt
# create a list file with all the kcf file names
kcf_list_file="${output_dir}/kcf_list.txt"
while IFS=$'\t' read -r sample fq1 fq2; do
    echo "${output_dir}/${ref_name}.${sample}.k${kmer_size}.w${window_length}.kcf" >> ${kcf_list_file}
done < "${data_file}"
# cohort the kcf files
run_cmd "{ time kcftools cohort \
    -l ${kcf_list_file} \
    -o ${output_dir}/${ref_name}.k${kmer_size}.w${window_length}.kcf ; } &> ${output_dir}/${ref_name}.k${kmer_size}.w${window_length}.kcf.log 2>&1"

# find IBS windows
print_log "Finding IBS windows"
run_cmd "{ time kcftools findIBS \
    -i ${output_dir}/${ref_name}.k${kmer_size}.w${window_length}.kcf \
    -o ${output_dir}/${ref_name}.k${kmer_size}.w${window_length}.ibs \
    --summary --min 6 ; } &> ${output_dir}/${ref_name}.k${kmer_size}.w${window_length}.ibs.log 2>&1"

# clean the output directory
print_log "Cleaning the output directory"
rm -f ${output_dir}/kcf_list.txt
# remove the kmc tmp files
rm -rf ${output_dir}/*_kmc_tmp
# fq.list files
rm -f ${output_dir}/*.fq.list

# end clock
end=$(date +%s)
runtime=$((end-start))

print_log "Finished in $(seconds_to_hhmmss $runtime)"
# EOF