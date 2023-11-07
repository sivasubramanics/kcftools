#!/bin/bash
# script to call variations using IBSpy
# contact: siva.selvanayagam[at]wur.nl
# version: 0.1
# date: 02-08-2023

set -u
set -e

start=$(date +%s)

# functin to print log messages
print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] $1"
}

# function to run commands
run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

# function to convert seconds to hh:mm:ss format
seconds_to_hhmmss(){
    local seconds=$1
    local hours=$(printf "%02d" $((seconds / 3600)))
    local minutes=$(printf "%02d" $(((seconds / 60) % 60)))
    local seconds=$(printf "%02d" $((seconds % 60)))
    echo "$hours:$minutes:$seconds"
}

# function to print usage
usage(){
    echo "Usage: $0 [-h] [-f fq.list] [-F fq_list] [-r ref.fa] [-o output] [-p threads] [-k kmer_size] [-w window_size]"
    echo "Options:"
    echo "  -f STR    list of fastq files"
    echo "  -F STR    format (fq/fq_list/fa/fa_list/)"
    echo "  -r STR    reference genome"
    echo "  -o STR    output file"
    echo "  -p INT    number of threads"
    echo "  -k INT    kmer size"
    echo "  -w INT    window size"
    echo "  -h        help"
    exit 1
}

# function to check if a command exists
check_exec(){
    if ! command -v $1 &> /dev/null; then
        echo "Error: $1 is not installed"
        exit 1
    fi
}

# parse arguments
while getopts "f:F:r:o:p:k:w:h" opt; do
    case $opt in
        f) input=$OPTARG;;
        F) format=$OPTARG;;
        r) ref=$OPTARG;;
        o) output=$OPTARG;;
        p) threads=$OPTARG;;
        k) kmer_size=$OPTARG;;
        w) window_size=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

# check arguments and input files
if [ $# -eq 0 ] || [ -z "$input" ] || [ -z "$format" ] || [ -z "$ref" ] || [ -z "$output" ] || [ -z "$threads" ] || [ -z "$kmer_size" ] || [ -z "$window_size" ]; then
    usage
fi

# check if input files exist
if [ "$format" == "fq" ]; then
    in_cmd="-fq $input"
elif [ "$format" == "fq_list" ]; then
    in_cmd="-fq @$input"
elif [ "$format" == "fa" ]; then
    in_cmd="-fa $input"
elif [ "$format" == "fa_list" ]; then
    in_cmd="-fa @$input"
else
    echo "Error: format should be fq/fq_list/fa/fa_list"
    exit 1
fi

# check if the kmersGWAS tools are installed
print_log "Checking tools"
check_exec kmc 
check_exec $TOOLS_PATH/kmersGWAS/bin/kmers_add_strand_information

# pipeline begins here
print_log "Creating output directory $output"
mkdir -p $output

print_log "KMC canonized counting"
run_cmd "kmc -k$kmer_size -t$threads -ci0 ${in_cmd} $output/canon $output" > $output/kmc.canon.log 2>&1

print_log "KMC no canonized counting (all k-mers)"
run_cmd "kmc -k$kmer_size -t$threads -ci0 -b ${in_cmd} $output/all $output" > $output/kmc.all.log 2>&1

print_log "Adding strand information"
run_cmd "$TOOLS_PATH/kmersGWAS/bin/kmers_add_strand_information -c $output/canon -n $output/all -k $kmer_size -o $output/kmc31" > $output/kmers_add_strand_information.log 2>&1

print_log "IBScpp variant calling"
run_cmd "$TOOLS_PATH/IBSpy/IBScpp/build/IBScpp -d $output/kmc31 -r $ref -p $threads -k $kmer_size -w $window_size > $output/variations.tsv" > $output/IBScpp.variations.log 2>&1

# stop clock
end=$(date +%s)

# print log
print_log "elapsed time: $(seconds_to_hhmmss $((end - start)))"

# END