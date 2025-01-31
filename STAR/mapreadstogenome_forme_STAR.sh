#!/bin/bash

# Author: Bernadette Johnson 
# Created: 13 August 2019
# Last Updated: 30 January 2024
# Purpose: Pipeline for aligning paired end RNAseq reads to a genome and quantifying, with options for user-defined or default configurations, and progress tracking.
# Pipeline: Trimmomatic > STAR
# Requirements (Pre-req programs): Trimmomatic, STAR
# Requirements (In directory): 	
  # a. Reads in .fastq format:
	#  a1. forward reads labeled with the ending "_R1.fastq"
	#  a2. reverse reads labeled with the ending "_R2.fastq"
	#  a3. forward and reverse reads should be labeled the exact same apart from R1 or R2.
 # b. The reference genome/chromosome level assembly with a "genome.fna" or "genome.fasta" extension.
 # c. Having a chromosome level assembly annotation file with a ".gff" or ".gtf" extension. 
 # d. The adapter file from Trimmomatic. This program currently is set up to use TruSeq3-PE-2.fa.
 # e. This program.
 # f. No other files with these extensions or labels!
 # * All of these files can be labeled with other information, so long as they also adhere to these conditions. 
 # Example contents of a current directory:
  	# mapreadstogenome_forme_STAR.sh
  	# Sample1_ATGCA_fish_R1.fastq
  	# Sample1_ATGCA_fish_R2.fastq
  	# Sample2_GTAGA_fish_R1.fastq
  	# Sample2_GTAGA_fish_R2.fastq
  	# genomefile_fish.fasta
  	# assemblyannotation_fish.gff
  	# TruSeq3-PE-2.fa
 # Important:
   	# This program calculates the number of processing units available, and uses 2 less than 50%. If you would like to changes this alter $PARALLEL_JOBS and/or $STAR_THREADS. 

# --- Default Configurations ---
DEFAULT_TRIMMOMATIC_JAR="./trimmomatic-0.39.jar"
DEFAULT_ILLUMINACLIP="TruSeq3-PE-2.fa:2:30:10:2:keepBothReads"
DEFAULT_TRIM_OPTS="LEADING:3 TRAILING:3 MINLEN:75"
DEFAULT_GENOME_DIR="./genome_index"
DEFAULT_OUTPUT_DIR="./results"
DEFAULT_TEMP_DIR="./temp"
DEFAULT_LOG_DIR="./logs"

# --- Parse Command-Line Arguments ---
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -t, --trimmomatic PATH    Path to Trimmomatic JAR file (default: $DEFAULT_TRIMMOMATIC_JAR)"
    echo "  -i, --illumina            Trimmomatic adapter clip parameters (default: $DEFAULT_ILLUMINACLIP)"
    echo "  -r, --trim                Trimmomatic trimming parameters (default: $DEFAULT_TRIM_OPTS)"
    echo "  -g, --genome-dir DIR      Directory (output) for STAR genome index (default: $DEFAULT_GENOME_DIR)"
    echo "  -o, --output-dir DIR      Directory (output) for output files (default: $DEFAULT_OUTPUT_DIR)"
    echo "  -l, --log-dir DIR         Directory for log files (default: $DEFAULT_LOG_DIR)"
    echo "  -h, --help                Display this help message"
}

# Default values
TRIMMOMATIC_JAR="$DEFAULT_TRIMMOMATIC_JAR"
ILLUMINACLIP="$DEFAULT_ILLUMINACLIP"
TRIM_OPTS="$DEFAULT_TRIM_OPTS"
GENOME_DIR="$DEFAULT_GENOME_DIR"
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
TEMP_DIR="$DEFAULT_TEMP_DIR"
LOG_DIR="$DEFAULT_LOG_DIR"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--trimmomatic)
            TRIMMOMATIC_JAR="$2"
            shift 2
            ;;
        -i|--illumina)
            ILLUMINACLIP="$2"
            shift 2
            ;;
		-r|--trim)
            TRIM_OPTS="$2"
            shift 2
            ;;
        -g|--genome-dir)
            GENOME_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -l|--log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# --- Dynamic Resource Allocation ---
calculate_resources() {
    echo "Calculating system resources..."
    TOTAL_CPUS=$(nproc)
    PARALLEL_JOBS=$((TOTAL_CPUS / 2)-2) 
    STAR_THREADS=$((TOTAL_CPUS / 2)-2) 
    if [[ $STAR_THREADS -lt 1 ]]; then STAR_THREADS=1; fi
    if [[ $PARALLEL_JOBS -lt 1 ]]; then PARALLEL_JOBS=1; fi
    echo "Using $PARALLEL_JOBS parallel jobs and $STAR_THREADS STAR threads."
}

# --- Functions ---

# Progress tracker
track_progress() {
    local total=$1    	# Total number of tasks to track
    local marker=$2 	# Completion marker to search for in the log file
    local completed=0 	# Initialize completed tasks to 0

    while [[ $completed -lt $total ]]; do
		# Count occurrences of the completion marker in the log file
        completed=$(ls "$LOG_DIR" | grep -E "${marker}.*log" | wc -l)
        local percentage=$(( (100 * completed) / total ))
        echo -ne "\r$marker: $percentage% completed ($completed/$total)"
        sleep 1
    done
    echo -e "\r$marker: 100% completed ($completed/$total)"
}

# Validate required tools and input files
validate_inputs() {
    echo "Validating inputs..."
    required_tools=("java" "STAR" "samtools" "stringtie" "gffread" "parallel")
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &>/dev/null; then
            echo "Error: $tool is not installed or not in PATH." >&2
            exit 1
        fi
    done

    if ! [ -f "$TRIMMOMATIC_JAR" ]; then
        echo "Error: Trimmomatic JAR not found at $TRIMMOMATIC_JAR." >&2
        exit 1
    fi

    ls *R1.fastq &>/dev/null || { echo "Error: No forward reads found (*.R1.fastq)." >&2; exit 1; }
    ls *R2.fastq &>/dev/null || { echo "Error: No reverse reads found (*.R2.fastq)." >&2; exit 1; }

	#Find the reference genome file 
    genome_file=$(find . -name "*genome.fna" -o -name "*genome.fasta" | head -n 1)

    if [[ -z "$genome_file" ]]; then
        echo "Error: Reference genome not found (*genome.fna or *genome.fasta)." >&2
        exit 1
    fi

    echo "All inputs validated successfully."
}

#Check to see if there is an annotation file
check_annotation() {
    echo "Checking for an annotation file in the current directory (.gff or .gtf)..."
	local annotation
	annotation=$(find . -name "*.gtf" -o -name "*.gff" | head -n 1)
    
if [[ -n "$annotation" ]]; then
        echo "Found annotation file."
        USE_ANNOTATION=true
    else
        echo "No annotation file found. I am currently not set up to run STAR without an annotation."
		exit 1
        USE_ANNOTATION=false
    fi
}

# Prepare input files
prepare_files() {
    echo "Preparing input files..."
    mkdir -p "$TEMP_DIR" "$OUTPUT_DIR" "$LOG_DIR"
    ls *R1.fastq > "$TEMP_DIR/forwardreads.txt"
    ls *R2.fastq > "$TEMP_DIR/reversereads.txt"
    sed 's/_R1.fastq//' "$TEMP_DIR/forwardreads.txt" > "$TEMP_DIR/basenamereads.txt"
}

# Trimmomatic wrapper for parallel execution
trimmomatic_process() {
    local base="$1"
    java -jar "$TRIMMOMATIC_JAR" PE \
        "${base}_R1.fastq" "${base}_R2.fastq" \
        "${OUTPUT_DIR}/${base}_R1_paired.fastq" "${OUTPUT_DIR}/${base}_R1_unpaired.fastq" \
        "${OUTPUT_DIR}/${base}_R2_paired.fastq" "${OUTPUT_DIR}/${base}_R2_unpaired.fastq" \
        ILLUMINACLIP:$ILLUMINACLIP $TRIM_OPTS &>> "$LOG_DIR/Trimmomatic_${base}.log"
}

# Run Trimmomatic in parallel
run_trimmomatic() {
    echo "Starting Trimmomatic..."
    local total_samples=$(wc -l < "$TEMP_DIR/basenamereads.txt")
    export -f trimmomatic_process
    export TRIMMOMATIC_JAR ILLUMINACLIP TRIM_OPTS OUTPUT_DIR LOG_DIR
    parallel -j "$PARALLEL_JOBS" trimmomatic_process ::: $(cat "$TEMP_DIR/basenamereads.txt")
    track_progress "$total_samples" "Trimmomatic"
    echo "Completed Trimmomatic."
}

# Run STAR in parallel, start by building index
run_STARindexbuild() {
    echo "Starting STAR..."
    mkdir -p "$GENOME_DIR"
	
    cp "$genome_file" "$GENOME_DIR/"
	
	#Find the reference annotation file and put it into a folder to build index
	if $USE_ANNOTATION; then
		echo Using annotation.
		local annotation
		annotation=$(find . -name "*.gtf" -o -name "*.gff" | head -n 1)
    	cp "$annotation" "$GENOME_DIR/"	
		
		#Build STAR genomic index using annotation
		STAR --runThreadN "$STAR_THREADS" --runMode genomeGenerate --genomeDir "$GENOME_DIR" --genomeFastaFiles "$GENOME_DIR/$(basename "$genome_file")" --sjdbGTFfile "$GENOME_DIR/$(basename "$annotation")" &>> "$LOG_DIR/STAR_build.log"
	
	else
		echo "Not using annotation. I am currently not set up to run STAR without an annotation."
		exit 1
	fi

    echo "Completed STAR index building."
}

# Wrapper for STAR alignment of reads and counting
countem_up() {
    local base="$1"
	mkdir -p "${OUTPUT_DIR}/${base}"
	STAR --genomeDir "$GENOME_DIR" --readFilesIn "${OUTPUT_DIR}/${base}_R1_paired.fastq" "${OUTPUT_DIR}/${base}_R2_paired.fastq" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix "${OUTPUT_DIR}/${base}/${base}" --runThreadN "$STAR_THREADS" &>> "$LOG_DIR/STAR_align_and_count.log"
	
	}

# Run STAR
run_countem_up() {
    echo "Starting the alignment and counting of reads with STAR..."
    #local total_samples=$(wc -l < "$TEMP_DIR/basenamereads.txt")
    export -f countem_up
    export GENOME_DIR OUTPUT_DIR STAR_THREADS LOG_DIR
    parallel -j "$PARALLEL_JOBS" countem_up ::: $(cat "$TEMP_DIR/basenamereads.txt")
    track_progress "$total_samples" "STAR align and count progress"
    echo "Completed STAR count."
}

# --- Main Script Execution ---
main() {
    calculate_resources
    validate_inputs
    check_annotation
    prepare_files
    run_trimmomatic
    run_STARindexbuild
    run_countem_up
    echo "Done. Results are in $OUTPUT_DIR."
}

main "$@"
