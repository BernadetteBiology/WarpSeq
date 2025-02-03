#!/bin/bash

# Author: Bernadette Johnson 
# Created: 13 August 2019
# Last Updated: 31 January 2024
# Purpose: Pipeline for aligning paired end RNAseq reads to a genome and quantifying, with options for user-defined or default configurations, and progress tracking.
# Pipeline: Trimmomatic > HISAT2 > STRINGTIE
# Requirements (Pre-req programs): Trimmomatic, HISAT2, STRINGTIE, SAMTOOLS, GFFREAD, PARALLEL, DOS2UNIX
# Requirements (In directory): 	
  # a. Reads in .fastq format:
	#  a1. forward reads labeled with the ending "_R1.fastq"
	#  a2. reverse reads labeled with the ending "_R2.fastq"
	#  a3. forward and reverse reads should be labeled the exact same apart from R1 or R2.
 # b. The reference genome/chromosome level assembly with a "genome.fna" or "genome.fasta" extension.
 # c. Optional: A chromosome level assembly annotation file with a ".gff" or ".gtf" extension. 
 # d. The adapter file from Trimmomatic. This program currently is set up to use TruSeq3-PE-2.fa.
 # e. The prepDE.py3 script from StringTie.
 # f. This program.
 # g. No other files with these extensions or labels!
 # * All of these files can be labeled with other information, so long as they also adhere to these conditions. 
 # Example contents of a current directory:
  	# mapreadstogenome_forme_HISAT2StringTie.sh
  	# Sample1_ATGCA_fish_R1.fastq
  	# Sample1_ATGCA_fish_R2.fastq
  	# Sample2_GTAGA_fish_R1.fastq
  	# Sample2_GTAGA_fish_R2.fastq
  	# genomefile_fish.fasta
  	# assemblyannotation_fish.gff (optional)
  	# TruSeq3-PE-2.fa
	# prepDE.py3
 # Important:
   	# This program calculates the number of processing units available, and uses 2 less than 50%. If you would like to changes this alter $PARALLEL_JOBS and/or $HISAT2_THREADS. 


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
    echo "  -g, --genome-dir DIR      Directory (output) for HISAT2 genome index (default: $DEFAULT_GENOME_DIR)"
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
    PARALLEL_JOBS=$((((TOTAL_CPUS / 2))-2)) 
    HISAT2_THREADS=$((((TOTAL_CPUS / 2))-2)) 
    if [[ $HISAT2_THREADS -lt 1 ]]; then HISAT2_THREADS=1; fi
    if [[ $PARALLEL_JOBS -lt 1 ]]; then PARALLEL_JOBS=1; fi
    echo "Using $PARALLEL_JOBS parallel jobs and $HISAT2_THREADS threads."
}

# --- Functions ---

# Progress tracker
track_progress() {
    local total=$1
    local message=$2
    local completed=0

    while [[ $completed -lt $total ]]; do
        completed=$(ls "$LOG_DIR" | grep -E "${message}.*log" | wc -l)
        local percentage=$(( (100 * completed) / total ))
        echo -ne "\r$message: $percentage% completed ($completed/$total)"
        sleep 1
    done
    echo -e "\r$message: 100% completed ($completed/$total)"
}

# Validate required tools and input files
validate_inputs() {
    echo "Validating inputs..."
    required_tools=("java" "hisat2" "samtools" "stringtie" "gffread" "parallel" "dos2unix")
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

    echo "All inputs validated successfully."
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

# HISAT2 wrapper for parallel execution
hisat2_process() {
    local base="$1"
    hisat2 --dta -x "$GENOME_DIR/genome_index" -1 "${OUTPUT_DIR}/${base}_R1_paired.fastq" -2 "${OUTPUT_DIR}/${base}_R2_paired.fastq" -S "${OUTPUT_DIR}/${base}.sam" -p "$HISAT2_THREADS" &>> "$LOG_DIR/HISAT2_${base}.log"
}

# Run HISAT2 in parallel, start by building index
run_hisat2() {
    echo "Starting HISAT2..."
    mkdir -p "$GENOME_DIR"
    local genome_file
    genome_file=$(find . -name "*genome.fna" -o -name "*genome.fasta" | head -n 1)

    if [[ -z "$genome_file" ]]; then
        echo "Error: Reference genome not found (*.fna or *.fasta)." >&2
        exit 1
    fi

    cp "$genome_file" "$GENOME_DIR/"
    
	#hisat2-build "$GENOME_DIR/$(basename "$genome_file")" "$GENOME_DIR/genome_index" &>> "$LOG_DIR/HISAT2_build.log"

    local total_samples=$(wc -l < "$TEMP_DIR/basenamereads.txt")
    export -f hisat2_process
    export GENOME_DIR OUTPUT_DIR LOG_DIR HISAT2_THREADS
    parallel -j "$PARALLEL_JOBS" hisat2_process ::: $(cat "$TEMP_DIR/basenamereads.txt")
	#track_progress "$total_samples" "HISAT2 alignments"
    echo "Completed HISAT2." 
}

# Reformat output
convert_process() {
	dos2unix "${OUTPUT_DIR}"/*.sam
}
	
# SAM to BAM and sorting using Samtools
bam_process() {
    local base="$1"
	samtools view -b -S "${OUTPUT_DIR}/${base}.sam" > "${OUTPUT_DIR}/${base}.bam"
    samtools sort "${OUTPUT_DIR}/${base}.bam" -o "${OUTPUT_DIR}/${base}_sorted.bam"
    #rm "${OUTPUT_DIR}/${base}.sam" "${OUTPUT_DIR}/${base}.bam"
	}

process_bam_files() {
    echo "Processing BAM files..."
    while read -r sample; do
        bam_process "$sample"
    done < "$TEMP_DIR/basenamereads.txt"
    echo "BAM processing completed."
}

#Check to see if there is an annotation file
check_gff_and_convert() {
    echo "Checking for .gff file in the current directory..."
    if ls *.gff &> /dev/null; then
        echo "Found .gff file. Converting to .gff3 format..."
        gffread -E *.gff -o- > "${GENOME_DIR}/annotatedgenome_gffread.gff3"
        echo "Conversion to .gff3 completed: ${GENOME_DIR}/annotatedgenome_gffread.gff3"
        USE_ANNOTATION=true
    else
        echo "No .gff file found. Proceeding without annotation."
        USE_ANNOTATION=false
    fi
}

# Run StringTie
run_stringtie() {
    echo "Starting StringTie assembly..."
    mkdir -p "${OUTPUT_DIR}/stringtie"
    while read -r sample; do
        echo "Processing sample: $sample"
        local output_dir="${OUTPUT_DIR}/stringtie/"
        mkdir -p "$output_dir"

        if $USE_ANNOTATION; then
            stringtie "${OUTPUT_DIR}/${sample}_sorted.bam" \
                -o "${output_dir}/${sample}_gtf.gtf" \
                -G "${GENOME_DIR}/annotatedgenome_gffread.gff3" \
                -A "${output_dir}/${sample}_abund.tab" \
                -C "${output_dir}/${sample}_cov_refs.gtf" \
                -p "$HISAT2_THREADS" &>> "$LOG_DIR/StringTie_${base}.log"
        else
            stringtie "${OUTPUT_DIR}/${sample}_sorted.bam" \
                -o "${output_dir}/${sample}_gtf.gtf" \
                -A "${output_dir}/${sample}_abund.tab" \
                -p "$HISAT2_THREADS" &>> "$LOG_DIR/StringTie_${base}.log"
        fi
    done < "$TEMP_DIR/basenamereads.txt"
    echo "StringTie assembly completed."
}

# Run StringTie merge
run_stringtie_merge() {
    echo "Starting StringTie merge..."
    mkdir -p "${OUTPUT_DIR}/stringtie"

    # Generate a list of GTF files
    cd "${OUTPUT_DIR}/stringtie"
	ls *.gtf > gtf_names.txt

	
    if $USE_ANNOTATION; then
		# Update `gtf_names.txt` with required formatting
		sed -i 's/$/ \\/' gtf_names.txt
		sed -i '1s/^/stringtie --merge -G "${GENOME_DIR}/annotatedgenome_gffread.gff3" \\\n-o merged_transcripts.gff \\\n/' gtf_names.txt
		sed -i '$ s/\\//g' gtf_names.txt

		# Run the merge command
		bash gtf_names.txt
		echo "Completed merging of transcripts."
		date

	else
		# Update `gtf_names.txt` with required formatting
		sed -i 's/$/ \\/' gtf_names.txt
		sed -i '1s/^/stringtie --merge \\\n-o merged_transcripts.gff \\\n/' gtf_names.txt
		sed -i '$ s/\\//g' gtf_names.txt

		# Run the merge command
		bash gtf_names.txt
		echo "Completed merging of transcripts."
		date
		
	fi
	
	cd ..
    echo "StringTie merge completed: merged_transcripts.gff"
}

# StringTie after merge for each sample
run_stringtie_after_merge() {
    dos2unix "$TEMP_DIR/basenamereads.txt"
    echo "Running StringTie after merge for each sample..."
    while read -r sample; do
        echo "Processing sample: $sample"
        local output_dir="${OUTPUT_DIR}/stringtie"    
		mkdir -p "${OUTPUT_DIR}/stringtie/${sample}"
		
        stringtie "${OUTPUT_DIR}/${sample}_sorted.bam" -e -B \
            -o "${output_dir}/${sample}/${sample}_merged.gtf" \
            -G "${output_dir}/merged_transcripts.gff" \
            -A "${output_dir}/${sample}/${sample}_abund_merged.tab" \
            -C "${output_dir}/${sample}/${sample}_cov_refs_merged.gtf"
    done < "$TEMP_DIR/basenamereads.txt"  &>> "$LOG_DIR/StringTieAfterMerge_${base}.log"
    echo "Completed StringTie after merge."
}

# Generate counts table using prepDE.py3
generate_counts_table() {
    echo "Generating counts matrix using prepDE.py3..."
    if [[ ! -f "prepDE.py3" ]]; then
        echo "Error: prepDE.py3 not found in the current directory." >&2
        exit 1
    fi
	
	cp "${TEMP_DIR}/basenamereads.txt" sample_lst.txt
	dos2unix sample_lst.txt
	
    sed "s/$/\ /" "${TEMP_DIR}/basenamereads.txt" > sample_lst.txt
	
    while read -r sample; do
	echo ${sample}
	sed -i "s|^${sample}|${sample} ${OUTPUT_DIR}/stringtie/${sample}/${sample}_merged.gtf|" sample_lst.txt

    done < "${TEMP_DIR}/basenamereads.txt"

    python3 prepDE.py3 -i sample_lst.txt
    echo "Counts matrix generated."

}

# Cleanup and organize outputs
cleanup_and_organize() {
    echo "Organizing outputs..."
    mkdir -p "$OUTPUT_DIR/sorted_bam"
    mv "$OUTPUT_DIR/"*_sorted.bam "$OUTPUT_DIR/sorted_bam/"
    rm -rf "$TEMP_DIR"
    echo "Output files organized."
}

# --- Main Script Execution ---
main() {
    calculate_resources
    validate_inputs
    prepare_files
	#run_trimmomatic
    run_hisat2
	convert_process
    process_bam_files
	check_gff_and_convert
    run_stringtie
	run_stringtie_merge
    run_stringtie_after_merge
    generate_counts_table
    cleanup_and_organize
    echo "Pipeline completed. Results are in $OUTPUT_DIR."
}

main "$@"
