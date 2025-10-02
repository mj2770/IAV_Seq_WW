#!/bin/bash
#SBATCH -t 00-24:00:00
#SBATCH -N 1
#SBATCH -o stdout_Twist_IAV_%j
#SBATCH -e stderr_Twist_IAV_%j
#SBATCH --job-name preview_seq
#SBATCH --mail-user mj2770@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -A careswwm

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if not submitting to SLURM

####################################################################

Combined_Twist_IAV.sh
created 12/28/2024 by Minxi Jiang

####################################################################
#
# This script performs the overall analysis of sequencing data from Twist_IAV.
# It starts with  QC trimming, bowtie2 reads mapping and filtering,and post-mapping analysis.
#
#   1. QC trimming of the raw data
#  
#    1.1 Read cleaning with fastp 
#        Trims adaptors, polyX, polyG tails, minimum length 70, and base quality 15 (default)
#    1.2 [SKIP THIS STEP] Human reads removal
#    1.3 Deduplication using seqkit
#    1.4 Pairing the deduplicated reads and storing the unpaired reads
#    1.5 Statistical analysis of the raw, cleaned, and deduplicated reads using seqkit stats
#
#   2. Bowtie2 mapping and filtering of cleaned reads
#
#    2.1 Build-up reference genome index
#    2.2 Paired reads mapping using bowtie2
#        only output mapped reads, no discordantly mapped reads, no single mapped reads (must be pairs)
#    2.3 Converted to sorted bam file
#    2.4 (Optional) Filter the mapped reads by reformat editfilter Notes: cannot add verifypaired=t, otherwise there's problem running reformat.sh
#        Filter by the sum of indels, deletions, insertions, substitutions
#    2.5 repaired, sorted, indexed the filtered bam file
#   
   3. Post-mapping analys
#
#    3.1 Index the reference genome and filtered sorted bam files from bowite2 mapping
#    3.2 Generate the position-by-position coverage depth file
#    3.3 Generate the overall summarized coverage breadth and depth file for each reference genome (reconfirm min length 75)
#    3.4 Generate the reads mapping summary
#    3.5 Generate the vcf file
#        (Optional) Filter the vcf file and generate consensus sequences
#    3.6 Generate the ivar variants calling file with filtering options
####################################################################
####################################################################
# Set the PATH 
export PATH=/p/lustre1/preview/seq_software/sambamba/:/p/lustre1/preview/seq_software/bbmap:/p/lustre1/preview/seq_software/:/p/lustre1/preview/seq_software/bowtie2-2.5.1/:/p/lustre1/preview/seq_software/samtools-1.17/:$PATH
####################################################################
# STEP 1: QC trim

# Define directories
INPUT_DIR=/p/vast1/preview/241226_twist_tiled_ucb_ww_IAV/sample_4
OUT_DIR=/p/lustre1/preview/IAV_seq_data/sample_4
# Optional: filter the human reads
# Ref=/p/lustre1/preview/seq_data/output/STEP_1_1_Human_scrubber/Human_genome_db

# Create subdirectories for each step
mkdir -p "$OUT_DIR/step_1_fastp"
# mkdir -p "$OUT_DIR/step_2_human_reads"
mkdir -p "$OUT_DIR/step_3_deduplication"
mkdir -p "$OUT_DIR/step_4_paired"
mkdir -p "$OUT_DIR/step_5_stats"

# Log file
LOG_FILE="$OUT_DIR/debug_twist_comp.log"
echo "Starting pipeline at $(date)." > "$LOG_FILE"

# Record the start time of the step
start_time=$(date +%s)

###################################
# STEP 1.1: QC Trim with fastp
echo "Checking fastp step..." >> "$LOG_FILE"
if [ ! "$(ls -A "$OUT_DIR/step_1_fastp")" ]; then
    echo "Running fastp step..." >> "$LOG_FILE"
    for f1 in "$INPUT_DIR"/*R1*fastq.gz; do
        f2=$(echo "$f1" | awk -F "R1" '{print $1 "R2" $2}')
        output_f1="$OUT_DIR/step_1_fastp/$(basename "$f1" .fastq.gz)_filt.fastq.gz"
        output_f2="$OUT_DIR/step_1_fastp/$(basename "$f2" .fastq.gz)_filt.fastq.gz"
        fastp \
            -g --poly_g_min_len=5 -x --poly_x_min_len=5 -y \
            -l 75 \
            -i "$f1" -I "$f2" \
            -o "$output_f1" \
            -O "$output_f2" \
            -h "$OUT_DIR/step_1_fastp/$(basename "$f1" .fastq.gz)_fastp.html" \
            -j "$OUT_DIR/step_1_fastp/$(basename "$f1" .fastq.gz)_fastp.json"
    done
    echo "fastp step completed." >> "$LOG_FILE"
else
    echo "fastp step already completed." >> "$LOG_FILE"
fi

###################################
# STEP 1.3: Remove duplicates with seqkit
echo "Checking seqkit deduplication step..." >> "$LOG_FILE"
if [ ! "$(ls -A "$OUT_DIR/step_3_deduplication")" ]; then
    echo "Running seqkit deduplication step..." >> "$LOG_FILE"
    for filt_file in "$OUT_DIR/step_1_fastp"/*R1_001_filt.fastq.gz; do
        filename=$(basename "$filt_file" _R1_001_filt.fastq.gz)
        seqkit rmdup -s "$filt_file" -o "$OUT_DIR/step_3_deduplication/$filename.forward.deduped.fastq.gz"
        reverse_file="${filt_file/R1/R2}"
        seqkit rmdup -s "$reverse_file" -o "$OUT_DIR/step_3_deduplication/$filename.reverse.deduped.fastq.gz"
    done
    echo "seqkit deduplication step completed." >> "$LOG_FILE"
else
    echo "seqkit deduplication step already completed." >> "$LOG_FILE"
fi

###################################
# STEP 1.4: Pair deduplicated reads and store unpaired reads
echo "Checking pairing step..." >> "$LOG_FILE"
if [ ! "$(ls -A "$OUT_DIR/step_4_paired")" ]; then
    echo "Running pairing step..." >> "$LOG_FILE"
    for forward_file in "$OUT_DIR/step_3_deduplication/"*.forward.deduped.fastq.gz; do
        reverse_file="${forward_file/forward/reverse}"
        seqkit pair -1 "$forward_file" -2 "$reverse_file" -O "$OUT_DIR/step_4_paired" -u
    done
    echo "Pairing step completed." >> "$LOG_FILE"
else
    echo "Pairing step already completed." >> "$LOG_FILE"
fi

###################################
# STEP 1.5: Generate statistics for raw, cleaned, and deduplicated reads
echo "Checking stats step..." >> "$LOG_FILE"
if [ ! -f "$OUT_DIR/step_5_stats/stats_raw.txt" ]; then
    echo "Running stats step..." >> "$LOG_FILE"
    seqkit stats -T -j 100 "$INPUT_DIR"/*fastq.gz -o "$OUT_DIR/step_5_stats/stats_raw.txt"
    seqkit stats -T -j 100 "$OUT_DIR/step_1_fastp"/*fastq.gz -o "$OUT_DIR/step_5_stats/stats_fastq.txt"
    seqkit stats -T -j 100 "$OUT_DIR/step_3_deduplication"/*fastq.gz -o "$OUT_DIR/step_5_stats/stats_dedupe.txt"
    seqkit stats -T -j 100 "$OUT_DIR/step_4_paired"/*fastq.gz -o "$OUT_DIR/step_5_stats/stats_paired.txt"
    echo "Stats step completed." >> "$LOG_FILE"
else
    echo "Stats step already completed." >> "$LOG_FILE"
fi

###################################
# Log the time taken for Step 1
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Step 1 completed in $elapsed_time seconds." >> "$LOG_FILE"
echo "Pipeline completed at $(date)." >> "$LOG_FILE"

exit 0


########################################################################
# STEP 2: bowtie2 alignment
# Define input and output paths
ref_genomes="/p/lustre1/preview/influenza/Spike_virus_genomes/same_length/remove_gaps/cut_ref_S1S2S3.fasta"
input="/p/lustre1/preview/IAV_seq_data/twist_IAV_full/step_1_fastp"
output="/p/lustre1/preview/IAV_seq_data/twist_IAV_full/STEP_2_bowtie2_nondedup"

# STEP 2.1 Create Bowtie2 index for the reference genome
bowtie2-build $ref_genomes ${input}/ref_genomes
samtools faidx $ref_genomes

# Loop over paired-end reads
for f in ${input}/*R1_001_filt.fastq.gz; do
  # Find the corresponding reverse read
  r=$(echo $f | sed 's/R1/R2/')

  # Extract the base name for the reads
  read_base=$(basename $f _L001_R1_001_filt.fastq.gz)

  # STEP 2.2 Align paired-end reads to the reference genome
  bowtie2 -p 100 -x ${input}/ref_genomes --no-unal --no-discordant --no-mixed -1 $f -2 $r -S ${output}/${read_base}.sam

  # STEp 2.3 Convert SAM to BAM, sort, and index BAM file
  samtools view -bS ${output}/${read_base}.sam | samtools sort -o ${output}/${read_base}.sorted.bam
  samtools index ${output}/${read_base}.sorted.bam

  # Clean up the SAM file to save space
  rm ${output}/${read_base}.sam

  # STEP 2.4 Filter mapped paired-end reads with fewer than 5 edits
  reformat.sh in=${output}/${read_base}.sorted.bam out=${output}/${read_base}.filtered.bam editfilter=5

 # STEP 2.5 Filter, discard unpaired reads, and sort the BAM file
  samtools view -h -f 1 -b ${output}/${read_base}.filtered.bam | \
  samtools sort -o ${output}/${read_base}.filtered.paired.sorted.bam

  # Clean up intermediate BAM file
  rm ${output}/${read_base}.filtered.bam
done
##############################################################################
# STEP 3. Post-mapping analysis
# Add software paths to the PATH environment variable
export PATH=/p/lustre1/preview/seq_software/bin/ivar/:/p/lustre1/preview/seq_software/bcftools/bcftools/:/p/lustre1/preview/seq_software/:/p/lustre1/preview/seq_software/bowtie2-2.5.1/:/p/lustre1/preview/seq_software/samtools-1.17/:$PATH

# Set up reference and output directories
input_dir="/p/lustre1/preview/IAV_seq_data/twist_IAV_full/STEP_2_bowtie2_nondedup/"
output_dir="/p/lustre1/preview/IAV_seq_data/twist_IAV_full/STEP_3_depth_breadth_nondedup/"
ref="/p/lustre1/preview/influenza/Spike_virus_genomes/same_length/remove_gaps/cut_ref_S1S2S3.fasta"

# Index the reference genome
samtools faidx "$ref"

# Loop through BAM files in the input directory
for sorted_bam_file in "$input_dir"*.filtered.paired.sorted.bam; do

    # STEP 3.1: Index the sorted BAM file
    samtools index "$sorted_bam_file" || { echo "Error: Failed to index BAM file $sorted_bam_file"; exit 1; }

    # Extract the basename of the BAM file for output file naming
    basename=$(basename "$sorted_bam_file" .filtered.paired.sorted.bam)

    # STEP 3.2: Generate the depth file
    echo "Generating coverage depth file for $basename..."
    samtools depth -a "$sorted_bam_file" > "${output_dir}${basename}_coverage_depth.txt" || { echo "Error: Failed to generate depth file for $basename"; exit 1; }

    # STEP 3.3: Generate the coverage breadth file with a minimum read length of 75
    echo "Generating coverage breadth file with minimum read length of 75 for $basename..."
    samtools coverage -l 75 "$sorted_bam_file" > "${output_dir}${basename}_coverage_breadth.txt" || { echo "Error: Failed to generate coverage breadth file for $basename"; exit 1; }

    # STEP 3.4: Generate the mapping summary
    echo "Generating mapping summary for $basename..."
    samtools stats "$sorted_bam_file" > "${output_dir}${basename}_stats.txt" || { echo "Error: Failed to generate stats for $basename"; exit 1; }

    # STEP 3.5: Run bcftools mpileup to generate a VCF file
    echo "Generating VCF file using bcftools for $basename..."
    vcf_output_file="${output_dir}${basename}_vcf.gz"

    bcftools mpileup -a FORMAT/AD,FORMAT/DP,FORMAT/SP -d 1000000 -f "$ref" "$sorted_bam_file" | \
    bcftools call -mv -Oz -o "$vcf_output_file" || { echo "Error: Failed to call variants with bcftools for $basename"; exit 1; }

    # (Optional) Filter the raw VCF file
    # echo "Filtering VCF file..."
    # bcftools filter -i 'INFO/AC / INFO/AN >= 0.1 & FMT/DP > 10 & QUAL > 10' \
    #     -o "$output_dir/filtered.vcf.gz" \
    #     "$vcf_output_file" || { echo "Error: Failed to filter VCF file"; exit 1; }

    # (Optional) Index the filtered VCF file
    # echo "Indexing the filtered VCF file..."
    # bcftools index "$output_dir/filtered.vcf.gz" || { echo "Error: Failed to index the filtered VCF file"; exit 1; }

    # (Optional) Generate a consensus sequence
    # echo "Generating consensus sequence..."
    # bcftools consensus -f "$ref" "$output_dir/filtered.vcf.gz" \
    #     > "$output_dir/consensus.fasta" || { echo "Error: Failed to generate consensus sequence"; exit 1; }

    # STEP 3.6: Run samtools mpileup and iVar for variant calling
    # echo "Calling variants with iVar for $basename..."
    # ivar_output_file="${output_dir}${basename}_ivar_call.txt"

    # samtools mpileup -aa -A -d 0 -B -Q 0 -f "$ref" "$sorted_bam_file" | \
    # ivar variants -p "$ivar_output_file" -q 20 -m 10 -r "$ref" || { echo "Error: Failed to call variants with iVar for $basename"; exit 1; }

done
