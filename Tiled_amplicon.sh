#!/bin/bash
#SBATCH -t 00-24:00:00
#SBATCH -N 1
#SBATCH -o stdout_trim_cutadapt_%j
#SBATCH -e stderr_trim_cutadapt_%j
#SBATCH --job-name qc-cutadapt

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if not submitting to SLURM

####################################################################

Combined_tiled_amplicon.sh
created by Rose Kantor and  Minxi Jiang

####################################################################
#
# This script performs the overall analysis of sequencing data from tiled-amplicon sequencing# of IAV. It starts with QC trim,primer and adaptor removal, reads statistic analysis to bowtie2 alignment and the post-mapping analysis
#
#   1. QC trimming of the raw data
#  
#    1.1 Read cleaning with fastp 
#        Trims adaptors, polyX, polyG tails, minimum length 70, and base quality 15 (default)
#    1.2 Primers and adaptors removal using cutadapt
#    1.3 Statistical analysis of the raw, filtered, and after-cutadapt reads using seqkit stats
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
STEP 1: QC Trim
export PATH=/p/lustre1/preview/seq_software/samtools-1.17/:/p/lustre1/preview/seq_software/bowtie2-2.5.1/:/p/lustre1/preview/seq_software/:/p/lustre1/preview/seq_software/bbmap/:/p/lustre1/preview/seq_software/bedtools2/bin/:$PATH

INPUT_DIR=/p/vast1/preview/241226_twist_tiled_ucb_ww_IAV/Tiled_amplicon
OUT_DIR=/p/lustre1/preview/IAV_seq_data/tiled_amp_full
# note: the primer sequences are currently in OUT_DIR, change as needed

# run fastp
for f1 in "$INPUT_DIR"/*R1*fastq.gz; do
  f2=$(echo "$f1" | awk -F "R1" '{print $1 "R2" $2}')
  output_f1="$OUT_DIR/step_1_qc/$(basename "$f1" _001.fastq.gz)_filt.fastq.gz"
  output_f2="$OUT_DIR/step_1_qc/$(basename "$f2" _001.fastq.gz)_filt.fastq.gz"
  fastp \
  -g --poly_g_min_len=5 -x --poly_x_min_len=5 -y \
  -l 70 \
  -i "$f1" -I "$f2" \
  -o "$output_f1" \
  -O "$output_f2" \
  -h "$OUT_DIR/step_1_qc/$(basename "$f1" .fastq.gz)_fastp.html" \
  -j "$OUT_DIR/step_1_qc/$(basename "$f1" .fastq.gz)_fastp.json"
done

# run cutadapt
conda activate cutadapt

for f1 in $OUT_DIR/step_1_qc/*R1*_filt.fastq.gz; do
  f2=$(echo "$f1" | awk -F "R1" '{print $1 "R2" $2}')

  output_f1="$OUT_DIR/step_1_qc/$(basename "$f1" _filt.fastq.gz)_cutadapt.fastq.gz"
  output_f2="$OUT_DIR/step_1_qc/$(basename "$f2" _filt.fastq.gz)_cutadapt.fastq.gz"

  cutadapt --pair-adapters --discard-untrimmed \
  -g^file:$OUT_DIR/primers_F.fasta \
  -G^file:$OUT_DIR/primers_R.fasta -m 50 \
  -o $output_f1 \
  -p $output_f2 \
  $f1 \
  $f2 \

done

conda deactivate

# get stats
seqkit stats -T -j 100 $INPUT_DIR/*fastq.gz -o $OUT_DIR/stats/stats_raw.txt
seqkit stats -T -j 100 $OUT_DIR/step_1_qc/*filt.fastq.gz -o $OUT_DIR/stats/stats_fastp.txt
seqkit stats -T -j 100 $OUT_DIR/step_1_qc/*cutadapt.fastq.gz -o $OUT_DIR/stats/stats_cutadapt.txt
###################################################################################################
STEP 2: Bowtie2 mapping

export PATH=/p/lustre1/preview/seq_software/:/p/lustre1/preview/seq_software/bowtie2-2.5.1/:/p/lustre1/preview/seq_software/samtools-1.17/:$PATH

# Define input and output paths
ref_genomes="/p/lustre1/preview/influenza/Spike_virus_genomes/same_length/remove_gaps/cut_ref_S1S2S3.fasta"
input=/p/lustre1/preview/IAV_seq_data/tiled_amp_test4/step_1_qc
output=/p/lustre1/preview/IAV_seq_data/tiled_amp_test4/step_2_bowtie/

# Create Bowtie2 index for the reference genome
bowtie2-build $ref_genomes ${output}/ref_genomes
samtools faidx $ref_genomes

# Loop over paired-end reads
for f in ${input}/*R1*_cutadapt.fastq.gz; do
  # Find the corresponding reverse read
  r=$(echo $f | sed 's/R1/R2/')

  # Extract the base name for the reads
  read_base=$(basename $f _cutadapt.fastq.gz)

  # Align paired-end reads to the reference genome
  bowtie2 -p 100 -x ${output}/ref_genomes --no-unal --no-discordant --no-mixed -1 $f -2 $r -S ${output}/${read_base}.sam
  
  # Convert SAM to BAM, sort, and index BAM file
  samtools view -bS ${output}/${read_base}.sam | samtools sort -o ${output}/${read_base}.sorted.bam
  samtools index ${output}/${read_base}.sorted.bam

  # Clean up the SAM file to save space
  rm ${output}/${read_base}.sam
done
##################################################################################################
STEP 3: Post-mapping analysis


# Set up reference and output directories
input_dir="/p/lustre1/preview/IAV_seq_data/tiled_amp_test4/step_2_bowite/"
output_dir="/p/lustre1/preview/IAV_seq_data/tiled_amp_test4/step_3_depth_breadth/"
ref="/p/lustre1/preview/influenza/Spike_virus_genomes/same_length/remove_gaps/cut_ref_S1S2S3.fasta"
mkdir -p "$output_dir"

# Index the reference genome
samtools faidx "$ref"

# Loop through BAM files in the input directory
for sorted_bam_file in "$input_dir"*.sorted.bam; do

    # STEP 1: Index the sorted BAM file
    samtools index "$sorted_bam_file" || { echo "Error: Failed to index BAM file $sorted_bam_file"; exit 1; }

    # Extract the basename of the BAM file for output file naming
    basename=$(basename "$sorted_bam_file" .sorted.bam)

    # STEP 2: Generate the depth file
    echo "Generating coverage depth file for $basename..."
    samtools depth -a "$sorted_bam_file" > "${output_dir}${basename}_coverage_depth.txt" || { echo "Error: Failed to generate depth file for $basename"; exit 1; }

    # STEP 3: Generate the coverage breadth file with a minimum read length of 75
    echo "Generating coverage breadth file with minimum read length of 75 for $basename..."
    samtools coverage -l 75 "$sorted_bam_file" > "${output_dir}${basename}_coverage_breadth.txt" || { echo "Error: Failed to generate coverage breadth file for $basename"; exit 1; }

    # STEP 4: Generate the mapping summary
    echo "Generating mapping summary for $basename..."
    samtools stats "$sorted_bam_file" > "${output_dir}${basename}_stats.txt" || { echo "Error: Failed to generate stats for $basename"; exit 1; }

    # STEP 5: Run bcftools mpileup to generate a VCF file
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

    # STEP 6: Run samtools mpileup and iVar for variant calling
    echo "Calling variants with iVar for $basename..."
    ivar_output_file="${output_dir}${basename}_ivar_call.txt"

    samtools mpileup -aa -A -d 0 -B -Q 0 -f "$ref" "$sorted_bam_file" | \
    ivar variants -p "$ivar_output_file" -q 20 -m 10 -r "$ref" || { echo "Error: Failed to call variants with iVar for $basename"; exit 1; }

done

