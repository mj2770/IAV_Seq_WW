#!/bin/bash
#SBATCH -t 00-24:00:00
#SBATCH -N 1
#SBATCH -o stdout_trim_nano_%j
#SBATCH -e stderr_trim_nano_%j
#SBATCH --job-name qc-map-nano

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if not submitting to SLURM

####################################################################

Combined_Universal_amplicon_nanopore.sh
created by Rose Kantor and Minxi Jiang

####################################################################
#
# This script performs the overall analysis of sequencing data from universal-amplicon sequencing followed by nanopore.
# It starts with  QC trimming, bowtie2 reads mapping and filtering,and post-mapping analysis.
#
#   1. QC trimming of the raw data
#  
#    1.1 Read cleaning with fastp 
#        Trims adaptors, polyX, polyG tails, minimum length 70, and base quality 15 (default)
#    1.2 [SKIP THIS STEP] Human reads removal
#    1.3 minimap2
#    1.4 Statistical analysis of the raw, cleaned, and deduplicated reads using seqkit stats
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
##########################################################################
STEP 1: QC trim

export PATH=/p/lustre1/preview/seq_software/samtools-1.17/:/p/lustre1/preview/seq_software/seq_software/minimap2-2.28_x64-linux/:/p/lustre1/preview/seq_software/:/p/lustre1/preview/seq_software/bbmap/:/p/lustre1/preview/seq_software/bedtools2/bin/:$PATH

INPUT_DIR=/p/vast1/preview/20241115_UCB_IAV_hac
OUT_DIR=/p/lustre1/preview/IAV_seq_data/nanopore_full

# run fastp
for f in $INPUT_DIR/*barcode*.fastq; do
	output_f=$OUT_DIR/step_1_qc/$(basename $f .fastq)_filt.fastq.gz
	fastp -l 500 -f 23 -t 23 --unqualified_percent_limit 10 --thread 100 \
		-g -x -y -q 10 --poly_g_min_len=20 --poly_x_min_len=20 \
		--n_base_limit 50 \
		-i $f -o $output_f --adapter_fasta $OUT_DIR/universal_primers.fasta \
		-h "$OUT_DIR/step_1_qc/$(basename $f .fastq)_fastp.html"
done

# minimap2
ref=/p/lustre1/preview/influenza/Spike_virus_genomes/same_length/remove_gaps/cut_ref_S1S2S3.fasta
for f in $OUT_DIR/step_1_qc/*_filt.fastq.gz; do
	samp=$(basename $f _filt.fastq.gz)
	out=$OUT_DIR/step_2_map/$samp
	minimap2 -t 100 -ax map-ont $ref $f > $out.sam
	samtools view -b $out.sam | samtools sort > $out.bam
	samtools index $out.bam
	rm $out.sam
done

# get stats
# seqkit stats -T -j 100 $INPUT_DIR/*barcode*.fastq -o $OUT_DIR/stats/stats_raw.txt
seqkit stats -T -j 100 $OUT_DIR/step_1_qc/*filt.fastq.gz -o $OUT_DIR/stats/stats_fastp.txt
#####################################################################################
STEP 2: bowtie2 mapping

# Define input and output paths
ref_genomes="/p/lustre1/preview/influenza/Spike_virus_genomes/same_length/remove_gaps/cut_ref_S1S2S3.fasta"
input="/p/lustre1/preview/IAV_seq_data/nanopore_full/step_1_fastp"
output="/p/lustre1/preview/IAV_seq_data/nanopore_full/STEP_2_bowtie"

# Create Bowtie2 index for the reference genome
bowtie2-build $ref_genomes ${input}/ref_genomes
samtools faidx $ref_genomes

# Loop over paired-end reads
for f in ${input}/*R1_001_filt.fastq.gz; do
  # Find the corresponding reverse read
  r=$(echo $f | sed 's/R1/R2/')

  # Extract the base name for the reads
  read_base=$(basename $f _R1_001_filt.fastq.gz)

  # Align paired-end reads to the reference genome
  bowtie2 -p 100 -x ${input}/ref_genomes --no-unal --no-discordant --no-mixed -1 $f -2 $r -S ${output}/${read_base}.sam

  # Convert SAM to BAM, sort, and index BAM file
  samtools view -bS ${output}/${read_base}.sam | samtools sort -o ${output}/${read_base}.sorted.bam
  samtools index ${output}/${read_base}.sorted.bam

  # Clean up the SAM file to save space
  rm ${output}/${read_base}.sam

  # Filter mapped paired-end reads with fewer than 5 edits
  reformat.sh in=${output}/${read_base}.sorted.bam out=${output}/${read_base}.filtered.bam editfilter=5

 # Filter, discard unpaired reads, and sort the BAM file
  samtools view -h -f 1 -b ${output}/${read_base}.filtered.bam | \
  samtools sort -o ${output}/${read_base}.filtered.paired.sorted.bam

  # Clean up intermediate BAM file
  rm ${output}/${read_base}.filtered.bam
done
######################################################################################
STEP 3: Post-mapping analysis
export PATH=/p/lustre1/preview/seq_software/bin/:/p/lustre1/preview/seq_software/bcftools/:/p/lustre1/preview/seq_software/:/p/lustre1/preview/seq_software/bowtie2-2.5.1/:/p/lustre1/preview/seq_software/samtools-1.17/:$PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/p/lustre1/preview/seq_software/htslib-1.18

# Set up reference and output directories
input_dir="/p/lustre1/preview/IAV_seq_data/nanopore_full/step_2_map"
output_dir="/p/lustre1/preview/IAV_seq_data/nanopore_full/step_3_depth_breadth/"
ref="/p/lustre1/preview/influenza/Spike_virus_genomes/same_length/remove_gaps/cut_ref_S1S2S3.fasta"
mkdir -p "$output_dir"

# Index the reference genome
# samtools faidx "$ref"

# Loop through BAM files in the input directory
for sorted_bam_file in $input_dir/*.bam; do

    # STEP 1: Index the sorted BAM file - skip, already done
    # samtools index "$sorted_bam_file" || { echo "Error: Failed to index BAM file $sorted_bam_file"; exit 1; }

    # Extract the basename of the BAM file for output file naming
    basename=$(basename "$sorted_bam_file" .bam)

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

