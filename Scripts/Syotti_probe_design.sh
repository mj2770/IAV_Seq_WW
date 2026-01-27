#!/bin/bash
#SBATCH -t 00-24:00:00
#SBATCH -N 1
#SBATCH -o stdout_syotti_v1_%j
#SBATCH -e stderr_syotti_v1_%j 
#SBATCH --job-name preview_seq
#SBATCH --mail-user mj2770@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -A careswwm

######################################################################
STEP 1: Probe design
export PATH=/p/lustre1/preview/seq_software/syotti/bin/:/p/lustre1/preview/seq_software/ncbi-blast-2.14.0+/bin/:/p/lustre1/preview/seq_software/vsearch-2.23.0/bin/:/p/lustre1/preview/seq_software/:/p/lustre1/preview/seq_software/bbmap/:$PATH

input_dir="/p/lustre1/preview/influenza/sequencing_assay_design/clustered100_fastas/others/"
output_dir="/p/lustre1/preview/influenza/probe_capture/Syotti/others"

for target in $input_dir/*_filtered.clustered100.fasta; do
    base_name=$(basename "$target" _filtered.clustered100.fasta)
    output_file="${output_dir}/${base_name}_probes.fa"
    syotti design --bait-len 120 --hamming-distance 40 -s "$target" -r -o "$output_file"
done

############################################################################
STEP 2: in-silico testing of the probes

export PATH=/p/lustre1/preview/seq_software/syotti/bin/:/p/lustre1/preview/seq_software/ncbi-blast-2.14.0+/bin/:/p/lustre1/preview/seq_software/vsearch-2.23.0/bin/:/p/lustre1/preview/seq_software/:/p/lustre1/preview/seq_software/bbmap/:$PATH

input="/p/lustre1/preview/influenza/Spike_virus_genomes/"
FM_output="/p/lustre1/preview/influenza/probe_capture/Syotti/Examining/"
bait="/p/lustre1/preview/influenza/probe_capture/Syotti/baits/"

for spike in $input/*.fasta; do
  base_name=$(basename "$spike" .fasta)
  output_file="${FM_output}/${base_name}.fmi"
  syotti index -s $spike -o "$output_file"
  syotti examine -d 40 -s $spike -f "$output_file" -b "${bait}/${base_name}_bait.fna" -o "${FM_output}/${base_name}.out"
done
