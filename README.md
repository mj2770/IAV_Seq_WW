# One Health Genomic Surveillance Tools for Influenza A Virus in Wastewater
This repository contains the Probe/Primer design pipelines, sequencing data analysis pipelines, and other documentation associated with the study:
“Developing and Benchmarking One Health Genomic Surveillance Tools for Influenza A Virus in Wastewater” [(Jiang et al., 2025)](https://doi.org/10.1101/2025.09.19.676942)
## 👥 Contributions
- **Rose Kantor, Minxi Jiang** – pipeline design and processing 
- **Co-contributor ** – visualization  

## 📄 Project Overview

Influenza A viruses (IAV) remain a persistent One Health threat, and whole-genome sequencing from wastewater offers a promising surveillance tool. However, IAV is at low abundance in wastewater, making it difficult to sequence. 
We benchmarked four targeted enrichment methods suited for whole-genome sequencing including:
- HA Tiled-amplicon (custom, subtype-specific)
- Probe-IAV (custom probe-capture panel)
- Probe-Twist (commercial Twist CVRP panel)
- Universal-amplicon (whole-genome PCR method Zhou et al. 2009)

We evaluated their performance across spike-in mixtures, wastewater extracts, and concentration/extraction methods.

## 🔬 Pipeline Overview

### 1. Probe & Primer Design
- **Probe panels**  
  - Custom IAV probes designed using [Syotti](https://github.com/jnalanko/syotti) and in-house scripts  
  - Commercial probe panel: Twist Comprehensive Viral Research Panel (CVRP)  
- **Primers**  
  - HA-segment tiled amplicons designed with [PriMux](https://pubmed.ncbi.nlm.nih.gov/25157264/)  
  - Universal whole-genome amplicon primers adapted from published protocols  
- Outputs: BED/FASTA files of probe/primer sets used in targeted sequencing assays  

### 2. Preprocessing
- Input: Illumina/Nanopore raw FASTQ files  
- Tools:  
  - [fastp](https://github.com/OpenGene/fastp): adapter trimming, polyX removal, quality filtering  
  - [cutadapt](https://cutadapt.readthedocs.io/): primer trimming (tiled-amplicon workflows)  
  - [seqkit](https://bioinf.shenwei.me/seqkit/): QC and summary statistics  

### 3. Read Mapping
- Aligners:  
  - [Bowtie2](https://github.com/BenLangmead/bowtie2) for Illumina  
  - ONT basecalled reads for Nanopore  
- Post-processing:  
  - Reformat ([BBMap](https://jgi.doe.gov/data-and-tools/bbtools/)): remove reads with >5 mismatches/indels  
  - [Samtools](http://www.htslib.org/): mapping stats, coverage depth/breadth  

### 4. Variant Calling & Consensus
- [Samtools mpileup](http://www.htslib.org/doc/samtools.html)  
- [iVar](https://github.com/andersen-lab/ivar): variant calling & consensus generation  
- [BCFtools](http://samtools.github.io/bcftools/): VCF processing  

### 5. Coverage & Sensitivity Analysis
- Compute: coverage breadth (%) and coverage depth (RPKM)  
- Normalization equations implemented in `analysis/coverage_metrics.py`  
- Visualization with custom R/Python plotting scripts  

### 6. Comparative Analyses
- Coverage comparison across spike-in levels (S1–S3)  
- Sensitivity: correlation of input copies (dPCR) vs sequencing depth  
- Segment-level bias: uniformity across IAV genome segments  
- Decay/extraction effects: IP vs PMG concentration methods  

### 7. Economic Analysis
- Scripts in `analysis/cost_analysis.py` summarize reagent, sequencing, and labor costs across methods.  
