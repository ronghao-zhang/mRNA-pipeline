# mRNA Seqencing Pre-Processing Pipelines
# Ronghao Zhang on 2022-07-17

# ---------- Basic Information ----------
# package & tools used
#   1. FastQC High Throughput Sequence QC Report - version 0.11.9
#   2. Trimmomatic - version 0.39
#   3. bwa - version 0.7.17-r1188
#   4. samtools - version 1.12
#   5. bcftools - version 1.13
# referencing genome information: RefSeq NCBI Retreived 2022-07-08.
#   https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/
#     latest_assembly_versions/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.fna.gz
#   https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/
#     latest_assembly_versions/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.gbff.gz

# ---------- Set-Up Directories ----------
# Desktop/mRNA_rz_2022 -> mRNA_data_raw (36 samples PE in .fastq format)
#                      -> ref_genome 
#                      -> seq_tools (contains samtools & trimmomatic)
#                      -> mRNA_data_processed -> fastqc_untirm
#                                             -> mRNA_data_trim
#                                             -> fastqc_trim
#                                             -> sam
#                                             -> bam
#                                             -> bcf
#                                             -> vc
#                                             -> docs -> qc_summary
cd ~/Desktop/mRNA_rz_2022/
mkdir mRNA_data_processed
cd ~/Desktop/mRNA_rz_2022/mRNA_data_processed
mkdir fastqc_untrim mRNA_data_trim fastqc_trim sam bam bcf vcf docs
mkdir ~/Desktop/mRNA_rz_2022/mRNA_data_processed/docs/qc_summary

# ---------- FASTQC for Untrimmed Data ----------

# ---------- Filter and Trim mRNA Data ----------

# ---------- FASRQC for Trimmed Data ----------

# ---------- Sequence Alignment w/z ref_genome ----------

# ---------- Compress into BAM & Sorting ----------

# ---------- Variant Calling Format ----------
