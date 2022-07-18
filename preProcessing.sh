# mRNA Seqencing Pre-Processing Pipelines
# Ronghao Zhang on 2022-07-17


# ---------- Instruction of Script ------------------------------
# This script aims to pre-process the Nex-Generation-Sequencing (mRNA) Data.
# Before using this scripts, there are several pre-req need to be done. 
#   1. Download the following packages & tools. It is OK to download the following or newer version. 
#       a. if any package need to be download manually, put them into '~/Desktop/mRNA_rz_2022/seq_tools'
#   2. Download the referencing genome from GenBank or RefSeq of NCBI, rename as 'ref_genome.fna' -> '~/Desktop/mRNA_rz_2022/ref_genome'


# ---------- Basic Information ------------------------------
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


# ---------- Set-Up Directories ------------------------------
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


# ---------- FASTQC for Untrimmed Data ------------------------------
cd ~/Desktop/mRNA_rz_2022/mRNA_data_raw 
fastqc ./*.fastq.gz #perform quality check on all mRNA samples (both forward and reverse)
mv ./*fastqc.html ~/Desktop/mRNA_rz_2022/RNA_data_processed/fastqc_untirm
mv ./*fastqc.zip ~/Desktop/mRNA_rz_2022/RNA_data_processed/fastqc_untirm 

cd ~/Desktop/mRNA_rz_2022/RNA_data_processed/fastqc_untirm
for filename in *.zip
do
  unzip $filename
done

cat */summary.txt > ~/Desktop/mRNA_rz_2022/RNA_data_processed/docs/qc_summary/fastqc_summary_untrim.txt
grep FAIL ~/Desktop/mRNA_rz_2022/RNA_data_processed/docs/qc_summary/fastqc_summary_untrim.txt > \ 
          ~/Desktop/mRNA_rz_2022/RNA_data_processed/docs/qc_summary/fail_fastqc_summary_untrim.txt


# ---------- Filter and Trim mRNA Data ------------------------------
cd ~/Desktop/mRNA_rz_2022/mRNA_data_raw
for infile in *_R1.fastq.gz
do
  base=$(basename ${infile} _R1.fastq.gz)
  java -jar ~/Desktop/seq_tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
  ${infile} ${base}_R2.fastq.gz \
  ${base}_R1_trim.fastq.gz ${base}_R1_untrim.fastq.gz \
  ${base}_R2_trim.fastq.gz ${base}_R2_untrim.fastq.gz \
  ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:25
done
mv *trim.fastq.gz ~/Desktop/mRNA_rz_2022/mRNA_data_processed/mRNA_data_trim


# ---------- FASRQC for Trimmed Data ------------------------------
cd ~/Desktop/mRNA_rz_2022/mRNA_data_processed/mRNA_data_trim
fastqc ./*.fastq.gz 
mv ./*fastqc.html ~/Desktop/mRNA_rz_2022/RNA_data_processed/fastqc_tirm
mv ./*fastqc.zip ~/Desktop/mRNA_rz_2022/RNA_data_processed/fastqc_tirm 

cd ~/Desktop/mRNA_rz_2022/RNA_data_processed/fastqc_tirm
for filename in *.zip
do
  unzip $filename
done

cat */summary.txt > ~/Desktop/mRNA_rz_2022/RNA_data_processed/docs/qc_summary/fastqc_summary_trim.txt
grep FAIL ~/Desktop/mRNA_rz_2022/RNA_data_processed/docs/qc_summary/fastqc_summary_trim.txt > \ 
          ~/Desktop/mRNA_rz_2022/RNA_data_processed/docs/qc_summary/fail_fastqc_summary_trim.txt

# ---------- Sequence Alignment w/z ref_genome ------------------------------

# ---------- Compress into BAM & Sorting ------------------------------

# ---------- Variant Calling Format ------------------------------
