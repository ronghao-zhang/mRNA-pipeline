# mRNA-pipeline

This repository stores `shell` and `R` programming languages pipelines that are used to process mRNA data. 
The mRNA fragments are collected from rata and are sequenced using Illumina Nest Generation Sequencing Platforms. 

## Pre-Processing Pipeline
### Editing File Directories
The Pre-Processing Pipeline is written in `shell` scripts. It can be used directly in Ubuntu. 
If you want to change directory of your raw sequencing data, please follow the following steps in a terminal:
1. Install vim using the command `sudo apt install vim` .
2. View the script using the command `vim YOUR_DIREC_TO_FILE/preProcessing.sh`.
3. Allow Eidting the script by hit the button `i` on your keyboard. 
4. Change the directory under section *FASTQC for Untrimmed Data* and *Filter and Trim mRNA Data*. They are currently set to `~/Desktop/mRNA_rz_2022/mRNA_data_raw`. 
5. Save and Exit by hitting `Esc` and typing `:wq` and hitting `Enter` or `Return`. 
### Run the Script
To run `preProcessing.sh`, open a terminal and type `bash YOUR_DIREC_TO_FILE/preProcessing.sh`. 
