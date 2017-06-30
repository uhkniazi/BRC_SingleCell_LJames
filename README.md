# BRC_SingleCell_LJames
Analysis of Single Cell sequencing data for L James  

# Description of Scripts.  
### header.R and utilities.R  
Scripts with helper functions and global variables, e.g. to calculate lpd, lppd, predictions on fitted models etc.  

### 01_createDbEntry.R  
Create database entries for samples, files and annotations.  
##Raw data QC and Processing steps  
### 02_fastQC.R to 07_bam_files_qa.R  
Mostly data QC and processing steps. Look at quality plots and diagnostic plots to identify samples with possible issues.  
### 07.2_bam_files_qa.R  
Modified version of the original bam_files_qa script, to index original bam files before trimming to count the number of reads and perform QA.  
### 08_counts_from_bams.R  
Generate count matrix, while accounting for spike-in ERCC transcripts. Additionally it also uses the count matrix to select transcripts with greater than one coverage, extracts those transcript positions, checking if the width of the transcripts is large enough, chooses a sample of those transcripts, calculates coverage over the transcripts - using the function **getTranscriptCoverage**, bins the coverage into K bins and returns a vector. Plots this coverage to show 5 prime to 3 prime coverage.  
### 09_clustering_counts_matrix.R  
Uses the scater library and the single cell data QC workflow to drop low quality samples and generate QC plots along with PCA plots for samples based on read count matrix. We also choose the normalization strategy i.e. gene based or spike-in based, by examining the quality plots. The transcript coverage data generated in previous script is loaded, and the remaining samples passing QC are used to subset this matrix and make plots. The GC content vs average coverage over a transcript are also calculated and a plot is made.  
### 10_BASIC_array_job.R  
Script to run the BASIC array job on the FASTQ files to assemble the B Cell Receptors.  
### 11_AlignAndAssign.R  
The assembled sequences are aligned to the references sequences to assign a light and heavy chain to the sequences from respective samples. 
### 11.2_AlignAndAssign_imtg.R  
Very similar to the previous script but uses the IMTG database of sequences instead of the UCSC sequences to align and assign classes. The script downloads the references sequences fasta file from the IMTG database, extracts the sequences of interest using some grep and gsub string parsing. The rest is very similar to previous script.  
### 12_SummarizeResults.R  
Creates a CSV summary report with collated information from various steps performed earlier.  
### 13_createDbEntry.R  
Creates the database entries for data from NanoString.  
### 14_normalizeNanoString.R  
Uses the nanostringnorm package to import the nanoString data to normalize it and perform QC. Eventually normalised the data using top expressed genes after subtracting background. 
### 15_clusterCountMatrixNanoString.R  
Diagnostics for the normalized data along with clustering and plots using CDiagnostics plots.  
### 16_signatureGenesCellTypes.R  
This is a pretty long script that performs variable selection on a subset of genes common with single cell matrix. We perform a binomial regression using a optimization and SIR sampling based approach. Furthermore, the fun part here is calculating the AIC, pWAIC and WAIC on the data as a strategy to model selection. 
### 17_nanoStringDE.R  
Using a different approach, as we wanted to reduce number of groups. A random effects linear model is fit to each gene in the data, and using the top genes (with lowest adjust p-values) the data is clustered and eventually 3 and 2 group strategy is used.  
### 17.2_nanoStringDEwith3Groups.R  
This is a basic approach, and using the 3 new groups (instead of original 6) we new random effects model is fit, to get the top genes. These are then used to cluster and nanoString and Single cell data to see if 3 clusters can be separated using the single cell data.  
### 18_classifyNanoStToSingleCell.R  
This script is similar to the previous variable selection script (16), but this time we have a 2 class problem. This problem is mainly tackled using the CCrossValidation library, and following the steps of Random Forest, Exhaustive subset selection and k-fold Cross validation to select a set of 3 genes. These three genes are used then to cluster the nanoString and single cell data.  The second section of the script uses the binomial regression to do the same thing but is redundant and can be ignored.  













