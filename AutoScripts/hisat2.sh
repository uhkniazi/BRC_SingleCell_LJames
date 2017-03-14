# Autogenerated script from write_hisat2_script.R
# date Tue Mar 14 13:26:23 2017
# make sure directory paths exist before running script
#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8
#$ -cwd
#$ -N hisat2-array
#$ -j y
#$ -l h_vmem=19G
#$ -t 1-71



module load bioinformatics/hisat2/2.0.4



# Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=hisat2_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`
inr2=`sed -n ${number}p $paramfile | awk '{print $2}'`
outsam=`sed -n ${number}p $paramfile | awk '{print $3}'`

# 9. Run the program.
hisat2 --trim5 10 -x /users/k1625253/brc_scratch/Data/MetaData/GenomeIndex/hg38_spikein_hisat2/hg38/genome -1 $inr1 -2 $inr2 -S $outsam



