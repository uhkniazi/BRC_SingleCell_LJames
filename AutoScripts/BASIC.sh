# Autogenerated script from BASIC_array_job.R 
# date Tue Apr  4 10:38:35 2017
# make sure directory paths exist before running script
#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 2
#$ -cwd
#$ -N BASIC-array
#$ -j y
#$ -l h_vmem=19G
#$ -t 1-71



module load bioinformatics/bowtie2/2.2.5
module load general/python/3.5.1



# Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=BASIC_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`
inr2=`sed -n ${number}p $paramfile | awk '{print $2}'`
outdir=`sed -n ${number}p $paramfile | awk '{print $3}'`

# create the output directory
mkdir $outdir

# 9. Run the program.
python3.5 BASIC.py -p 2 -b /opt/apps/bioinformatics/bowtie2/2.2.5/ -PE_1 $inr1 -PE_2 $inr2 -o $outdir



