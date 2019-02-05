#!/bin/bash
#PBS -l walltime=120:00:00
#PBS -o extract_all_mr_paper_snps.txt
#PBS -me
#PBS -l nodes=1:ppn=4
#
#---------------------------------------------
#
cd "${HOME}/"
#
#---------------------------------------------

export DIR1="*"
export DIR2="**"
export DIR3="***"


# The loop below is based on Kaitlins script and uses qctool to extract SNPs of interest. Data is written to a .gen file which plink can read. 

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23
do
qctool -g $DIR1/data.chr"$i".bgen -og $DIR2/"$i"_snp.dosage.gen -s $DIR3/data.chr01.sample -omit-chromosome -incl-rsids $DIR2/AAM389_snps.txt -excl-samples $DIR2/biobank_exclusion_list_unique_ids.txt
done

rm $DIR2/snps-out.gen
touch $DIR2/snps-out.gen

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23
do
cat $DIR2/${i}_snp.dosage.gen >> $DIR2/snps-out.gen
done
 
