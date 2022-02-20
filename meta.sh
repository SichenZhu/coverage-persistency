#!/bin/bash

#BATCH --job-name="metagenomics"
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=aquila,parallel,hetao
#SBATCH --array=0-5
#SBATCH --output=meta.%A_%a.out

data_dir="raw_pair_ended_data" 
for dir in $( ls ${data_dir} )
do
    dir_name+=(${dir})
done

current_dir_name=${dir_name[SLURM_ARRAY_TASK_ID]}


runDir="${data_dir}/${current_dir_name}/${current_dir_name}-meta" 
rawdataDir="${data_dir}/${current_dir_name}/raw_data" # !!need to change SPAdes & Bowtie2 input command!!
codeDir="code" 
dbDir="database"

# -------------------------------------------- Step1:Assemble raw pair ended raw reads (SPAdes) --------------------------------------------
module purge
module load spades/3.14.1

mkdir ${runDir}/SPAdes

now=$(date +"%T")
echo "SPAdes Start (Step1): $now"

spades.py -1 ${rawdataDir}/${current_dir_name}_1.fastq.gz -2 ${rawdataDir}/${current_dir_name}_2.fastq.gz -o ${runDir}/SPAdes -t 20 --meta

now=$(date +"%T")
echo "SPAdes Finish (Step1): $now"


# -------------------------------------------- Step2:Diamond process -------------------------------------------------------------------------
module purge
module load diamond/0.9.36

mkdir ${runDir}/DIAMOND

# ----------------------------------- Step2.1: Diamond build database ----------------------------------------------
find path_to_database_eggNOG/bacteria/2_raw/ | xargs zcat | diamond makedb -d 2_raw_algs_db -p 40
echo "Step2.1: DIAMOND database built!"

# ----------------------------------- Step2.2: Diamond find alignment ----------------------------------------------
now=$(date +"%T")
echo "DIAMOND alignment Start (Step2.2): $now"

diamond blastx -d ${dbDir}/2_raw_algs_db.dmnd -q ${runDir}/SPAdes/scaffolds.fasta -o ${runDir}/DIAMOND/DIAMOND_raw_output.m8 -p 20
# -e 0.01 

now=$(date +"%T")
echo "DIAMOND alignment Finish (Step2.2): $now"

module purge
module load python/intel/3.6.6

python ${codeDir}/DIAMOND_output_find_OG_PI_interval.py   ${dbDir}/COGFC_GN.tsv.gz  ${runDir}/DIAMOND/DIAMOND_raw_output.m8   ${runDir}/DIAMOND/qseq_sseq_GNs_num.tsv

# ----------------------------------- Step2.3: Reformat DIAMOND output (add OG; PI; Node_idx; Contig_length)---------

awk '{print $2"\t"$3}' ${runDir}/DIAMOND/qseq_sseq_GNs_num.tsv | sort | uniq -c > ${runDir}/DIAMOND/D_temp1.tsv  

awk '{print $2"\t"$3}' ${runDir}/DIAMOND/D_temp1.tsv > ${runDir}/DIAMOND/sorted_nondup_COGFunC_OG.tsv   

awk '{print $1}' ${runDir}/DIAMOND/D_temp1.tsv > ${runDir}/DIAMOND/add_dnum.tsv

python ${codeDir}/find_highest_OG_PI.py  ${runDir}/DIAMOND/sorted_nondup_COGFunC_OG.tsv  ${dbDir}/COG_sort_PI.tsv  ${dbDir}/nCOG_PI.tsv  ${runDir}/DIAMOND/nondup_COGFunC_OG_PI.tsv

awk '{print $1"\t"$2"\t"$3}' ${runDir}/DIAMOND/nondup_COGFunC_OG_PI.tsv > ${runDir}/DIAMOND/D_temp2.tsv
paste -d $'\t' ${runDir}/DIAMOND/nondup_COGFunC_OG_PI.tsv ${runDir}/DIAMOND/add_dnum.tsv > ${runDir}/DIAMOND/nondup_COGFunC_OG_PI_dnum_unsorted.tsv


sort -t $'\t' -k1,1 ${runDir}/DIAMOND/nondup_COGFunC_OG_PI_dnum_unsorted.tsv > ${runDir}/DIAMOND/nondup_COGFunC_OG_PI_dnum.tsv
awk '{print NR"\t"$2}' ${runDir}/DIAMOND/DIAMOND_raw_output.m8 | sort -t $'\t' -k2,2 > ${runDir}/DIAMOND/sorted_COGFC_diamond.tsv


python ${codeDir}/duplicated_OG_PI.py  ${runDir}/DIAMOND/sorted_COGFC_diamond.tsv  ${runDir}/DIAMOND/nondup_COGFunC_OG_PI_dnum.tsv   ${runDir}/DIAMOND/idx_COGFunC_OG_PI.tsv


awk '{print NR"\t"$1"\t"$2"\t"$7"\t"$8"\t"$3}' ${runDir}/DIAMOND/DIAMOND_raw_output.m8 > ${runDir}/DIAMOND/D_temp3.tsv
sort -t $'\t' -k3,3 ${runDir}/DIAMOND/D_temp3.tsv > ${runDir}/DIAMOND/sorted_idx_diamond_output.tsv

sort -k1 -n ${runDir}/DIAMOND/idx_COGFunC_OG_PI.tsv > ${runDir}/DIAMOND/D_temp4.tsv
sort -k1 -n ${runDir}/DIAMOND/sorted_idx_diamond_output.tsv > ${runDir}/DIAMOND/D_temp5.tsv
paste -d $'\t' ${runDir}/DIAMOND/D_temp5.tsv ${runDir}/DIAMOND/D_temp4.tsv > ${runDir}/DIAMOND/D_temp6.tsv
awk '{printf "%d\t%s\t%s\t%s\t%8.5f\t%d\t%d\t%4.1f\n", $1,$2,$3,$9,$10,$4,$5,$6}' ${runDir}/DIAMOND/D_temp6.tsv > ${runDir}/DIAMOND/D_temp7.tsv
awk "{print $2}" ${runDir}/DIAMOND/D_temp7.tsv | awk -F "_" '{print $4}' > ${runDir}/DIAMOND/len.tsv
paste -d $'\t' ${runDir}/DIAMOND/D_temp7.tsv ${runDir}/DIAMOND/len.tsv > ${runDir}/DIAMOND/final_format_diamond_output.tsv

sort -t $'\t' -k9nr -k6nr -k1n ${runDir}/DIAMOND/final_format_diamond_output.tsv | awk '{printf "%d\t%d\t%s\t%8.5f\t%d\t%d\t%4.1f\t%s\t%s\n", $1, $9, $4, $5, $6, $7, $8, $2, $3}' > ${runDir}/DIAMOND/reformat_diamond_output.tsv

# --------------------------------- Step2.4: Change start and end point in the reformat DIAMOND output ------------
python ${codeDir}/change_s_e.py   ${runDir}/DIAMOND/reformat_diamond_output.tsv ${runDir}/DIAMOND/change_s_e.tsv

awk "{print $7}" ${runDir}/DIAMOND/change_s_e.tsv | awk -F "_" '{print $2}' > ${runDir}/DIAMOND/node.tsv 
paste -d $'\t' ${runDir}/DIAMOND/change_s_e.tsv ${runDir}/DIAMOND/node.tsv > ${runDir}/DIAMOND/change_s_e_node_unsorted.tsv

sort -t $'\t' -k2n -k10n -k5n -k6n ${runDir}/DIAMOND/change_s_e_node_unsorted.tsv | awk '{printf "%d\t%d\t%s\t%8.5f\t%d\t%d\t%4.1f\t%s\t%s\t%d\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > ${runDir}/DIAMOND/change_s_e_node.tsv

# -------------------------------- Step2.5: Find the longest segment for each OG ----------------------------------
python ${codeDir}/find_overlap.py  ${runDir}/DIAMOND/change_s_e_node.tsv ${runDir}/DIAMOND/longest_per_og.tsv


awk '{print $1}' ${runDir}/DIAMOND/longest_per_og.tsv | awk -F "_" '{printf  "%d\t%d\n", $4, $2}' > ${runDir}/DIAMOND/len_node.tsv 
paste -d $'\t' ${runDir}/DIAMOND/longest_per_og.tsv ${runDir}/DIAMOND/len_node.tsv > ${runDir}/DIAMOND/longest_per_og_len_node.tsv

sort -k6n -k7n -k4n -k5n ${runDir}/DIAMOND/longest_per_og_len_node.tsv > ${runDir}/DIAMOND/sorted_longest_per_og_len_node.tsv


# -------------------------------- Step2.6: find overlap & nonoverlap part between each OG -------------------------
python ${codeDir}/find_overlap2.py  ${runDir}/DIAMOND/sorted_longest_per_og_len_node.tsv ${runDir}/DIAMOND/overlap.tsv ${runDir}/DIAMOND/nonoverlap.tsv

awk '{printf "%d\t%d\t%d\t%s\t%8.5f\n", $7,$4,$5,$2,$3}' ${runDir}/DIAMOND/overlap.tsv | sort -k1n -k2n -k3n | awk '{printf "%d\t%d\t%d\t%s\t%8.5f\n", $1,$2,$3,$4,$5}' > ${runDir}/DIAMOND/overlap_node_s_e_og_pi.tsv

# -------------------------------- Step2.7: choose PI for overlapping segment --------------------------------------
python ${codeDir}/number_of_overlap.py ${runDir}/DIAMOND/overlap_node_s_e_og_pi.tsv  ${runDir}/DIAMOND/record_overlap_interval_len.tsv  ${runDir}/DIAMOND/record_nonverlap_interval_len.tsv

now=$(date +"%T")
echo "DIAMOND temporarily Finish (Step5): $now"

# -------------------------------- Step3.1: Bowtie2 --------------------------------------------------------------------------------------------------------
module purge
module load bowtie2/intel/2.3.4.1

mkdir ${runDir}/Bowtie2

now=$(date +"%T")
echo "Bowtie2 alignment Start (Step3.1): $now"

time bowtie2-build ${runDir}/SPAdes/scaffolds.fasta ${runDir}/Bowtie2/data
time bowtie2 -x ${runDir}/Bowtie2/data \
-1 ${rawdataDir}/${current_dir_name}_1.fastq.gz \
-2 ${rawdataDir}/${current_dir_name}_2.fastq.gz \
-t -p 20 \
--no-unal \
-S ${runDir}/Bowtie2/bowtie2_output.sam
# !!!Attention: -U is only for this specific case

now=$(date +"%T")
echo "Bowtie2 Finish (Step3.1): $now"

# -------------------------------- Step3.2: Bowtie2 output conversion ---------------------------------------------------------------
module purge
module load samtools/intel/1.9
module load bedtools2/2.27.1

now=$(date +"%T")
echo "Bowtie2 output conversion Start (Step3.2): $now"

samtools sort -o ${runDir}/Bowtie2/bowtie2_output_sorted.bam -@ 20 ${runDir}/Bowtie2/bowtie2_output.sam

now=$(date +"%T")
echo "Bowtie2 output conversion Finish (Step3.2): $now"

# -------------------------------- Step4: Calculate coverage ----------------------------------------------------------------------------------------------
mkdir ${runDir}/coverage

now=$(date +"%T")
echo "Coverage calculation Start (Step4): $now"

sed -n '/@SQ/p' ${runDir}/Bowtie2/bowtie2_output.sam | awk '{print $2}' | awk -F ":" '{print $2}' > ${runDir}/coverage/seq.tsv
awk -F "_" '{print $4}' ${runDir}/coverage/seq.tsv > ${runDir}/coverage/len.tsv
paste -d $'\t' ${runDir}/coverage/seq.tsv ${runDir}/coverage/len.tsv > ${runDir}/coverage/Bowtie2_output.genome

awk '{print $1}' ${runDir}/coverage/Bowtie2_output.genome | awk -F '_' '{printf "%d\n", $2}' > ${runDir}/coverage/node.tsv
paste -d $'\t' ${runDir}/coverage/node.tsv ${runDir}/coverage/len.tsv > ${runDir}/coverage/uniq_node_len.tsv  

bamToBed -i ${runDir}/Bowtie2/bowtie2_output_sorted.bam > ${runDir}/coverage/bowtie2_output.bed

genomeCoverageBed -i ${runDir}/coverage/bowtie2_output.bed -g ${runDir}/coverage/Bowtie2_output.genome -d > ${runDir}/coverage/positionalwise_coverage_from_bedtools.tsv # original name: cov_cal.tsv

awk -F "_" '{printf "%d\t%d\n", $2,$4}' ${runDir}/coverage/positionalwise_coverage_from_bedtools.tsv > ${runDir}/coverage/node_len.tsv

awk '{printf "%d\t%d\n",$2,$3}' ${runDir}/coverage/positionalwise_coverage_from_bedtools.tsv > ${runDir}/coverage/pos_cov.tsv

paste -d $'\t' ${runDir}/coverage/node_len.tsv ${runDir}/coverage/pos_cov.tsv > ${runDir}/coverage/node_len_pos_cov.tsv


now=$(date +"%T")
echo "Coverage calculation Finish (Step4): $now"

# -------------------------------- Step5: Expand DIAMOND & Bowtie2 output to each position -----------------------------------------------------------------
module purge
module load python/intel/3.6.6

mkdir ${runDir}/expand

now=$(date +"%T")
echo "Expanding to each position Start (Step5): $now"

awk '{printf "%d\t%d\t%d\t%s\t%8.5f\n",$7,$4,$5,$2,$3}' ${runDir}/DIAMOND/nonoverlap.tsv | sort -t $'\t' -k1,1n -k2,2n -s > ${runDir}/expand/sorted_no_overlap.tsv
sed '$d' ${runDir}/DIAMOND/record_overlap_interval_len.tsv | awk '{printf "%d\t%d\t%d\t%s\t%8.5f\n",$1,$2,$3,$5,$6}' > ${runDir}/expand/overlap_choose_higher_PI.tsv
sed '$d' ${runDir}/DIAMOND/record_nonverlap_interval_len.tsv | awk '{printf "%d\t%d\t%d\t%s\t%8.5f\n",$1,$2,$3,$5,$6}' > ${runDir}/expand/nonoverlap_in_overlap.tsv
cat ${runDir}/expand/sorted_no_overlap.tsv ${runDir}/expand/overlap_choose_higher_PI.tsv ${runDir}/expand/nonoverlap_in_overlap.tsv > ${runDir}/expand/merged_3_files.tsv

sort -t $'\t' -k1,1n -k2,2n -s ${runDir}/expand/merged_3_files.tsv > ${runDir}/expand/interval_result.tsv

python ${codeDir}/expand_to_each_pos.py ${runDir}/coverage/uniq_node_len.tsv ${runDir}/expand/interval_result.tsv ${runDir}/expand/pos_PI.tsv ${runDir}/expand/cov_no_pi.tsv ${runDir}/expand/pi_no_cov.tsv
#!!!REMARK!!!: need to check pi_no_cov.tsv: whether it is empty

paste -d $'\t' ${runDir}/coverage/positionalwise_coverage_from_bedtools.tsv ${runDir}/expand/pos_PI.tsv > ${runDir}/expand/node_pos_cov_posinfo_PI_og.tsv 

now=$(date +"%T")
echo "Expanding to each position Finish (Step5): $now"

# -------------------------------- Step6: TPM calculation --------------------------------------------------------------------------------------------------
module purge
module load anaconda/2.3.0

mkdir ${runDir}/TPM
mkdir ${runDir}/TPM/TSV

now=$(date +"%T")
echo "TPM Start (Step6): $now"

awk -F "_" '$4>100' ${runDir}/expand/node_pos_cov_posinfo_PI_og.tsv > ${runDir}/expand/node_pos_cov_info_PI_og_100.tsv

python ${codeDir}/PI_on_HPC.py  ${runDir}/expand/node_pos_cov_info_PI_og_100.tsv  ${runDir}/SPAdes/scaffolds.fasta ${runDir}/TPM/TSV ${runDir}/TPM/

now=$(date +"%T")
echo "TPM Finish (Step6) & Concatenation (Step7) Start: $now"

# -------------------------------- Step7: final step: concatenate the whole hashtable -------------------------------------------------------------------------
touch ${runDir}/TPM/all_contigs_node_pos_og_pi_tpm_nt.tsv

for f in ${runDir}/TPM/TSV/*
do
cat "$f" >> ${runDir}/TPM/all_contigs_node_pos_og_pi_tpm_nt.tsv
done

now=$(date +"%T")
echo "TPM Concatenation Finish (Step7) & Plotting (Step8) Start: $now"


# -------------------------------- Step8: Plotting TPM vs PI figure -------------------------------------------------------------------------------------------
module purge
module load python/intel/3.6.6

python ${codeDir}/TPM_PI_plot.py  ${runDir}/TPM/serial.out  ${runDir}/TPM/TPM_PI.pdf

now=$(date +"%T")
echo "Plotting (Step8) Finish: $now"
echo "ALL DONE!"


