#!/bin/bash
#$ -P myers.prjc -q long.qc -j yes -V
##this code is used to process the bam files, which include merging, duplication marking and sort the reads 

dir_picard=/well/myers/ranli/software/picard-tools-1/picard-tools-1.115
bam_dir=/well/myers/ranli/data/MouseSequencing_bam/F5/original_bam_F5_10_Jun_2016
bam_dir2=/well/myers/ranli/data/MouseSequencing_bam/F5/bamfiles
dir_ref=/well/myers/ranli/data/mm10

SGE_TASK_ID=$1
java -Xmx8G -jar ${dir_picard}/MergeSamFiles.jar INPUT=${bam_dir}/WTCHG_274207_${SGE_TASK_ID}.bam  INPUT=${bam_dir}/WTCHG_277894_${SGE_TASK_ID}.bam INPUT=${bam_dir}/WTCHG_277895_${SGE_TASK_ID}.bam INPUT=${bam_dir}/WTCHG_277896_${SGE_TASK_ID}.bam INPUT=${bam_dir}/WTCHG_277897_${SGE_TASK_ID}.bam INPUT=${bam_dir}/WTCHG_277898_${SGE_TASK_ID}.bam OUTPUT=${bam_dir2}/WTCHG_${SGE_TASK_ID}.bam SORT_ORDER=coordinate

temp=/well/myers/ranli/temp
dir_bam=${bam_dir2}
java -Xmx8g -jar ${dir_picard}/SortSam.jar INPUT=${dir_bam}/WTCHG_${SGE_TASK_ID}.bam OUTPUT=${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted.bam SORT_ORDER=coordinate TMP_DIR=${temp}


######################mark duplicate ######################################################
java -Xmx8g -jar ${dir_picard}/MarkDuplicates.jar INPUT=${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted.bam OUTPUT=${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted_dedup.bam METRICS_FILE=${dir_bam}/WTCHG_${SGE_TASK_ID}.matrics TMP_DIR=${temp}

############# recorder sam ####################################

java -Xmx8g -jar ${dir_picard}/ReorderSam.jar INPUT=${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted_dedup.bam OUTPUT=${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted_dedup_OrderByRef.bam  R=${dir_ref}/mm10.fa CREATE_INDEX=true TMP_DIR=${temp}


