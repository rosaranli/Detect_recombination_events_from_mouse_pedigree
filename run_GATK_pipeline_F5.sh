#!/bin/bash
#$ -P myers.prjc -q long.qc -j yes -V
##this code is used to call variants from bam files of mice samples
SGE_TASK_ID=$1


dir_GATK=/well/myers/ranli/software/GenomeAnalysisTK-3.3-0
key_file=${dir_GATK}/ranli_well.ox.ac.uk.key
dir_ref=/well/myers/ranli/data/mm10
dir_bam=/well/myers/ranli/data/MouseSequencing_bam/F5/bamfiles
dir_known=/well/myers/ranli/data/MouseSequencing_bam/reference_SNPs_Indels/MGP_v4/MGP_V4_chr

java -Xmx8g -jar ${dir_GATK}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${dir_ref}/mm10.fa -I ${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted_dedup_OrderByRef.bam -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr1.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr2.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr3.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr4.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr5.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr6.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr7.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr8.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr9.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr10.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr11.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr12.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr13.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr14.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr15.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr16.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr17.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr18.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr19.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chrX.vcf -o ${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted_dedup_OrderByRef.intervals

########################### local realign #################################
java -Xmx8g -jar ${dir_GATK}/GenomeAnalysisTK.jar -T IndelRealigner -R ${dir_ref}/mm10.fa -I ${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted_dedup_OrderByRef.bam -targetIntervals ${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted_dedup_OrderByRef.intervals -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr1.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr2.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr3.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr4.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr5.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr6.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr7.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr8.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr9.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr10.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr11.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr12.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr13.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr14.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr15.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr16.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr17.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr18.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chr19.vcf -known ${dir_known}/mgp.v4.indels.dbSNP_B6_CAST_chrX.vcf -o ${dir_bam}/WTCHG_${SGE_TASK_ID}_sorted_dedup_OrderByRef_realigned.bam

#################################################################################

dir_vcf=/well/myers/ranli/data/MouseSequencing_bam/F5/VCFfiles
dir_BQSR=/well/myers/ranli/data/MouseSequencing_bam/F5/bamfiles
dir_knownSites=/well/myers/ranli/data/MouseSequencing_bam/reference_SNPs_Indels/MGP_v4/MGP_V4_chr
middle_bamfile="sorted_dedup_OrderByRef_realigned"
#middle_bamfile "sorted_dedup_OrderByRef_without_KnownSites_realigned"
java -Xmx8g -jar ${dir_GATK}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 4 -R ${dir_ref}/mm10.fa -I ${dir_bam}/WTCHG_${SGE_TASK_ID}_${middle_bamfile}.bam  -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr1.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr2.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr3.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr4.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr5.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr6.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr7.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr8.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr9.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr10.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr11.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr12.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr13.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr14.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr15.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr16.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr17.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr18.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chr19.vcf -knownSites ${dir_knownSites}/mgp.v4.snps.dbSNP_B6_CAST_chrX.vcf -o ${dir_BQSR}/WTCHG_${SGE_TASK_ID}_${middle_bamfile}_recal_MGPv4Known.table -et NO_ET -K ${key_file}

######################print reads#####################

java -Xmx8g -jar ${dir_GATK}/GenomeAnalysisTK.jar -T PrintReads -nct 2 -R ${dir_ref}/mm10.fa -I ${dir_bam}/WTCHG_${SGE_TASK_ID}_${middle_bamfile}.bam  -BQSR ${dir_BQSR}/WTCHG_${SGE_TASK_ID}_${middle_bamfile}_recal_MGPv4Known.table -o ${dir_bam}/WTCHG_${SGE_TASK_ID}_${middle_bamfile}_recal_MGPv4Known.bam -et NO_ET -K ${key_file}

##############################calling variants ##################################################


