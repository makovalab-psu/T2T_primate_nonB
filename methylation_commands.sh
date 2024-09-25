################## METHYLATION IN SATELLITES AND OTHER REPEATS #################
# written by Linnéa Smeds, September 2024
# Code for analysing the methylation in satellites and repeats enriched for 
# G-quadruplexes. This analysis is run only for the primary haplotype

# Download methylation data
mkdir ref/methylation
cd ref/methylation
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/regulation/chm13v2.0_hg002_CpG_ont_guppy6.1.2.bed
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/regulation/chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw
# Convert bigwig to bed 
~/software/bigWigToWig chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.wig
# (bigWigToBed returns a 0-based format, so use special flag for wig2bed)
wig2bed --zero-indexed <chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.wig >chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bed
# remove id column4 to make the column order similar to the H002 file
 cut -f1,2,3,5 chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bed >chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.formatted.bed
cd ../.. 
# The two methylation files are listed in:
# T2T_primate_nonB/helpfiles/human_methylation.txt


# We will focus on the following repeats, that showed enrichment of G4s in 
# humans:
# LSAU, retro/SVA, HSAT5, ACRO1, GSAT, GSATII, GSATX, MSR1, SSTI, TAR1, rRNA,
# COMP-subunit_LSAU-BSAT_rnd-1_family-1, COMP-subunit_LSAU-BSAT_rnd-1_family-4,
# COMP-subunit_LSAU-BSAT_rnd-6_family-5403, SAT-VAR_rnd-6_family-3554
# They are listed in T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt

mkdir -p methylation/human/

# Save the interesting repeats as bed
python3 T2T_primate_nonB/python/extract_from_bed.py -b repeats/human/RepeatMasker.bed -l T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt -c 4 -o methylation/human/repeats_of_interest.bed

# Get average methylation for each of these repeat categories
module load bedtools/2.31.0
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo "Repeat Avg_meth CpG_sites" | sed 's/ /\t/g' >methylation/human/methylation_in_roi_$ind.txt
    for rep in $(cat T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt)
    do
    meth=`awk -v r=$rep '($4==r){print}' methylation/human/repeats_of_interest.bed |intersectBed -a $filename -b - |awk -v sum=0 '{sum+=$4}END{m=sum/NR; print m, NR}'`
    echo $rep" "$meth | sed 's/ /\t/g' >>methylation/human/methylation_in_roi_$ind.txt
    done
done 

# And the individual scores
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    intersectBed -a '$filename' -b methylation/human/repeats_of_interest.bed -wa -wb >methylation/human/repeats_with_scores_'$ind'.bed
    #simplify names for plotting 
    cat methylation/human/repeats_with_scores_'$ind'.bed  |sed 's/-subunit_/*/'|sed 's/_family-/*/'  >methylation/human/repeats_with_scores_'$ind'_renamed.bed
    '| sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done 

# Get the average methylation for G4s overlapping with any of those
# Here we are interested in each G4s separately, so we use the original
# unmerged bed file 
nonB_path=/path/to/original/GQfile
# inhouse: nonB_path="../../shared/nonB_annotation/T2Tv2_primates/"
mkdir methylation/human/overlap
# To get full G4s that overlaps we need to use the flag -wa (save Quadron scores as well)
intersectBed -a <(cat $nonB_path/chm13v2.0.GQ.bed) -b methylation/human/repeats_of_interest.bed  -wao |awk -v OFS="\t" '($11!=0){print $1,$2,$3,$10,$4}' >methylation/human/overlap/G4s_in_roi.bed
# And methylation scores
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    intersectBed -a '$filename' -b methylation/human/overlap/G4s_in_roi.bed -wa -wb >methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed
    # simplify names for plotting
    cat methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed |sed 's/-subunit_/*/' |sed 's/_family-/*/' >methylation/human/overlap/G4s_in_roi_with_scores_'$ind'_renamed.bed
    '| sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done 

# Go through the repeat classes and calculate average methylation for G4s 
# in repeats of interest
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Repeat Avg_meth CpG_sites" | sed "s/ /\t/g" >methylation/human/methylation_in_G4_roi_'$ind'.txt
    for rep in $(cat T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt)
    do
        echo $rep
        meth=`awk -v r=$rep '"'"'($8==r){print}'"'"' methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed | awk -v sum=0 '"'"'{sum+=$4}END{m=sum/NR; print m, NR}'"'"'`
        echo $rep" "$meth | sed "s/ /\t/g" >>methylation/human/methylation_in_G4_roi_'$ind'.txt
    done
    '| sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done

# Average over full genome
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    awk -v sum=0 -v i=$ind '{sum+=$4}END{m=sum/NR; print i, m, NR}' $filename
done 
# H002 0.672558 33328344
# chm13 0.443678 32498737

# Average in all G4s
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    intersectBed -a $filename -b <(cat nonB_annotation/human_pri/*_GQ.bed) -wa -wb \
    |awk -v sum=0 -v i=$ind '{sum+=$4}END{m=sum/NR; print i, m, NR}'
done
# H002 0.44667 751010
# chm13 0.31006 794497


# Plot with 
Rscript plot_fig6_methylation.R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Look at G4s with low methylation in interesting repeats.
# We will focus on the satellites: LSAU, GSAT, GSATII, SAT-VAR_rnd-6*3554,
# ACRO1 and TAR1. These are found in the file: 
# T2T_primate_nonB/helpfiles/human_satellites_of_interest.txt

# First count G4s with methylation per satellite, and how many of them that has
# more than one or more than 3 CpG sites 
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Satellite G4s_with_CpG G4_with_more_than_1_CpG G4_with_more_than_3_CpG" |sed "s/ /\t/g" >methylation/human/Satellite_G4_summary_'$ind'.txt
    for rep in $(cut -f1 T2T_primate_nonB/helpfiles/human_satellites_of_interest.txt)
    do
        echo $rep
        G4=`awk -v r=$rep '"'"'($8==r){print}'"'"' methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed |cut -f6 |uniq -c |wc -l`
        G4_1=`awk -v r=$rep '"'"'($8==r){print}'"'"' methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed |cut -f6 |uniq -c |awk '"'"'($1>1){print}'"'"' |wc -l`
        G4_3=`awk -v r=$rep '"'"'($8==r){print}'"'"' methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed |cut -f6 |uniq -c |awk '"'"'($1>3){print}'"'"' |wc -l`
        echo $rep $G4 $G4_1 $G4_3| sed "s/ /\t/g" >>methylation/human/Satellite_G4_summary_'$ind'.txt
    done 
    ' | sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done 

# And make lists of all G4s, with their min, max and mean methylation
# along with their quadron score 
module load bedtools/2.31.0
mkdir methylation/human/G4_tables
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo "looking at ind $ind"
    cat T2T_primate_nonB/helpfiles/human_satellites_of_interest.txt |while read -r rep alias
    do
        echo "looking at repeat $rep with alias $alias"
        echo "Chr Start Stop N Mean Min Max" |sed "s/ /\t/g" >methylation/human/G4_tables/${alias}_$ind.txt
        awk -v r=$rep -v chr=0 '($8==r){if(chr==0){chr=$5; s=$6; e=$7; sum=$4; min=$4; max=$4; n=1}else{if(chr==$5 && s==$6){if(min>$4){min=$4}; if(max<$4){max=$4}; sum=sum+$4; n++}else{m=sum/n; print chr":"s"-"e,n,m,min,max; chr=$5; s=$6; e=$7; sum=$4; min=$4; max=$4; n=1}}}END{print chr":"s"-"e,n,m,min,max}' methylation/human/overlap/G4s_in_roi_with_scores_$ind.bed >>methylation/human/G4_tables/${alias}_$ind.txt
    done 
done 

#Making combined table to find low methylated G4s, and add Quadron scores 
 cat T2T_primate_nonB/helpfiles/human_satellites_of_interest.txt |while read -r rep alias
 do
    join <(tail -n+2 methylation/human/G4_tables/${alias}_H002.txt|sort) <(tail -n+2 methylation/human/G4_tables/${alias}_chm13.txt|sort) |sort -k1,1 |join -  <(awk -v OFS="\t" '{s=$2-1; print $1":"$2"-"$3, $5, $1,s,$3}' methylation/human/overlap/G4s_in_roi.bed |sort -k1,1) >methylation/human/G4_tables/${alias}_combined.txt
done

# INHOUSE: ADD SEQUENCES
## Add sequences for the G4 motifs to find motifs to verify experimentally
# From original Quadron files (locally)
for i in {1..22} "X" "Y"
do
    grep "^DATA" ../Kaivans_annotation/non-B-DNA-Annotations/Homo_sapiens/Quadron/chr${i}_out.txt | \
    awk -v c="chr"$i '($5!="NA"){s=$2-1; e=s+4; print c":"s"-"e"\t"$3"\t"$6}' >>methylation/human/G4_with_sequence.txt
done 
# Then, merge with the other tables
cat T2T_primate_nonB/helpfiles/human_satellites_of_interest.txt |while read -r rep alias;
do     
    join methylation/human/G4_tables/${alias}_combined.txt <(sort methylation/human/G4_with_sequence.txt) |sed 's/ /\t/g' >methylation/human/G4_tables/${alias}_seq_combined.txt; 
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ALIGNMENT OF SATELLITE COPIES 

intersectBed -wo -a methylation/human/repeats_of_interest.bed -b methylation/human/overlap/G4s_in_roi.bed |grep LSAU |grep -v COMP |awk '($1=="chr1"){print}' |cut -f1,2,3 |uniq |bedtools getfasta -fi new_sequence/human_pri/t2t.chr1.fa -bed - >methylation/human/test.chr1.LSAU.fasta