################################################################################
# CODE FOR INVESTIGATING METHYLATION PATTERNS IN REPEATS AND G4 MOTIFS

# # # # # TABLE OF CONTENTS 
# - METHYLATION IN SATELLITES AND OTHER REPEATS
# - DEEPER INVESTIGATION OF CERTAIN REPEATS 
# - MORE DETAILED STUDY OF LSAU 
# - INCORPORATING ANOTHER EXPERIMENTAL DATASET 

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
# Update this to add strands 
awk -v OFS="\t" '{if(/Satellite/){name=$11"_"$10}else{name=$11}; s=$6-1; if($9=="C"){$9="-"}; print $5,s,$7,name,"0",$9}' ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out >repeats/human/RepeatMasker_strand.bed
python3 T2T_primate_nonB/python/extract_from_bed.py -b repeats/human/RepeatMasker_strand.bed -l T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt -c 4 -o methylation/human/repeats_of_interest_strand.bed

# Get average methylation for each of these repeat categories
module load bedtools/2.31.0
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo "Repeat Avg_meth CpG_sites" | sed 's/ /\t/g' >methylation/human/methylation_in_roi_$ind.txt
    cat T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt |while read -r rep repshort
    do
    echo $rep
    meth=`awk -v r=$rep '($4==r){print}' methylation/human/repeats_of_interest.bed |intersectBed -a $filename -b - |awk -v sum=0 '{sum+=$4}END{m=sum/NR; print m, NR}'`
    echo $rep" "$meth | sed 's/ /\t/g' >>methylation/human/methylation_in_roi_$ind.txt
    done
done 

# And the individual scores
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    #intersectBed -a '$filename' -b methylation/human/repeats_of_interest.bed -wa -wb >methylation/human/repeats_with_scores_'$ind'.bed
    #simplify names for plotting 
    #cat methylation/human/repeats_with_scores_'$ind'.bed  |sed 's/-subunit_/*/'|sed 's/_family-/*/'  >methylation/human/repeats_with_scores_'$ind'_renamed.bed
    cat methylation/human/repeats_with_scores_'$ind'.bed  |sort -k8,8 |join -1 8 -2 1 - <(sort T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt)|awk -v OFS="\t" '"'"'{print $2,$3,$4,$5,$6,$7,$8,$1,$9}'"'"' >methylation/human/repeats_with_scores_'$ind'_shortnames.bed
    '| sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done 

# Get the average methylation for G4s overlapping with any of those
# Here we are interested in each G4s separately, so we use the original
# unmerged bed file 
nonB_path=/path/to/original/GQfile
# inhouse: nonB_path="../../shared/nonB_annotation/T2Tv2_primates/"
mkdir methylation/human/overlap
# To get full G4s that overlaps we need to use the flag -wa (save Quadron scores as well)
# Sometimes there are overlaps in the repeat anntotation that lead to duplicated
# lines after intersecting -> only save each position once!
intersectBed -a <(cat $nonB_path/chm13v2.0.GQ.bed) -b methylation/human/repeats_of_interest.bed  -wao |awk -v OFS="\t" '($11!=0){print $1,$2,$3,$10,$4}' |uniq  >methylation/human/overlap/G4s_in_roi.bed
# And methylation scores
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    intersectBed -a '$filename' -b methylation/human/overlap/G4s_in_roi.bed -wa -wb >methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed
    # simplify names for plotting
    cat methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed |sort -k8,8 |join -1 8 -2 1 - <(sort T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt)|awk -v OFS="\t" '"'"'{print $2,$3,$4,$5,$6,$7,$8,$9,$1,$10}'"'"' >methylation/human/overlap/G4s_in_roi_with_scores_'$ind'_shortnames.bed
    '| sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done 

# Go through the repeat classes and calculate average methylation for G4s 
# in repeats of interest
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Repeat Avg_meth CpG_sites" | sed "s/ /\t/g" >methylation/human/methylation_in_G4_roi_'$ind'.txt
    cat T2T_primate_nonB/helpfiles/human_repeats_of_interest.txt |while read -r rep repshort
    do
        echo $rep
        meth=`awk -v r=$rep '"'"'($8==r){print}'"'"' methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed | awk -v sum=0 '"'"'{sum+=$4}END{m=sum/NR; print m, NR}'"'"'`
        echo $rep" "$meth | sed "s/ /\t/g" >>methylation/human/methylation_in_G4_roi_'$ind'.txt
    done
    '| sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done

# Numbers for Table S3 are found in the files:
#methylation/human/methylation_in_G4_roi_chm13.txt  
#methylation/human/methylation_in_roi_chm13.txt
#methylation/human/methylation_in_G4_roi_HG002.txt   
#methylation/human/methylation_in_roi_HG002.txt

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
    intersectBed -a $filename -b <(cat nonB_annotation/human_pri/genome_GQ.bed) -wa -wb |cut -f1-4 \
    |awk -v sum=0 -v i=$ind '{sum+=$4}END{m=sum/NR; print i, m, NR}'
done
# H002 0.44667 751010
# chm13 0.31006 794497


# Plot with 
Rscript plot_fig6ABC_repeats_methylation.R

# Noted that for some repeats (like SVA), CHM13 has many more methylated sites
# than HG002. Intersected them with -v to get everything NOT in HG002, and a 
# great part of it came from a region on a single chromosome (chr2). This 
# region doesn't have any methylated sites annotated for HG002, so most likely
# this is due to true differences (for example a deletion in HG002, or an 
# expansion in chm13)
grep SVA methylation/human/repeats_with_scores_H002.bed |intersectBed -v -a <(grep SVA methylation/human/repeats_with_scores_chm13.bed) -b - >temp_meth_in_chm13_nt_in_HG002_SVA.txt
cut -f1 temp_meth_in_chm13_nt_in_HG002_SVA.txt |uniq -c



################### DEEPER INVESTIGATION OF CERTAIN REPEATS ###################
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

# ADD SEQUENCES 
## Add sequences for the G4 motifs to find motifs to verify experimentally
# From original Quadron files 
for i in {1..22} "X" "Y"
do
    grep "^DATA" path/to/Quadron/chr${i}_out.txt | \
    awk -v c="chr"$i '($5!="NA"){s=$2-1; e=s+4; print c":"s"-"e"\t"$3"\t"$6}' >>methylation/human/G4_with_sequence.txt
done 
# Then, merge with the other tables
cat T2T_primate_nonB/helpfiles/human_satellites_of_interest.txt |while read -r rep alias;
do     
    join methylation/human/G4_tables/${alias}_combined.txt <(sort methylation/human/G4_with_sequence.txt) |sed 's/ /\t/g' >methylation/human/G4_tables/${alias}_seq_combined.txt; 
done


######################### MORE DETAILED STUDY OF LSAU #########################

# LSAU is the satellite with most G4s, so we focus on this. 

# Save full information before merging with G4s 
grep "LSAU" ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |grep -v "COMP" |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$10,$12":"$13":"$14,$9}' >methylation/human/satellites/LSAU.fullRepMask.bed

# Merge with G4 (only full overlaps)
intersectBed -wo -a methylation/human/satellites/LSAU.fullRepMask.bed -b methylation/human/overlap/G4s_in_roi.bed |awk -v OFS="\t" '{l=$3-$2; s=$8-$2; print $1":"$2"-"$3,l,$5,$6,$7":"$8"-"$9,$12,s}' |uniq |awk '($7>0 && $2-7>$6){print}' > methylation/human/satellites/LSAU.lengths.full.txt

# Align full LSAUs (only those that overlap with G4)
module load bedtools/2.31.0
module load mafft/7.481 
cut -f2 -d" " methylation/human/satellites/LSAU.lengths.full.txt |tail -n+2 |sed 's/:/\t/' |sed 's/-/\t/' |uniq |bedtools getfasta -s -fi ref/assemblies/chm13v2.0.fa -bed - >methylation/human/satellites/LSAU.fasta
mafft methylation/human/satellites/LSAU.fasta >methylation/human/satellites/LSAU.aligned.fasta

# We see several groups, the satellites on chr10 and 4 (on the minus strand) 
# are very similar to eachother, while LSAU on the other chromosomes varies a
# bit more (copies look more diverged, they don't always contain the full consensus
# sequence, etc)

# We first separate the G4s per strand
awk '($4=="+"){print}' methylation/human/satellites/LSAU.lengths.full.txt >methylation/human/satellites/LSAU.G4s.fwd.txt
awk '($4=="C"){print}' methylation/human/satellites/LSAU.lengths.full.txt >methylation/human/satellites/LSAU.G4s.rev.txt

# #Add methylation and sequence
for type in "fwd" "rev"
do
   echo "G4ID SatID Sat_len Rep Sat_strand G4_len G4_offset H002_CpG H002_mean H002_min H002_max CHM13_CpG CHM13_mean CHM13_min CHM13_max G4_score G4_strand G4_seq" >methylation/human/satellites/LSAU.G4s.$type.meth.txt
    join -1 5 -2 1 <(sort -k5,5 methylation/human/satellites/LSAU.G4s.$type.txt) <(sort -k1,1 methylation/human/G4_tables/LSAU_seq_combined.txt) | cut -d " " -f1-16,20-21 >>methylation/human/satellites/LSAU.G4s.$type.meth.txt
done 

# Cluster G4s
mkdir methylation/human/cluster
for type in "fwd" "rev"
do  
    cut -f18 -d" " methylation/human/satellites/LSAU.G4s.$type.meth.txt| tail -n+2 |~/software/starcode/starcode -s -d 8 --seq-id >methylation/human/cluster/LSAU.$type.startcode.out
done

# Get the methylation for the top five 
module load python/3.11.2
echo "Type Cluster G4ID SatID Sat_len Rep Sat_strand G4_len G4_offset H002_CpG H002_mean H002_min H002_max CHM13_CpG CHM13_mean CHM13_min CHM13_max G4_score G4_strand G4_seq" >methylation/human/cluster/LSAU.G4s.meth.fwdrev.txt
for type in "fwd" "rev"
do  
    for i in {1..5}
    do 
        line=`head -n $i methylation/human/cluster/LSAU.$type.startcode.out | tail -n1`
        #echo $line 
        seq=`echo $line |cut -f1 -d" "`
        echo $seq
        rows=`echo $line |cut -f3 -d" "`
        echo $rows
        python T2T_primate_nonB/python/print_lines.py -i methylation/human/satellites/LSAU.G4s.$type.meth.txt -s -r "$rows" | awk -v c=$i -v type=$type '{print type,c,$0}' >>methylation/human/cluster/LSAU.G4s.meth.fwdrev.txt
    done
done 

# This makes quite large clusters, some with quite diverged sequences
# Since chr4 and 10 show very distinct G4s, we run the clustering 
# separately on these chromosomes 

# Separate chr10, chr4, rest 
for type in "fwd" "rev"
do  
    awk '(NR==1 || /chr10:/){print}' methylation/human/satellites/LSAU.G4s.$type.meth.txt > methylation/human/satellites/LSAU.chr10.G4s.$type.meth.txt
    awk '(NR==1 || /chr4:/){print}' methylation/human/satellites/LSAU.G4s.$type.meth.txt > methylation/human/satellites/LSAU.chr4.G4s.$type.meth.txt
    grep -v "chr10:" methylation/human/satellites/LSAU.G4s.$type.meth.txt |grep -v "chr4:" >methylation/human/satellites/LSAU.excl4-10.G4s.$type.meth.txt
done 

module load python/3.11.2
for type in "fwd" "rev"
do
    echo "Type Cluster G4ID SatID Sat_len Rep Sat_strand G4_len G4_offset HG002_CpG HG002_mean HG002_min HG002_max CHM13_CpG CHM13_mean CHM13_min CHM13_max G4_score G4_strand G4_seq" >methylation/human/cluster/LSAU.sepchr.$type.G4s.meth.txt
    for chr in "chr10" "chr4" "excl4-10" 
    do  
        cut -f18 -d" " methylation/human/satellites/LSAU.$chr.G4s.$type.meth.txt| tail -n+2 |~/software/starcode/starcode -s -d 8 --seq-id >methylation/human/cluster/LSAU.$chr.$type.starcode.out
        # Extract first three clusters from each type 
        for i in {1..5}
        do 
            line=`head -n $i methylation/human/cluster/LSAU.$chr.$type.starcode.out | tail -n1`
            #echo $line 
            seq=`echo $line |cut -f1 -d" "`
            echo $seq
            rows=`echo $line |cut -f3 -d" "`
            echo $rows
            python T2T_primate_nonB/python/print_lines.py -i methylation/human/satellites/LSAU.$chr.G4s.$type.meth.txt -s -r "$rows" | awk -v c=$i -v type=$chr '{print type,c,$0}' >>methylation/human/cluster/LSAU.sepchr.$type.G4s.meth.txt
        done 
    done 
done

# Plot methylation score and quadron score distribution 
# plot_FigS2_methylation_LSAU.R 


################## INCORPORATING ANOTHER EXPERIMENTAL DATASET ##################
# TRY INCORPORATING OTHER EXPERIMENTAL DATA, BG PEAKS FROM HU ET AL. 2021

mkdir G4peaks
cd G4peaks

# Download the data from
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181373


# Translate peaks from hg38 to chm13
~/software/liftOver GSE181373_K562_100k_bg4_bulk_1000.bio.bed ~/software/liftOver_chains/hg38ToHs1.over.chain GSE181373_K562.hs1liftOver.bed GSE181373_K562.hs1liftOver.unmapped.bed
~/software/liftOver GSE181373_U2OS_100k_bg4_bulk_1000.bio.bed ~/software/liftOver_chains/hg38ToHs1.over.chain GSE181373_U2OS.hs1liftOver.bed GSE181373_U2OS.hs1liftOver.unmapped.bed

# Check numbers
wc GSE181373_K562_100k_bg4_bulk_1000.bio.bed
#   18365   55095  437847 GSE181373_K562_100k_bg4_bulk_1000.bio.bed
wc GSE181373_K562.hs1liftOver.bed
#   18290   54870  436084 GSE181373_K562.bio.hs1liftOver.bed
wc GSE181373_U2OS_100k_bg4_bulk_1000.bio.bed
#   35452  106356  847398 GSE181373_U2OS_100k_bg4_bulk_1000.bio.bed
wc GSE181373_U2OS.hs1liftOver.bed
#   35355  106065  845381 GSE181373_U2OS.bio.hs1liftOver.bed
# Almost all peaks can be lift over.

# Length:
awk '{sum+=$3-$2}END{m=sum/NR; print m}' GSE181373_K562_100k_bg4_bulk_1000.bio.bed
#1348.11
awk '{sum+=$3-$2}END{m=sum/NR; print m}' GSE181373_K562.hs1liftOver.bed
#1367.1
awk '{sum+=$3-$2}END{m=sum/NR; print m}' GSE181373_U2OS_100k_bg4_bulk_1000.bio.bed
#1013.13
awk '{sum+=$3-$2}END{m=sum/NR; print m}' GSE181373_U2OS.hs1liftOver.bed
#1020.9
# Peaks are shorter after liftover. 

# Overlap peaks with Quadron
intersectBed -a GSE181373_K562.hs1liftOver.bed -b <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) -wo |cut -f1-3 |uniq |wc
#   13377   40131  318786
intersectBed -a GSE181373_U2OS.hs1liftOver.bed -b <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) -wo |cut -f1-3 |uniq |wc
#   18366   55098  438592
# => Clearly more than half for the first cell, and ~half for the second cell line.

# How many G4s are included?
intersectBed -a GSE181373_K562.hs1liftOver.bed -b <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) -wo |cut -f5-7 |uniq |wc
#   50795  152385 1079781
intersectBed -a GSE181373_U2OS.hs1liftOver.bed -b <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) -wo |cut -f5-7 |uniq |wc
#   57194  171582 1220448
#Quite a lot! Meaning each peak contains many G4s

# # # ~~~~~~~~~~~~~~~~~~~~
# Overlap with repeats of interests
# ~~~~~~~~~~~~~~ Cell type K562
intersectBed -a GSE181373_K562.hs1liftOver.bed -b ../methylation/human/repeats_of_interest.bed -wo |cut -f1-3 |uniq |wc
#   143     429    3340
# With how many repeat regions?
intersectBed -a GSE181373_K562.hs1liftOver.bed -b ../methylation/human/repeats_of_interest.bed -wo |cut -f4-6 |uniq |wc
#   198     594    4601
# From how many Repeat types?
intersectBed -a GSE181373_K562.hs1liftOver.bed  -b ../methylation/human/repeats_of_interest.bed -wo |cut -f7 |sort |uniq
#Retroposon/SVA
#Satellite/acro_ACRO1
#Satellite/centr_GSATII
#Satellite/centr_SST1
#Satellite/subtelo_TAR1
#Satellite_COMP-subunit_LSAU-BSAT_rnd-1_family-1
#Satellite_COMP-subunit_LSAU-BSAT_rnd-1_family-4
#Satellite_COMP-subunit_LSAU-BSAT_rnd-6_family-5403
#Satellite_HSAT5
#Satellite_LSAU
#Satellite_MSR1
#Satellite_SAT-VAR_rnd-6_family-3554
#rRNA

# But does any of these (overlapping regions, not full peaks!) have G4s in them?
# (This gives the overlap between overlapping regions and G4s, not the number
# of regions)
intersectBed -a GSE181373_K562.hs1liftOver.bed -b ../methylation/human/repeats_of_interest.bed |intersectBed -a - -b <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) |uniq |wc
#     146     438    3358

# How many G4 per repeat type?
intersectBed -b GSE181373_K562.hs1liftOver.bed -a ../methylation/human/repeats_of_interest.bed |intersectBed -a - -b <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) |cut -f4 |sort |uniq |wc
# 2 Retroposon/SVA
# 10 Satellite/acro_ACRO1
# 7 Satellite/centr_GSATII
# 12 Satellite/centr_SST1
# 6 Satellite/subtelo_TAR1
# 19 Satellite_COMP-subunit_LSAU-BSAT_rnd-1_family-1
# 7 Satellite_COMP-subunit_LSAU-BSAT_rnd-1_family-4
# 5 Satellite_HSAT5
# 12 Satellite_LSAU
# 47 Satellite_MSR1
# 3 Satellite_SAT-VAR_rnd-6_family-3554
# 16 rRNA


# ~~~~~~~~~~~~~~ Cell type U20S
intersectBed -a GSE181373_U2OS.hs1liftOver.bed -b ../methylation/human/repeats_of_interest.bed -wo |cut -f1-3 |uniq |wc
#     164     492    3861
# With how many repeat regions?
intersectBed -a GSE181373_U2OS.hs1liftOver.bed -b ../methylation/human/repeats_of_interest.bed -wo |cut -f4-6 |uniq |wc
#      206     618    4850
# From how many Repeat types?
intersectBed -a GSE181373_U2OS.hs1liftOver.bed  -b ../methylation/human/repeats_of_interest.bed -wo |cut -f7 |sort |uniq
#Retroposon/SVA
#Satellite/acro_ACRO1
#Satellite/centr_GSATII
#Satellite/centr_SST1
#Satellite/subtelo_TAR1
#Satellite_COMP-subunit_LSAU-BSAT_rnd-1_family-1
#Satellite_COMP-subunit_LSAU-BSAT_rnd-1_family-4
#Satellite_COMP-subunit_LSAU-BSAT_rnd-6_family-5403
#Satellite_HSAT5
#Satellite_LSAU
#Satellite_MSR1
#Satellite_SAT-VAR_rnd-6_family-3554
#rRNA

# Do they have G4s in them?
intersectBed -a GSE181373_U2OS.hs1liftOver.bed -b ../methylation/human/repeats_of_interest.bed |intersectBed -a - -b <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) |uniq |wc
#     162     486    3805


# ~~~~~~~~~~ COMPARING CHM13 WITH HG38
# How many of Repeats of interest were actually annotated in Hg38?
~/software/liftOver ../methylation/human/repeats_of_interest.bed ~/software/liftOver_chains/hs1ToHg38.over.chain repeats_of_interests.Hg38liftOver.bed repeats_of_interests.Hg38liftOver.unmapped.bed

# In CHM13:
wc ../methylation/human/repeats_of_interest.bed
#15592   62368  630421 ../methylation/human/repeats_of_interest.bed
sort -k1,1 -k2,2n ../methylation/human/repeats_of_interest.bed |uniq |awk '{sum+=$3-$2}END{m=sum/NR; print m}'
#658.368
sort -k1,1 -k2,2n ../methylation/human/repeats_of_interest.bed |uniq |awk '{sum+=$3-$2}END{print sum}'
#10265273

# In HG38
wc repeats_of_interests.Hg38liftOver.bed
# 13730   54920  570085 repeats_of_interests.Hg38liftOver.bed

sort -k1,1 -k2,2n repeats_of_interests.Hg38liftOver.bed |uniq |awk '{sum+=$3-$2}END{m=sum/NR; print m}'
# 555.618
sort -k1,1 -k2,2n repeats_of_interests.Hg38liftOver.bed |uniq |awk '{sum+=$3-$2}END{print sum}'
# 6269592
# Check regions on Un chromosomes
cut -f1-3 repeats_of_interests.Hg38liftOver.bed |grep "_"|wc
# 1572    4716   55532
# That sum up to
cut -f1-3 repeats_of_interests.Hg38liftOver.bed |grep "_"|awk '{sum+=$3-$2}END{print sum}'
# 1074571
 # Meaning about half of the total repeat length is not found in the anchored
 # chromosomes of Hg38.


# Checking lengths before and after liftover - only include anchored chromosomes
# and remove overlap after liftover
for i in $(cat repeats_of_interest.txt)
do
  before=`awk -v rep=$i '($4==rep){print}' ../methylation/human/repeats_of_interest.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '{sum+=$3-$2}END{print sum}'`
  after=`awk -v rep=$i '($1!~/_/ && $4==rep){print}' repeats_of_interests.Hg38liftOver.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '{sum+=$3-$2}END{print sum}'`
  echo $i" "$before" "$after |sed 's/ /\t/g'
done
# Quite a lot of overlapping sequences, I'll save the new merged bed files:
for i in $(cat repeats_of_interest.txt)
do
  newname=`echo $i |sed 's/\//-/g'`
  echo $newname
  awk -v rep=$i '($1!~/_/ && $4==rep){print}' repeats_of_interests.Hg38liftOver.bed |sort -k1,1 -k2,2n |mergeBed -i - >$newname.Hg38liftOver.bed
done


# Check number of G4 before and after liftover - need to liftover G4s first!!
~/software/liftOver <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) ~/software/liftOver_chains/hs1ToHg38.over.chain G4s.Hg38liftOver.bed G4s.Hg38liftOver.unmapped.bed

# This overlaps the chm13 repeats with our own G4 annotation, and the translated
# hg38 repeats with the translated G4 (that potentially can be many-to-one), and
# returns the number of unique G4s.
for i in $(cat repeats_of_interest.txt)
do
  before=`awk -v rep=$i '($4==rep){print}' ../methylation/human/repeats_of_interest.bed |intersectBed -a - -b <(cat ../nonB_annotation/Homo_sapiens/autosomes_GQ.bed ../nonB_annotation/Homo_sapiens/chrX_GQ.bed ../nonB_annotation/Homo_sapiens/chrY_GQ.bed) | sort |uniq |wc -l`
  newname=`echo $i |sed 's/\//-/g'`
  after=`intersectBed -a $newname.Hg38liftOver.bed -b G4s.Hg38liftOver.bed |sort |uniq |wc -l`
  echo $i" "$before" "$after |sed 's/ /\t/g'
done

# The last numbers that would be good to have is the number of peaks for each
# repeat type (and cell type).
for i in $(cat repeats_of_interest.txt)
do
  c1=`awk -v rep=$i '($4==rep){print}' ../methylation/human/repeats_of_interest.bed |intersectBed -a - -b GSE181373_K562.hs1liftOver.bed |grep -v "chrM" |sort |uniq |wc -l`
  c2=`awk -v rep=$i '($4==rep){print}' ../methylation/human/repeats_of_interest.bed |intersectBed -a - -b GSE181373_U2OS.hs1liftOver.bed |grep -v "chrM" |sort |uniq |wc -l`
    echo $i"    "$c1"    "$c2
done

for i in $(cat repeats_of_interest.txt)
do
  c1=`awk -v rep=$i '($4==rep){print}' repeats_of_interests.Hg38liftOver.bed |intersectBed -a - -b  GSE181373_K562_100k_bg4_bulk_1000.bio.bed |grep -v "chrM" |sort |uniq |wc -l`
  c2=`awk -v rep=$i '($4==rep){print}' repeats_of_interests.Hg38liftOver.bed |intersectBed -a - -b  GSE181373_U2OS_100k_bg4_bulk_1000.bio.bed|grep -v "chrM" |sort |uniq |wc -l`
    echo $i"    "$c1"    "$c2
done
