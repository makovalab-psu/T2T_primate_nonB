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
intersectBed -a <(cat $nonB_path/chm13v2.0.GQ.bed) -b methylation/human/repeats_of_interest.bed  -wao |awk -v OFS="\t" '($11!=0){print $1,$2,$3,$10,$4}' >methylation/human/overlap/G4s_in_roi.bed
# And methylation scores
cat T2T_primate_nonB/helpfiles/human_methylation.txt |while read -r ind filename
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    #intersectBed -a '$filename' -b methylation/human/overlap/G4s_in_roi.bed -wa -wb >methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed
    # simplify names for plotting
    #cat methylation/human/overlap/G4s_in_roi_with_scores_'$ind'.bed |sed 's/-subunit_/*/' |sed 's/_family-/*/' >methylation/human/overlap/G4s_in_roi_with_scores_'$ind'_renamed.bed
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
mkdir methylation/human/sequences/
mkdir methylation/human/satellites/
module load bedtools/2.31.0

# INHOUSE Test: extracting sequences
# Start by taking LSAUs on chr1 (we decided not to mix chromosomes)
intersectBed -wo -a methylation/human/repeats_of_interest.bed -b methylation/human/overlap/G4s_in_roi.bed |grep LSAU |grep -v COMP |awk '($1=="chr1"){print}' |cut -f1,2,3 |uniq |bedtools getfasta -fi new_sequence/human_pri/t2t.chr1.fa -bed - >methylation/human/test.chr1.LSAU.fasta
# When aligning (locally with MUSCLE in MEGA, or on ROAR with mafft), there are 
# clearly (at least) two groups of sequences with very different lengths
module load mafft/7.481 
mafft methylation/human/test.chr1.LSAU.fasta >testmafft

# Locally in MEGA: MANUALLY assign the sequences into two groups and save as
# new fastas, LSAU1.fas and LSAU2.fas, and sequences in LSAU*.txt


# Getting lengths of satellites and G4s, plus starting positions 
intersectBed -wo -a methylation/human/repeats_of_interest.bed -b methylation/human/overlap/G4s_in_roi.bed |grep LSAU |grep -v COMP  |awk '{l=$3-$2; s=$6-$2; print $1":"$2"-"$3,l,$5":"$6"-"$7,$10,s}' >methylation/human/satellites/LSAU.lengths.txt

# Looking at them in R with plot_G4_satellite_test.R, shows there are no 
# differences between the two Satellite types in terms of lengths, and that 
# for each of them there are three common starting points of G4 within the 
# satellites. 

# trying all LSAUs 
intersectBed -wo -a methylation/human/repeats_of_interest_strand.bed -b methylation/human/overlap/G4s_in_roi.bed |grep LSAU |grep -v COMP |cut -f1,2,3,4,5,6 |uniq |bedtools getfasta -s -fi ref/assemblies/chm13v2.0.fa -bed - >methylation/human/satellites/LSAU.strand.fasta
module load mafft/7.481 
mafft methylation/human/satellites/LSAU.strand.fasta >methylation/human/satellites/LSAU.strand.aligned.fasta

# opening them in MEGA, I find at least three-four distinct types 
# get lengths and methylation for each of these types 
# They are saved in LSAU.Type*.fas

# Merge ALL info we have 
grep ">" methylation/human/satellites/LSAU.Type1-1.fas |sed 's/>//' |sort |join - <(sort methylation/human/satellites/LSAU.lengths.txt) |sort -k3,3 |join -1 3 -2 1 - <(sort methylation/human/G4_tables/LSAU_seq_combined.txt) |sed 's/ /\t/g' | cut -f1-14,18-19 >methylation/human/satellites/LSAU.Type1-1.allInfo.txt
grep ">" methylation/human/satellites/LSAU.Type1-2.fas |sed 's/>//' |sort |join - <(sort methylation/human/satellites/LSAU.lengths.txt) |sort -k3,3 |join -1 3 -2 1 - <(sort methylation/human/G4_tables/LSAU_seq_combined.txt) |sed 's/ /\t/g'| cut -f1-14,18-19 >methylation/human/satellites/LSAU.Type1-2.allInfo.txt
grep ">" methylation/human/satellites/LSAU.Type2.fas |sed 's/>//' |sort |join - <(sort methylation/human/satellites/LSAU.lengths.txt) |sort -k3,3 |join -1 3 -2 1 - <(sort methylation/human/G4_tables/LSAU_seq_combined.txt) |sed 's/ /\t/g'| cut -f1-14,18-19 >methylation/human/satellites/LSAU.Type2.allInfo.txt
grep ">" methylation/human/satellites/LSAU.Type3.fas |sed 's/>//' |sort |join - <(sort methylation/human/satellites/LSAU.lengths.txt) |sort -k3,3 |join -1 3 -2 1 - <(sort methylation/human/G4_tables/LSAU_seq_combined.txt) |sed 's/ /\t/g'| cut -f1-14,18-19 >methylation/human/satellites/LSAU.Type3.allInfo.txt


  1   62 C       LSAU    (0)     498     135     366
  2   58 C       LSAU    (359)   139     2       139
  3   58 +       LSAU    1       139     (359)   150
  4   53 +       LSAU    1       498     (0)     496
  5   49 +       LSAU    1       498     (0)     497
  6   43 C       LSAU    (359)   139     2       149
  7   40 C       LSAU    (0)     498     1       497
  8   31 C       LSAU    (359)   139     1       150
  9   31 +       LSAU    16      139     (359)   122
 10   29 +       LSAU    2       139     (359)   139
 11   25 C       LSAU    (0)     498     1       496
 12   21 +       LSAU    2       139     (359)   149
 13   20 C       LSAU    (359)   139     16      122
 14   12 C       LSAU    (0)     498     1       495
 15   11 +       LSAU    1       498     (0)     493
 16   11 +       LSAU    1       498     (0)     494
 17   11 +       LSAU    1       498     (0)     495

 #Type 1  (MEGA Type3)
2221    0.8     0.8     1.6     chr10   134623348       134623714       (134420)        C       LSAU    Satellite       (0)     498     135    366
# Type 2-6-8-13 (MEGA None - might not have G4s)
618     13.2    0.7     9.3     chr1    128193494       128193644       (120193684)     C       LSAU    Satellite       (359)   139     1       150
# Type 3-9-10-12 (MEGA None - might not have G4s)
672     10.4    2.0     9.3     chr1    128223783       128223931       (120163397)     +       LSAU    Satellite       1       139     (359)   148
# Type 4-5-15-16-17 (MEGA Type 1-1)
2514    16.3    1.2     0.6     chr1    128118759       128119253       (120268075)     +       LSAU    Satellite       1       498     (0)     494
# Type 7 (MEGA None?)
2538    15.4    0.8     0.6     chr1    128131313       128131808       (120255520)     C       LSAU    Satellite       (0)     498     2       495

# AFTER MERGING WITH G4
 1   274 1:498:(0)       +
 2   194 (0):498:135     C
 3   168 (0):498:1       C
 4    44 2:498:(0)       +
 5    36 (0):498:3       C
 6    34 (359):139:2     C
 7    30 (0):498:2       C
 8    24 3:498:(0)       +
 9    23 139:498:(0)     +
 10   19 1:497:(1)       +

 # Type 1-4-8-10 (MEGA Type1-1) fullfwd
chr1    128118758       128119253       LSAU    1:498:(0)       +       chr1    128119087       128119135       Satellite_LSAU  27.33   48 
# Type 2 (MEGA Type3) startrev
chr10   134633274       134633641       LSAU    (0):498:135     C       chr10   134633393       134633442       Satellite_LSAU  31.44   49
# Type 3-5-7 (MEGA Type2) fullrev
chr1    128122820       128123316       LSAU    (0):498:2       C       chr1    128122814       128122837       Satellite_LSAU  7.69    17
# Type 6 (MEGA None) endrev
chr1    128217037       128217187       LSAU    (359):139:2     C       chr1    128217147       128217174       Satellite_LSAU  11.81   27
# Type 9 (MEGA None) endfwd
chr14   3183303         3183662          LSAU    139:498:(0)     +       chr14   3183497 3183538 Satellite_LSAU  33.18   41




# Save full information before merging with G4s 
grep "LSAU" ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |grep -v "COMP" |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$10,$12":"$13":"$14,$9}' >methylation/human/satellites/LSAU.fullRepMask.bed

# Merge with G4 
intersectBed -wo -a methylation/human/satellites/LSAU.fullRepMask.bed -b methylation/human/overlap/G4s_in_roi.bed |awk -v OFS="\t" '{l=$3-$2; s=$8-$2; print $1":"$2"-"$3,l,$5,$6,$7":"$8"-"$9,$12,s}' |uniq |awk '($7>0 && $2-7>$6){print}' > methylation/human/satellites/LSAU.lengths.full.txt

# Look into the biggest 5 groups (full-fwd, start-rev, full-rev, end-rev, and end-fwd)
# fullfwd 1:498:(0) +
sed 's/(//g' methylation/human/satellites/LSAU.lengths.full.txt |sed 's/)//g' | awk '{split($3,s,":"); if($4=="+" && s[1]<5 && s[2]>490 && s[3]<5){print}}' >methylation/human/satellites/LSAU.G4s.fullfwd.txt
#startrev (0):498:135 C and fullrev (0):498:1 try take together
sed 's/(//g' methylation/human/satellites/LSAU.lengths.full.txt |sed 's/)//g' | awk '{split($3,s,":"); if($4=="C" && s[1]<5 && s[2]>490){print}}' >methylation/human/satellites/LSAU.G4s.startrev.txt

# #Add methylation 
for type in "fullfwd" "startrev"
do
   echo "G4ID SatID Sat_len Rep Sat_strand G4_len G4_offset H002_CpG H002_mean H002_min H002_max CHM13_CpG CHM13_mean CHM13_min CHM13_max G4_score G4_strand G4_seq" >methylation/human/satellites/LSAU.G4s.$type.meth.txt
    join -1 5 -2 1 <(sort -k5,5 methylation/human/satellites/LSAU.G4s.$type.txt) <(sort -k1,1 methylation/human/G4_tables/LSAU_seq_combined.txt) | cut -d " " -f1-16,20-21 >>methylation/human/satellites/LSAU.G4s.$type.meth.txt
done 


# Check alignment again 
module load bedtools/2.31.0
module load mafft/7.481 
for type in "fullfwd" "startrev"
do  
    cut -f2 -d" " methylation/human/satellites/LSAU.G4s.$type.meth.txt |tail -n+2 |sed 's/:/\t/' |sed 's/-/\t/' |uniq |bedtools getfasta -fi ref/assemblies/chm13v2.0.fa -bed - >methylation/human/satellites/LSAU.$type.fasta
    mafft methylation/human/satellites/LSAU.$type.fasta >methylation/human/satellites/LSAU.$type.aligned.fasta
done 

# Cluster G4s
mkdir methylation/human/cluster
for type in "fullfwd" "startrev"
do  
    cut -f18 -d" " methylation/human/satellites/LSAU.G4s.$type.meth.txt| tail -n+2 |~/software/starcode/starcode -s -d 8 --seq-id >methylation/human/cluster/LSAU.$type.startcode.out
done
# Get the methylation for the top five 
module load python/3.11.2
echo "Type Cluster G4ID SatID Sat_len Rep Sat_strand G4_len G4_offset H002_CpG H002_mean H002_min H002_max CHM13_CpG CHM13_mean CHM13_min CHM13_max G4_score G4_strand G4_seq" >methylation/human/cluster/LSAU.G4s.meth.txt
for type in "fullfwd" "startrev"
do  
    for i in {1..5}
    do 
        line=`head -n $i methylation/human/cluster/LSAU.$type.startcode.out | tail -n1`
        #echo $line 
        seq=`echo $line |cut -f1 -d" "`
        echo $seq
        rows=`echo $line |cut -f3 -d" "`
        echo $rows
        python T2T_primate_nonB/python/print_lines.py -i methylation/human/satellites/LSAU.G4s.$type.meth.txt -s -r "$rows" | awk -v c=$i -v type=$type '{print type,c,$0}' >>methylation/human/cluster/LSAU.G4s.meth.txt
    done
done 


echo ACCTTGAAA | tr ACGTacgt TGCAtgca | rev

# LOKKING AT SOME OTHER TYPES 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GSAT 
grep "GSAT" ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |grep -v "GSATII" |grep -v "GSATX" |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$10,$12":"$13":"$14,$9}' >methylation/human/satellites/GSAT.fullRepMask.bed
# 39 copies in total, 30 on chr8

# GET FASTA 
module load bedtools/2.31.0
intersectBed -wo -a methylation/human/repeats_of_interest_strand.bed -b methylation/human/overlap/G4s_in_roi.bed |grep GSAT |grep -v GSATII  |grep -v GSATX |cut -f1,2,3,4,5,6 |uniq |bedtools getfasta -s -fi ref/assemblies/chm13v2.0.fa -bed - >methylation/human/satellites/GSAT.strand.fasta
module load mafft/7.481 
mafft methylation/human/satellites/GSAT.strand.fasta >methylation/human/satellites/GSAT.strand.aligned.fasta
# ALIGNMENT FAILS! SEQ ARE TOO DIFFERENT I THINK, VERY DIFFEREN LENGHTS - from 11bp to 26000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GSATII
grep "GSATII" ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out  |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$10,$12":"$13":"$14,$9}' >methylation/human/satellites/GSATII.fullRepMask.bed
# 174 in total, 102 on chr12. Length vary from 30 to 7000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACRO1
grep "ACRO1" ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out  |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$10,$12":"$13":"$14,$9}' >methylation/human/satellites/ACRO1.fullRepMask.bed
# 488 in total, 158 on chr22, 88 on chr15, 85 on chr14, 65 on chr21.  Length vary from 160 to 5000
module load bedtools/2.31.0
intersectBed -wo -a methylation/human/repeats_of_interest_strand.bed -b methylation/human/overlap/G4s_in_roi.bed |grep "ACRO1" |cut -f1,2,3,4,5,6 |uniq |bedtools getfasta -s -fi ref/assemblies/chm13v2.0.fa -bed - >methylation/human/satellites/ACRO1.strand.fasta
module load mafft/7.481 
mafft methylation/human/satellites/ACRO1.strand.fasta >methylation/human/satellites/ACRO1.strand.aligned.fasta
# Seems like most align well, apart from a handful or so sequences that are much longer than everything else. 

# Merge with G4 
intersectBed -wo -a methylation/human/satellites/ACRO1.fullRepMask.bed -b methylation/human/overlap/G4s_in_roi.bed |awk -v OFS="\t" '{l=$3-$2; s=$8-$2; print $1":"$2"-"$3,l,$5,$6,$7":"$8"-"$9,$12,s}' |uniq |awk '($7>0 && $2-7>$6){print}' > methylation/human/satellites/ACRO1.lengths.full.txt

# Divide into forward and reverse
rep="ACRO1"  
sed 's/(//g' methylation/human/satellites/$rep.lengths.full.txt |sed 's/)//g' | awk '{split($3,s,":"); if($4=="+"){print}}' >methylation/human/satellites/$rep.G4s.fwd.txt
sed 's/(//g' methylation/human/satellites/$rep.lengths.full.txt |sed 's/)//g' | awk '{split($3,s,":"); if($4=="C"){print}}' >methylation/human/satellites/$rep.G4s.rev.txt

# #Add methylation 
rep="ACRO1" 
for type in "fwd" "rev"
do
   echo "G4ID SatID Sat_len Rep Sat_strand G4_len G4_offset HG002_CpG HG002_mean HG002_min HG002_max CHM13_CpG CHM13_mean CHM13_min CHM13_max G4_score G4_strand G4_seq" >methylation/human/satellites/$rep.G4s.$type.meth.txt
    join -1 5 -2 1 <(sort -k5,5 methylation/human/satellites/$rep.G4s.$type.txt) <(sort -k1,1 methylation/human/G4_tables/${rep}_seq_combined.txt) | cut -d " " -f1-16,20-21 >>methylation/human/satellites/$rep.G4s.$type.meth.txt
done 

# Cluster G4s
rep="ACRO1" 
for type in "fwd" "rev"
do  
    cut -f18 -d" " methylation/human/satellites/$rep.G4s.$type.meth.txt| tail -n+2 |~/software/starcode/starcode -s -d 8 --seq-id >methylation/human/cluster/$rep.$type.startcode.out
done
# Get the methylation for the top five 
rep="ACRO1" 
module load python/3.11.2
echo "Type Cluster G4ID SatID Sat_len Rep Sat_strand G4_len G4_offset HG002_CpG HG002_mean HG002_min HG002_max CHM13_CpG CHM13_mean CHM13_min CHM13_max G4_score G4_strand G4_seq" >methylation/human/cluster/$rep.G4s.meth.txt
for type in "rev" #"fwd" 
do  
    for i in {1..2}
    do 
        line=`head -n $i methylation/human/cluster/$rep.$type.startcode.out | tail -n1`
        #echo $line 
        seq=`echo $line |cut -f1 -d" "`
        echo $seq
        rows=`echo $line |cut -f3 -d" "`
        echo $rows
        python T2T_primate_nonB/python/print_lines.py -i methylation/human/satellites/$rep.G4s.$type.meth.txt -s -r "$rows" | awk -v c=$i -v type=$type '{print type,c,$0}' >>methylation/human/cluster/$rep.G4s.meth.txt
    done
done 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SAT-VAR_rnd-6*3554
grep "SAT-VAR" ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |grep "family-3554"  |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$10,$12":"$13":"$14,$9}' >methylation/human/satellites/SAT-VAR.fullRepMask.bed
# 115 in total, 23 on chr19.  Length vary from 35 to 4754

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TAR1
grep "TAR1" ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$10,$12":"$13":"$14,$9}' >methylation/human/satellites/TAR1.fullRepMask.bed
# 201 in total, 25 on chr9.  Length vary from 13 to 2747

#############################################################
# UP NEXT
# DO TWO SEPARATE RUNS FOR LSAU WITH THE CHORMOSOMES THAT LOOK DIFFERENT: 
# chr10 and chr4. 
awk '(NR==1 || /chr10:/){print}' methylation/human/satellites/LSAU.G4s.startrev.meth.txt > methylation/human/satellites/LSAU.chr10.G4s.startrev.meth.txt
awk '(NR==1 || /chr4:/){print}' methylation/human/satellites/LSAU.G4s.startrev.meth.txt > methylation/human/satellites/LSAU.chr4.G4s.startrev.meth.txt
grep -v "chr10:" methylation/human/satellites/LSAU.G4s.startrev.meth.txt |grep -v "chr4:" >methylation/human/satellites/LSAU.excl4-10.G4s.startrev.meth.txt

for chr in "chr10" "chr4" "excl4-10" 
do  
    cut -f18 -d" " methylation/human/satellites/LSAU.$chr.G4s.startrev.meth.txt| tail -n+2 |~/software/starcode/starcode -s -d 8 --seq-id >methylation/human/cluster/LSAU.$chr.starcode.out
done
# Chr4 and 10 had three main clusters, the rest looks messier as expected. 
module load python/3.11.2
echo "Type Cluster G4ID SatID Sat_len Rep Sat_strand G4_len G4_offset HG002_CpG HG002_mean HG002_min HG002_max CHM13_CpG CHM13_mean CHM13_min CHM13_max G4_score G4_strand G4_seq" >methylation/human/cluster/LSAU.sepchr.G4s.meth.txt
for chr in  "chr10" "chr4" "excl4-10" 
do  
    for i in {1..3}
    do 
        line=`head -n $i methylation/human/cluster/LSAU.$chr.starcode.out | tail -n1`
        #echo $line 
        seq=`echo $line |cut -f1 -d" "`
        echo $seq
        rows=`echo $line |cut -f3 -d" "`
        echo $rows
        python T2T_primate_nonB/python/print_lines.py -i methylation/human/satellites/LSAU.$chr.G4s.startrev.meth.txt -s -r "$rows" | awk -v c=$i -v type=$chr '{print type,c,$0}' >>methylation/human/cluster/LSAU.sepchr.G4s.meth.txt
    done
done 

# RERUN ANALYSIS FOR LSAU AND ACRO1 IN A SINGLE LOOP, USING ALL REV AND FWD

