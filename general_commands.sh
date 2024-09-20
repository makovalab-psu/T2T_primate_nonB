# CODE FOR THE PAPER 
# Non-canonical DNA in human and other ape telomere-to-telomere genomes
# https://www.biorxiv.org/content/10.1101/2024.09.02.610891v1

# Non-B annotation can be downloaded from:




# CONTENTS 

####################### NOTES AND SOFTWARE REQUIREMENTS ######################## 
# Some code is written to be run as parallel batch jobs on a slurm cluster for 
# speed and efficiency. These could be run in the terminal istead by removing 
# the initial "echo '#!/bin/bash" and the trailing  "| sbatch [..] "

# Software used:
# bedtools v2.31.0
# meryl v1.4.1
# Winnowmap v2.03
# Bedops v2.4.41
# UCSC tools:
#   bigWigToWig


############ CREATE NON-OVERLAPPING BED FILES FROM NON-B ANNOTATION ############
nb_path=/path/to/nonB/tracks/converted/to/bed
#INHOUSE nb_path="/storage/group/kdm16/default/shared/nonB_annotation/T2Tv2_primates/"
#INHOUSE module load bedtools/2.31.0
mkdir nonB_annotation 

# Primary haplotype assembly 
cat T2T_primate_nonB/helpfiles/primary_species_list.txt |grep -v "chimp" | while read -r trivial latin filename;
do
    echo "Creating directory for $trivial.."
    mkdir -p nonB_annotation/${trivial}_pri
    prefix=`echo $filename | sed 's/.fasta//' | sed 's/.fa//'`
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
    do
        #echo $prefix
        #ls $nb_path/$prefix.$nb.bed
        grep "chrX" $nb_path/$prefix.$nb.bed |mergeBed -i - >nonB_annotation/${trivial}_pri/chrX_$nb.bed
        grep "chrY" $nb_path/$prefix.$nb.bed |mergeBed -i - >nonB_annotation/${trivial}_pri/chrY_$nb.bed
        grep -v "chrX" $nb_path/$prefix.$nb.bed |grep -v "chrY" \
        |mergeBed -i - >nonB_annotation/${trivial}_pri/autosomes_$nb.bed
    done 
   for chr in "autosomes" "chrX" "chrY"
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        cat nonB_annotation/'${trivial}'_pri/'${chr}'_*.bed |sort -k1,1 -k2,2n \
        |mergeBed -i - >nonB_annotation/'${trivial}'_pri/'${chr}'_all.bed
        ' | sbatch -J ${chr} --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G
    done
done 

# Secondary haplotype assembly 
cat T2T_primate_nonB/helpfiles/alt_species_list.txt |grep siamang | while read -r trivial latin filename;
do
    echo "Creating directory for $trivial.."
    mkdir -p nonB_annotation/${trivial}_alt
     prefix=`echo $filename | sed 's/.fasta//' | sed 's/.fa//'`
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
    do
        grep "chrX" $nb_path/$prefix.$nb.bed |mergeBed -i - >nonB_annotation/${trivial}_alt/chrX_$nb.bed
        grep "chrY" $nb_path/$prefix.$nb.bed |mergeBed -i - >nonB_annotation/${trivial}_alt/chrY_$nb.bed
        grep -v "chrX" $nb_path/$prefix.$nb.bed |grep -v "chrY" \
        |mergeBed -i - >nonB_annotation/${trivial}_alt/autosomes_$nb.bed
    done
    for chr in "autosomes" "chrX" "chrY"
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        cat nonB_annotation/'${trivial}'_alt/'${chr}'_*.bed |sort -k1,1 -k2,2n \
        |mergeBed -i - >nonB_annotation/'${trivial}'_alt/'${chr}'_all.bed
        ' | sbatch -J ${chr} --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G
    done
done 

# Fasta index files (.fai) are placed in ref/



######################### SPACER LENGTH DISTRIBUTION ###########################
# Get all spacers from inverted repeats and mirror repeats to produce figure S1
# Spacers taken from original GFA output (tsv files), only primary assemblies

mkdir spacer
for n in "IR" "MR"
do
  cut -f10 non-B-DNA-Annotations/Homo_sapiens/chr*_$n.tsv |grep -v "Spacer" >spacer/$n.spacers.txt
done

# Plot distributions with T2T_primate_nonB/R/plot_fig1_spacer.R



############################## NON-B STATISTICS ################################

# First make a file with the totals (use primary assembly, skip random chr)
rm -f T2T_primate_nonB/helpfiles/condensed_lengths.txt
cat T2T_primate_nonB/helpfiles/primary_species_list.txt | while read -r trivial latin filename;
do
  echo "Check lengths for $trivial"
  a=`grep -v "random" ref/$filename.fai |grep -v "chrX" |grep -v "chrY" |cut -f2 |awk '{sum+=$1}END{print sum}'`
  x=`grep "chrX" ref/$filename.fai |cut -f2`
  y=`grep "chrY" ref/$filename.fai |cut -f2`
  echo "$trivial autosomes $a" >>T2T_primate_nonB/helpfiles/condensed_lengths.txt
  echo "$trivial chrX $x" >>T2T_primate_nonB/helpfiles/condensed_lengths.txt
  echo "$trivial chrY $y" >>T2T_primate_nonB/helpfiles/condensed_lengths.txt
done

# Then sum up the non-B DNA and print numbers for table 2
for chr in "autosomes" "chrX" "chrY"
do
  echo "================= $chr =================="
  echo "species ====APR=== ===DR=== ===GQ=== ===IR=== ===MR=== ===STR=== ===Z-DNA=== === all ==="
  cat T2T_primate_nonB/helpfiles/primary_species_list.txt | while read -r trivial latin filename;
  do
    out=""
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
    do
      tot=`grep $trivial T2T_primate_nonB/helpfiles/condensed_lengths.txt |grep $chr |cut -f3 -d" "`
      nonB=`awk -v t=$tot '{sum+=$3-$2}END{frac=sum/t; mb=sum/1000000; print mb" "frac}' \
      nonB_annotation/${trivial}_pri/${chr}_$nb.bed`
      tmp=$out" "$nonB
      out=$tmp
    done
    echo $trivial" "$out
  done
done

# Detailed per chromosome statistics for each species
cat T2T_primate_nonB/helpfiles/primary_species_list.txt |grep "chimp" | while read -r trivial latin filename;
do
  echo "================== $trivial ========================"
  echo "Chrom ====APR=== ===DR=== ===GQ=== ===IR=== ===MR=== ===STR=== ===Z-DNA=== === all ==="
    cut -f1,2 ref/$filename.fai | while read -r chr tot;
    do
      out=""
      for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
      do
        nonB=`cat nonB_annotation/${trivial}_pri/*_${nb}.bed \
        | awk -v t=$tot -v c=$chr '($1==c){sum+=$3-$2}END{frac=sum/t; mb=sum/1000000; print mb" "frac}'`
        tmp=$out" "$nonB
        out=$tmp
      done
    echo $chr" "$out
  done
done


######################### OVERLAP BETWEEN NON-B TYPES ##########################
# Check overlap between all non-B types and visulize with an upset plot


# For the upset plot, one need a list with all bases and non-B types  
# ~~~~~~~~~ THIS PAST WAS NOT USED AND CAN BE REMOVED
mkdir overlap
for sp in "human"
do 
    echo '#!/bin/bash      
    echo "Chr Pos Type" |sed "s/ /\t/g" >overlap/primary_'$sp'.txt
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
    do
        cat nonB_annotation/'${sp}'_pri/*_$nb.bed \
        |awk -v OFS="\t" -v n="$nb" '"'"'{for(i=$2; i<$3; ++i){print $1,i,n}}'"'"' >>overlap/primary_'$sp'.txt
    done 
    ' | sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G
done 
# This creates a 14Gb file which is quite hard to work with (takes ages to load into R)
# ~~~~~~~~~~


# Get the overlap with a python script (still takes some time and memory 
# (67Gb, 42min for human autosomes))
# For all species except bonobo & gorilla, I used 2 cores with 40Gb ram
# For bonobo & gorilla I used 3 cores with 40Gb ram
for sp in "gorilla"  #"bonobo" #"chimp" "sorang" "borang" "siamang"
do 
  # Needed more memory for autosomes than X and Y
  echo '#!/bin/bash      
module load python/3.11.2
python3 T2T_primate_nonB/python/upset_summary.py -b nonB_annotation/'$sp'_pri/autosomes_ -o overlap/'$sp'.summary.autosomes.txt
 ' | sbatch -J $sp.A --ntasks=1 --cpus-per-task=1 --mem-per-cpu=120G

  echo '#!/bin/bash      
  module load python/3.11.2
  python3 T2T_primate_nonB/python/upset_summary.py -b nonB_annotation/'$sp'_pri/chrY_ -o overlap/'$sp'.summary.chrY.txt
  ' | sbatch -J $sp.Y --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G

  echo '#!/bin/bash      
  module load python/3.11.2
  python3 T2T_primate_nonB/python/upset_summary.py -b nonB_annotation/'$sp'_pri/chrX_ -o overlap/'$sp'.summary.chrX.txt
  ' | sbatch -J $sp.X --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G
done

# THIS IS NOT DONE YET!!!!
# Autosomes was too big to run as one file for gorilla and bonobo 
# do per chromosome 
for sp in "gorilla" "bonobo" #"chimp" "sorang" "borang" "siamang"
do 
  for i in {1..23}
  do 
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
    do
      grep "chr"$i"_" nonB_annotation/${sp}_pri/autosomes_$nb.bed >nonB_annotation/${sp}_pri/chr${i}_$nb.bed
      echo '#!/bin/bash      
      module load python/3.11.2
      python3 T2T_primate_nonB/python/upset_summary.py -b nonB_annotation/'$sp'_pri/chr'$i'_ -o overlap/'$sp'.summary.chr'$i'.txt
      ' | sbatch -J $sp.$i --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
    done
  done 
done 

# Merge autosomes 

# Plot with
scripts/R/plot_fig2_upset.R



################### ENRICHMENT OF NON-B DNA IN NEW SEQUENCE ####################
# This code was also used for the pig primates autosome paper (Yoo et al 2024) 
# Code adapted from Bob Harris (from Makova et al. Nature 2024)

# Old and new assemblies divided into separate files for each chromosomes,
# named after the new chromosome number, example: "chimp/old.chr2.fa" 

# Run meryl and winnowmap 
for sp in "bonobo" "gorilla" "chimp" "sorang" "human"
do
 for i in {1..23} "X" "Y" 
 do
   job=$(echo '#!/bin/bash
   ~/software/meryl-1.4.1/bin/meryl count k=19 output '$sp'/old.chr'$i'.meryldb '$sp'/old.chr'$i'.fa
   ~/software/meryl-1.4.1/bin/meryl print greater-than distinct=0.9998 '$sp'/old.chr'$i'.meryldb > '$sp'/old.chr'$i'.repeats
   ' | sbatch -J meryl.$sp.chr$i --ntasks=1 --cpus-per-task=1 --time=05:00 --partition=open |cut -f4 -d" ")
     job2=$(echo '#!/bin/bash
   ~/software/Winnowmap/bin/winnowmap -x asm20 -c --eqx -t 4 -W '$sp'/old.chr'$i'.repeats '$sp'/old.chr'$i'.fa '$sp'/t2t.chr'$i'.fa >'$sp'/chr'$i'.winnowmap.paf
   ' | sbatch -J winnowmap.$sp.chr$i --ntasks=1 --cpus-per-task=4 --time=5:00:00 --partition=open -d afterany:$job|cut -f4 -d" ")
 done
done

# I AM HERE 


# Fasta files were taken from Edmundo (Copied from Brubeck) and Bob Harris' code
# for winnomap was used for the alignment.
# Noted that there were very little unaligned for human chr X, but way too much
# for chr Y (three times more than the length difference between the assemblies)
# so I should discuss this with Kateryna.

# Count 'new' basepairs:
cd /storage/group/kdm16/default/lbs5874/new_vs_old_alignments
cat human/merged_unaliged.bed |awk '{sum+=$3-$2}END{print sum}'
# 88569799
cat human/chrX.unaligned.bed |awk '{sum+=$3-$2}END{print sum}'
# 136783
cat human/chrY.unaligned.bed |awk '{sum+=$3-$2}END{print sum}'
# 15097642

# Make a table that can be used for statistic tests
# Take the content from the summary files on Roar and place them in
ls new_sequence/summary_human_*.txt
# Then extract the numbers for each type
echo "Chr nonB Region bp_outside_nonB bp_inside_nonB" | sed 's/ /\t/' >new_sequence/human_stat_table.tsv
for chr in "autosomes" "chrX" "chrY"
do
   awk -v c=$chr '(NR>1){out_new=$2-$3; out_old=$5-$6; print c,$1,"New",out_new,$3,"\n",c,$1,"Old",out_old,$6}' new_sequence/summary_human_${chr}.txt |sed 's/^ //' |sed 's/ /\t/g' >>new_sequence/human_stat_table.tsv
done

# With fractions
echo "Chr nonB Region bp_outside_nonB bp_inside_nonB" | sed 's/ /\t/' >new_sequence/human_stat_table_fractions.tsv
for chr in "autosomes" "chrX" "chrY"
do
   awk -v c=$chr '(NR>1){out_new=($2-$3)/$2; in_new=$3/$2; out_old=($5-$6)/$5; in_old=$6/$5; print c,$1,"New",out_new,in_new,"\n",c,$1,"Old",out_old,in_old}' new_sequence/summary_human_${chr}.txt |sed 's/^ //' |sed 's/ /\t/g' >>new_sequence/human_stat_table_fractions.tsv
done

# More species
# First make the summary files
cat species_list.txt |grep -v "siamang" |grep -v "human" |while read -r trivial latin filename;
do
  for chr in "autosomes" "chrX" "chrY"
  do
    touch new_sequence/summary_${trivial}_${chr}.txt
  done
done

# Maybe a Mann Whitney U test is better than a goodness of fit/chisquare test,
# for that I need densities instead, for ALL species
echo "Species Chr nonB Region density" | sed 's/ /\t/' >new_sequence/5sp_stat_table_densities.tsv
cat species_list.txt |grep -v "siamang" |grep -v "borang" |while read -r trivial latin filename;
do
  for chr in "autosomes" "chrX" "chrY"
  do
     awk -v c=$chr  -v sp=$trivial '(NR>1){print sp,c,$1,"New",$4,"\n",sp,c,$1,"Old",$7}' new_sequence/summary_${trivial}_${chr}.txt |sed 's/^ //' |sed 's/ /\t/g' >>new_sequence/5sp_stat_table_densities.tsv
  done
done




