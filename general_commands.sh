################################################################################
# CODE FOR THE PAPER 
# Non-canonical DNA in human and other ape telomere-to-telomere genomes
# https://www.biorxiv.org/content/10.1101/2024.09.02.610891v3
# Written by Linnéa Smeds Aug 2024 - March 2025

# Non-B annotation in bb format can be found here:
# Human
# https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/browser/CHM13/bbi/nonB_*.bb
# Other apes:
# https://genomeark.s3.amazonaws.com/species/*/*/assembly_curated/repeats/*_v2.*.nonB_*.bb


# # # # # TABLE OF CONTENTS 
# - NOTES AND SOFTWARE REQUIREMENTS
# - CREATE NON-OVERLAPPING BED FILES FROM NON-B ANNOTATION
# - SPACER LENGTH DISTRIBUTION
# - NON-B STATISTICS
# - OVERLAP BETWEEN NON-B TYPES
# - TESTING I-MOTIF SEARCH WITH IM-SEEKER

####################### NOTES AND SOFTWARE REQUIREMENTS ######################## 
# Some code is written to be run as parallel batch jobs on a slurm cluster for 
# speed and efficiency. These could be run in the terminal instead by removing 
# the initial "echo '#!/bin/bash" and the trailing  "| sbatch [..] "

# Software used:
# bedtools v2.31.0
# samtools v1.19.2
# meryl v1.4.1
# Winnowmap v2.03
# Bedops v2.4.41
# circos v0.69-9
# starcode v1.4 
# python v3.11.2
# UCSC tools:
#   bigWigToWig
#   bigBedToBed
#   liftOver

############ CREATE NON-OVERLAPPING BED FILES FROM NON-B ANNOTATION ############
nb_path=/path/to/nonB/tracks/converted/to/bed
mkdir nonB_annotation 

# Primary haplotype assembly 
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep human | while read -r sp latin filename;
do
    echo "Creating directory for $sp.."
 #   mkdir -p nonB_annotation/${sp}_pri
    prefix=`echo $filename | sed 's/.fasta//' | sed 's/.fa//'`
    for nb in "TRI" #"APR" "DR" "GQ" "IR" "MR" "STR" "Z" 
    do
        #echo $prefix
        #ls $nb_path/$prefix.$nb.bed
        grep "chrX" $nb_path/$prefix.$nb.bed |mergeBed -i - >nonB_annotation/${sp}_pri/chrX_$nb.bed
        grep "chrY" $nb_path/$prefix.$nb.bed |mergeBed -i - >nonB_annotation/${sp}_pri/chrY_$nb.bed
        grep -v "chrX" $nb_path/$prefix.$nb.bed |grep -v "chrY" \
        |mergeBed -i - >nonB_annotation/${sp}_pri/autosomes_$nb.bed
    done 
   for chr in "autosomes" "chrX" "chrY"
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        cat nonB_annotation/'${sp}'_pri/'${chr}'_*.bed |sort -k1,1 -k2,2n \
        |mergeBed -i - >nonB_annotation/'${sp}'_pri/'${chr}'_all.bed
        ' | sbatch -J ${chr} --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G
    done
done 

nb_path="/storage/group/kdm16/default/shared/nonB_annotation/T2Tv2_primates/"
# Secondary haplotype assembly 
cat T2T_primate_nonB/helpfiles/alt_species_list.txt | grep siamang | while read -r sp latin filename;
do
    echo "Creating directory for $sp.."
    mkdir -p nonB_annotation/${sp}_alt
     prefix=`echo $filename | sed 's/.fasta//' | sed 's/.fa//'`
    for nb in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z"
    do
        grep "chrX" $nb_path/$prefix.$nb.bed |mergeBed -i - >nonB_annotation/${sp}_alt/chrX_$nb.bed
        grep "chrY" $nb_path/$prefix.$nb.bed |mergeBed -i - >nonB_annotation/${sp}_alt/chrY_$nb.bed
        grep -v "chrX" $nb_path/$prefix.$nb.bed |grep -v "chrY" \
        |mergeBed -i - >nonB_annotation/${sp}_alt/autosomes_$nb.bed
    done
    for chr in "autosomes" "chrX" "chrY"
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        cat nonB_annotation/'${sp}'_alt/'${chr}'_*.bed |sort -k1,1 -k2,2n \
        |mergeBed -i - >nonB_annotation/'${sp}'_alt/'${chr}'_all.bed
        ' | sbatch -J ${chr} --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G
    done
done 

# Fasta index files (.fai) are placed in ref/

# Merged files to simplify the code for some analyses 
for hap in "pri"  "alt" 
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |while read -r sp latin filename;
  do
    for nb in "TRI" #"all" "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
    do
      cat nonB_annotation/${sp}_$hap/autosomes_$nb.bed nonB_annotation/${sp}_$hap/chrX_$nb.bed nonB_annotation/${sp}_$hap/chrY_$nb.bed >nonB_annotation/${sp}_$hap/genome_$nb.bed
    done 
  done
done


######################### SPACER LENGTH DISTRIBUTION ###########################
# Get all spacers from inverted repeats and mirror repeats to produce figure S1
# Spacers taken from original GFA output (tsv files), only primary assemblies

mkdir spacer
for n in "DR" "IR" "MR"
do
  cut -f10 non-B-DNA-Annotations/Homo_sapiens/chr*_$n.tsv |grep -v "Spacer" >spacer/$n.spacers.txt
done

# Plot distributions with R/plot_figS1_spacer.R



############################## NON-B STATISTICS ################################

# First make a file with the totals (use primary assembly, skip random chr)
rm -f T2T_primate_nonB/helpfiles/condensed_lengths.txt
cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
do
  echo "Check lengths for $sp"
  a=`grep -v "random" ref/$filename.fai |grep -v "chrX" |grep -v "chrY" |cut -f2 |awk '{sum+=$1}END{print sum}'`
  x=`grep "chrX" ref/$filename.fai |cut -f2`
  y=`grep "chrY" ref/$filename.fai |cut -f2`
  echo "$sp autosomes $a" >>T2T_primate_nonB/helpfiles/condensed_lengths.txt
  echo "$sp chrX $x" >>T2T_primate_nonB/helpfiles/condensed_lengths.txt
  echo "$sp chrY $y" >>T2T_primate_nonB/helpfiles/condensed_lengths.txt
done

# Then sum up the non-B DNA and print numbers for table S1
for chr in "autosomes" "chrX" "chrY"
do
  echo "================= $chr =================="
  echo "species ====APR=== ===DR=== ===GQ=== ===IR=== ===MR=== ===TRI=== ===STR=== ===Z-DNA=== === all ==="
  cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
  do
    out=""
    for nb in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z" "all"
    do
      tot=`grep $sp T2T_primate_nonB/helpfiles/condensed_lengths.txt |grep $chr |cut -f3 -d" "`
      nonB=`awk -v t=$tot '{sum+=$3-$2}END{frac=sum/t; mb=sum/1000000; print mb" "frac}' \
      nonB_annotation/${sp}_pri/${chr}_$nb.bed`
      tmp=$out" "$nonB
      out=$tmp
    done
    echo $sp" "$out
  done
done

# SAME BUT PRINT TO FILE (added TRI, subset of MR) 
echo "Species Chr Type APR DR GQ IR MR TRI STR Z all" |sed 's/ /\t/g' >nonB_annotation/7sp_summary.txt
for chr in "autosomes" "chrX" "chrY"
do
  echo "================= $chr =================="
  echo "species ====APR=== ===DR=== ===GQ=== ===IR=== ===MR=== ===TRI=== ===STR=== ===Z-DNA=== ===all==="
  cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
  do
    out1=""
    out2=""
    for nb in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z" "all"
    do
      tot=`grep $sp T2T_primate_nonB/helpfiles/condensed_lengths.txt |grep $chr |cut -f3 -d" "`
      nonB=`awk -v t=$tot '{sum+=$3-$2}END{frac=sum/t; mb=sum/1000000; print mb" "frac}' \
      nonB_annotation/${sp}_pri/${chr}_$nb.bed`
      nonB_Mb=`echo $nonB |cut -f1 -d" "`
      nonB_frac=`echo $nonB |cut -f2 -d" "`
      tmp1=$out1" "$nonB_Mb
      tmp2=$out2" "$nonB_frac
      out1=$tmp1
      out2=$tmp2
    done
    echo $sp" "$chr" Mb "$out1 |sed 's/ /\t/g' >>nonB_annotation/7sp_summary.txt
    echo $sp" "$chr" frac "$out2 |sed 's/ /\t/g' >>nonB_annotation/7sp_summary.txt
  done
done

# TABLE S1B, ONLY TRI 
echo "Species Chr Type TRI" |sed 's/ /\t/g' >nonB_annotation/7sp_TRI_summary.txt
for chr in "autosomes" "chrX" "chrY"
do
  echo "================= $chr =================="
  echo "species ====TRI===="
  cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
  do
    out1=""
    out2=""
    for nb in "TRI"
    do
      tot=`grep $sp T2T_primate_nonB/helpfiles/condensed_lengths.txt |grep $chr |cut -f3 -d" "`
      nonB=`awk -v t=$tot '{sum+=$3-$2}END{frac=sum/t; mb=sum/1000000; print mb" "frac}' \
      nonB_annotation/${sp}_pri/${chr}_$nb.bed`
      nonB_Mb=`echo $nonB |cut -f1 -d" "`
      nonB_frac=`echo $nonB |cut -f2 -d" "`
      tmp1=$out1" "$nonB_Mb
      tmp2=$out2" "$nonB_frac
      out1=$tmp1
      out2=$tmp2
    done
    echo $sp" "$chr" Mb "$out1 |sed 's/ /\t/g' >>nonB_annotation/7sp_TRI_summary.txt
    echo $sp" "$chr" frac "$out2 |sed 's/ /\t/g' >>nonB_annotation/7sp_TRI_summary.txt
  done
done


# NOT ADDED TO THE PAPER 
# Detailed per chromosome statistics for each species
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "chimp" | while read -r sp latin filename;
do
  echo "================== $sp ========================"
  echo "Chrom ====APR=== ===DR=== ===GQ=== ===IR=== ===MR=== ===STR=== ===Z-DNA=== === all ==="
    cut -f1,2 ref/$filename.fai | while read -r chr tot;
    do
      out=""
      for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
      do
        nonB=`cat nonB_annotation/${sp}_pri/*_${nb}.bed \
        | awk -v t=$tot -v c=$chr '($1==c){sum+=$3-$2}END{frac=sum/t; mb=sum/1000000; print mb" "frac}'`
        tmp=$out" "$nonB
        out=$tmp
      done
    echo $chr" "$out
  done
done


# TABLE WITH TOTAL MOTIF NUMBER AND MEDIAN AND LARGEST MOTIF FOR EACH SPECIES 
# (Data is shown for the primary haplotypes)
echo "Species Type APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g' >nonB_annotation/7sp_length_and_number.txt
cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
do
  prefix=`echo $filename | sed 's/.fasta//' | sed 's/.fa//'`
  echo $prefix 
  out1="$sp Total"
  out2="$sp Median"
  out3="$sp Largest"
  for nb in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z"
  do
    stats=`python3 T2T_primate_nonB/python/stats_from_bed.py -b ../../shared/nonB_annotation/T2Tv2_primates/$prefix.$nb.bed -p all`
    tot=`echo $stats |cut -f1 -d" "`
    median=`echo $stats |cut -f2 -d" "` 
    largest=`echo $stats |cut -f5 -d" "`
    tmp1=$out1" "$tot
    tmp2=$out2" "$median
    tmp3=$out3" "$largest
    out1=$tmp1
    out2=$tmp2
    out3=$tmp3
  done 
  echo $out1 |sed 's/ /\t/g' >>nonB_annotation/7sp_length_and_number.txt
  echo $out2 |sed 's/ /\t/g' >>nonB_annotation/7sp_length_and_number.txt
  echo $out3 |sed 's/ /\t/g' >>nonB_annotation/7sp_length_and_number.txt
done 



######################### OVERLAP BETWEEN NON-B TYPES ##########################
# Check overlap between all non-B types and visulize with an upset plot
# This analysis is run on the primary haplotype ONLY

# Get the overlap with a python script (takes some time and memory -
# (67Gb, 42min for human autosomes))
# For all species except bonobo & gorilla, I used 2 cores with 40Gb ram
# For bonobo & gorilla I used 3 cores with 40Gb ram
for sp in "chimp" "sorang" "borang" "siamang" #"human"
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


# If there isn't enough memory to run all autosomes together, each chromosome
# can be run separately and merged afterwards.
for sp in "gorilla" #"bonobo"
do 
  for i in 18 #{1..23}
  do 
    #for nb in "TRI" #"APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z"
    #do
    #  grep "chr"$i"_" nonB_annotation/${sp}_pri/autosomes_$nb.bed >nonB_annotation/${sp}_pri/chr${i}_$nb.bed
    #done
    echo '#!/bin/bash
    echo "Running upset for '$sp' chr'$i'"      
    module load python/3.11.2
    python3 T2T_primate_nonB/python/upset_summary.py -b nonB_annotation/'$sp'_pri/chr'$i'_ -o overlap/'$sp'.summary.chr'$i'.txt
    ' | sbatch -J $sp.$i --ntasks=1 --mem-per-cpu=20G --cpus-per-task=1 --time=10:00:00 --partition=open
  done 
done 
# Merge autosomes 
for sp in "gorilla" #"bonobo"
do 
  cat overlap/$sp.summary.chr{1..23}.txt |sort | awk -v OFS="\t" '{if(NR==1){type=$1; sum=$2}else{if($1==type){sum+=$2}else{print type,sum; type=$1; sum=$2}}}END{print type,sum}' > overlap/$sp.summary.autosomes.txt
done 

# PAIRWISE OVERLAP WITH FRACTION
# Note that we cannot just take the pairwise rows (APR-DR etc) from the files
# created above, since threeway or more overlap also contributes to the pairwise
# fraction. 
module load bedtools/2.31.0 
for sp in "gorilla" "bonobo" "chimp" "sorang" "borang" "siamang" "human"
do 
  echo '#!/bin/bash 
  module load bedtools/2.31.0 
  echo "Chr NB1 NB2 Ovl Frac" | sed "s/ /\t/g" >overlap/'$sp'.pairwise.summary.txt
  for chr in "autosomes" "chrX" "chrY"
  do 
    # Get the total number of bp for each type (save for later in a summary file, can be skipped if files are already generated)
    echo "#NonB bp" |sed "s/ /\t/g" >nonB_annotation/'${sp}'_pri/annotated_bp_${chr}.txt
    for nb in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z"
    do 
      awk -v nb=$nb '"'"'{sum+=$3-$2}END{print nb"\t"sum}'"'"' nonB_annotation/'${sp}'_pri/${chr}_$nb.bed >>nonB_annotation/'${sp}'_pri/annotated_bp_${chr}.txt
    done 
  # Go through all combinations 
    for nb in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z"
    do 
      nbbp=`grep $nb nonB_annotation/'${sp}'_pri/annotated_bp_${chr}.txt |cut -f2`
      echo "$nb has $nbbp number of bp!" 
      for nb2 in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z"
      do
          if [[ $nb != $nb2 ]]
          then 
            intersectBed -a nonB_annotation/'${sp}'_pri/${chr}_$nb.bed -b nonB_annotation/'${sp}'_pri/${chr}_$nb2.bed -wo | awk -v tot=$nbbp -v chr=$chr -v nb1=$nb -v nb2=$nb2 -v OFS="\t" '"'"'{sum+=$7}END{f=sum/tot; print chr,nb1,nb2,sum,f}'"'"' >>overlap/'$sp'.pairwise.summary.txt
          fi 
      done
    done 
  done 
   ' | sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=1:00:00 
done 

# Plot with R/plot_fig3_overlap.R and R/plot_figS2_overlap.R


##################### TESTING I-MOTIF SEARCH WITH IM-SEEKER ####################
# Test to see if there is a point searching for i-motifs in addition to G4.
# IM seeker is a web based tools so I started by testing it on human chr1 
# (uploaded fasta to https://im-seeker.org/predict and downloaded the results)

mkdir test_imotif
# Convert the csv to bed format
cut -f2,3,4 -d"," ~/Downloads/f32f9f55-ddd1-42fb-b7b8-583f5fbe7b72_imotif_prediction.csv |sed 's/,/\t/g' |tail -n+2 |awk -v OFS="\t" '{$2=$2-1; print}' >test_imotif/Human_chr1.IM.bed

# Number of i-motifs found:
wc test_imotif/Human_chr1.IM.bed
#   35993  107979  847167 test_imotif/Human_chr1.IM.bed

# Number of i-motifs that overlap with G4 annotations
intersectBed -a test_imotif/Human_chr1.IM.bed -b nonB_annotation/human_pri/genome_GQ.bed |uniq |wc
#   35909  107727  846203
# 35909/35993=0.998 => Well over 99%! It doesn't make sense to add i-motifs as 
# a separate category in the paper.

 

