################################################################################
# CODE FOR THE PAPER 
# Non-canonical DNA in human and other ape telomere-to-telomere genomes
# https://www.biorxiv.org/content/10.1101/2024.09.02.610891v1
# Written by Linnéa Smeds Aug-Sept 2024

# Non-B annotation in bb format can be found here:
# Human
# https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/browser/CHM13/bbi/nonB_*.bb
# Other apes:
# https://genomeark.s3.amazonaws.com/species/*/*/assembly_curated/repeats/*_v2.*.nonB_*.bb


# CONTENTS 

####################### NOTES AND SOFTWARE REQUIREMENTS ######################## 
# Some code is written to be run as parallel batch jobs on a slurm cluster for 
# speed and efficiency. These could be run in the terminal istead by removing 
# the initial "echo '#!/bin/bash" and the trailing  "| sbatch [..] "

# Software used:
# bedtools v2.31.0
# samtools v1.19.2
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
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep -v "chimp" | while read -r sp latin filename;
do
    echo "Creating directory for $sp.."
    mkdir -p nonB_annotation/${sp}_pri
    prefix=`echo $filename | sed 's/.fasta//' | sed 's/.fa//'`
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
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

# Secondary haplotype assembly 
cat T2T_primate_nonB/helpfiles/alt_species_list.txt |grep siamang | while read -r sp latin filename;
do
    echo "Creating directory for $sp.."
    mkdir -p nonB_annotation/${sp}_alt
     prefix=`echo $filename | sed 's/.fasta//' | sed 's/.fa//'`
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
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
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
    do
      cat nonB_annotation/${sp}_$hap/autosomes_$nb.bed nonB_annotation/${sp}_$hap/chrX_$nb.bed nonB_annotation/${sp}_$hap/chrY_$nb.bed >nonB_annotation/${sp}_$hap/genome_$nb.bed
    done 
  done
done


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

# Then sum up the non-B DNA and print numbers for table 2
for chr in "autosomes" "chrX" "chrY"
do
  echo "================= $chr =================="
  echo "species ====APR=== ===DR=== ===GQ=== ===IR=== ===MR=== ===STR=== ===Z-DNA=== === all ==="
  cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
  do
    out=""
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
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
 #   for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
 #   do
      #grep "chr"$i"_" nonB_annotation/${sp}_pri/autosomes_$nb.bed >nonB_annotation/${sp}_pri/chr${i}_$nb.bed
  #  done
    echo '#!/bin/bash      
    module load python/3.11.2
    python3 T2T_primate_nonB/python/upset_summary.py -b nonB_annotation/'$sp'_pri/chr'$i'_ -o overlap/'$sp'.summary.chr'$i'.txt
    ' | sbatch -J $sp.$i --ntasks=1 --mem-per-cpu=8G --cpus-per-task=1 --time=10:00:00 --partition=open
  done 
done 

# Merge autosomes 


# Plot with
scripts/R/plot_fig2_upset.R


################################# COMPARE WITH SASWAT 

cd /storage/group/kdm16/default/skm6640/work/g4set_fork/Homo_sapiens/

cat nonB_annotation/human_pri/*_GQ.bed |awk '{sum+=$3-$2; n++}END{print n, sum}'
739112 27433138

cat /storage/group/kdm16/default/skm6640/work/g4set_fork/Homo_sapiens/chrG.pqsfinder.filtered.bed |awk '{sum+=$3-$2; n++}END{print n, sum}'
769188 19899699

 cat nonB_annotation/human_pri/*_GQ.bed | intersectBed -a - -b <(cat /storage/group/kdm16/default/skm6640/work/g4set_fork/Homo_sapiens/chrG.pqsfinder.filtered.bed) |awk '{sum+=$3-$2; n++}END{print n, sum}'
568692 15505883