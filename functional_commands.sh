################## ENRICHMENT OF NON-B IN FUNCTIONAL ANNOTATION #################
# written by Linnéa Smeds, October 2024
# Code for analysing non-B DNA enrichment in funtional regions of the genome. 


# Regions annotated by Sasway Mohanty, taken from the preprint 
# https://www.biorxiv.org/content/10.1101/2024.11.05.621973v1.full
# Code used for generating the annotation files are found in 
# https://github.com/makovalab-psu/GreatApeT2T-G4s/tree/main/datasets/dataForAnalysis/geneAnnotations

# Abbreviations for functional classes and file names used inhouse our found in 
# T2T_primate_nonB/helpfiles/functional_classes.txt 
func_path=/path/to/annotation/files
mkdir ref/annotation/human/
for file in $(ls $func_path)
do 
    ln -s $func_path/$file ref/annotation/human/$file
    awk '($1!="chrX" && $1!="chrY"){print}' $func_path/$file |grep -v "chrY" >ref/annotation/human/autosomes.$file
    awk '($1=="chrX"){print}' $func_path/$file >ref/annotation/human/chrX.$file
    awk '($1=="chrY"){print}' $func_path/$file >ref/annotation/human/chrY.$file
done 


# INTERSECT FUNCTIONAL ANNOTATION WITH NON-B
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Region Class APR DR GQ IR MR STR Z" |sed "s/ /\t/g" >functional/'${sp}'/enrichment_fullgenome.tsv
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        echo "Looking at $class" 
        tmp="genome "$class
        len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' ref/annotation/'$sp'/$bedfile`
        cat densities/'${sp}'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
        do
            echo "looking at $non_b"
            d=`intersectBed -a <(cut -f1-4 ref/annotation/human/$bedfile) -b <(cat nonB_annotation/'${sp}'_pri/genome_${non_b}.bed) -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$8}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
            tmp=`echo $tmp" "$d`
            echo $tmp >tmp
        done
        cat tmp |sed "s/ /\t/g" >>functional/'${sp}'/enrichment_fullgenome_new.tsv
    done
    ' | sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --partition=open
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CORRECT FOR GC CONTENT 

# Use method from Guiblet et al 2021 (https://github.com/makovalab-psu/G4_Selection)
# Dividing by GC_cont_region/GC_cont_genome

# Calculate GC for each element class 
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        echo '#!/bin/bash
        module load bedtools/2.31.0
        newfile=`echo '$bedfile' |sed "s/.bed/.withGC.bed/"`
        bedtools nuc -fi ref/assemblies/'$filename' -bed <(cut -f1-3 ref/annotation/'$sp'/'$bedfile') >ref/annotation/'$sp'/$newfile
        '| sbatch -J $class --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
    done 
done 

### CORRECT GQ WITH A SIMPLE FACTOR 
# Divide fold enrichment by (GC_content_in_region/GC_content_genome) 
# First count AT and GC content
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Region Class GCnum GCcont ATnum ATcont" |sed "s/ /\t/g" >ref/annotation/'${sp}'/all_classes_GC_AT.txt
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        gcfile=`echo $bedfile |sed "s/.bed/.withGC.bed/"`
        gc_region=`awk '"'"'{gc+=$7+$8; at+=$6+$9; tot+=$12}END{gccont=gc/tot; atcont=at/tot; print gc,gccont,at,atcont}'"'"' ref/annotation/'$sp'/$gcfile`
        echo "genome $class $gc_region" | sed "s/ /\t/g" >>ref/annotation/'${sp}'/all_classes_GC_AT.txt
    done 
    ' | sbatch --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done 

# Then go through the enrichment file and correct the GQ value 
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    gc_genome=`awk  '{gc=$1; tot=$1+$2}END{gccont=gc/tot; print gccont}' ref/assemblies/$filename.nuc`
    paste <(cut -f1,2,5 functional/$sp/enrichment_fullgenome.tsv) ref/annotation/$sp/all_classes_GC_AT.txt | \
    awk -v gcg=$gc_genome -v OFS="\t" '{if(NR==1){print $1,$2,$3}else{new=$3/($7/gcg); print $1,$2,new}}'  >functional/$sp/enrichment_GQ_MYcorr_fullgenome.tsv
done

# Get GC for each region and plot against nonB-density
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        echo "nonB Element Chr GCcont Density" |sed 's/ /\t/g' >functional/$sp/tables/$class.GC.nonB.txt   
        cat densities/human_pri_nonB_genome_wide.txt |grep -v "all"  | while read -r non_b tot dens;
        do 
#            echo '#!/bin/bash
 #           module load bedtools/2.31.0
  #          newfile=`echo '$bedfile' |sed "s/.bed/.'$non_b'.bed/"`
      #      intersectBed -a <(cut -f1-4 ref/annotation/human/'$bedfile') -b <(cat nonB_annotation/'$sp'_pri/genome_'$non_b'.bed) -wao |awk -v OFS="\t" '"'"'{if(NR==1){chr=$1; s=$2; e=$3; n=$4; sum=$8}else{if(chr==$1 && s==$2 && e==$3){sum+=$8}else{dens=sum/(e-s); print chr,s,e,n,sum,dens; chr=$1; s=$2; e=$3; n=$4; sum=$8}}}END{dens=sum/(e-s); print chr,s,e,n,sum,dens}'"'"'  >functional/human/overlap/$newfile

   #         ' | sbatch -J $ind --ntasks=1 --cpus-per-task=1  --mem-per-cpu=8G --time=1:00:00 --partition=open 
            gcfile=`echo $bedfile |sed "s/.bed/.withGC.bed/"`
            nonbfile=`echo $bedfile |sed "s/.bed/.$non_b.bed/"`
            join -1 4 -2 4 <(sort -k4,4 ref/annotation/$sp/$gcfile) <(sort -k4,4  functional/human/overlap/$nonbfile) |cut -f1,2,8,20 -d" " |awk -v nb=$non_b -v OFS="\t" '{print nb,$0}' >>functional/$sp/tables/$class.GC.nonB.txt
         done
    done
done
# There are no clear relationships with GC for the other non-B types! 

