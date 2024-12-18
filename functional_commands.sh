################## ENRICHMENT OF NON-B IN FUNCTIONAL ANNOTATION #################
# written by Linnéa Smeds, October 2024
# Code for analysing non-B DNA enrichment in funtional regions of the genome. 

# INHOUSE to separate annotation files into regions 
mkdir ref/annotation
# Use Saswat's annotation files, make list of them and 
# add category manually 
echo "CDS_ProtCode cds_protCode_singleTrans.bed
CpG_Islands chrG_cpgislands.bed
Origin_of_Repl core_stochastic_origins_chm13v2.0_hglft.ucsc.bed
Enhancers enhancer_protCode_singleTrans.bed
Exons_ProtCode exon_protCode_singleTrans.bed
Gene_NonProtCode gene_nonprotCode.bed
Gene_ProtCode gene_protCode_singleTrans.bed
Introns intron_protCode_singleTrans.bed
NGNR ngnr_protCode_singleTrans.bed
Promoters promoter_protCode_singleTrans.bed
Repeats repeats_protCode_singleTrans.bed
RNA_ProtCode rna_protCode_singleTrans.bed
3UTR utr3_protCode_singleTrans.bed
5UTR utr5_protCode_singleTrans.bed"|sed 's/ /\t/g' >T2T_primate_nonB/helpfiles/functional_classes.txt
#
#ls /storage/group/kdm16/default/shared/T2Tv2.functionalAnnotations/human/ |awk '{split($1,s,"_"); print s[1]"\t"$0}' >ref/annotation/human/classes.txt
mkdir ref/annotation/human/
for file in $(ls /storage/group/kdm16/default/shared/T2Tv2.functionalAnnotations/human/)
do 
    ln -s /storage/group/kdm16/default/shared/T2Tv2.functionalAnnotations/human/$file ref/annotation/human/$file
    awk '($1!="chrX" && $1!="chrY"){print}' /storage/group/kdm16/default/shared/T2Tv2.functionalAnnotations/human/$file |grep -v "chrY" >ref/annotation/human/autosomes.$file
    awk '($1=="chrX"){print}' /storage/group/kdm16/default/shared/T2Tv2.functionalAnnotations/human/$file >ref/annotation/human/chrX.$file
    awk '($1=="chrY"){print}' /storage/group/kdm16/default/shared/T2Tv2.functionalAnnotations/human/$file >ref/annotation/human/chrY.$file
done 


# NOT USED FOR THE PAPER!!!! 
# Go through the functional classes and intersect with non-B, one chromosome 
# category at the time 
module load bedtools/2.31.0
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    mkdir -p functional/${sp}/

     echo "Region Class APR DR GQ IR MR STR Z" |sed 's/ /\t/g' >functional/${sp}/enrichment_region.tsv
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        for chr in "autosomes" "chrX" "chrY" 
        do
            tmp=$chr" "$class
            len=`awk '{sum+=$3-$2}END{print sum}' ref/annotation/$sp/$chr.$bedfile`
            echo "$chr $class $len"
            cat densities/${sp}_pri_nonB_$chr.txt |grep -v "all" | while read -r non_b tot dens;
            do
                d=`intersectBed -a <(cut -f1-4 ref/annotation/human/$chr.$bedfile) -b nonB_annotation/${sp}_pri/${chr}_${non_b}.bed -wo |awk -v l=$len -v dtot=$dens '{sum+=$8}END{d=sum/l; frac=d/dtot; print frac}'`
                tmp=`echo $tmp" "$d`
                echo $tmp >tmp.$sp.$class
            done
        cat tmp.$sp.$class |sed 's/ /\t/g'
        cat tmp.$sp.$class |sed 's/ /\t/g' >>functional/${sp}/enrichment_region.tsv
        done
    done
done



# And for full genome
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

#Not used
#awk '(!/^>/){gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"")}END{print gc"\t"at}' >$nucfile


#Calculate GC corrected enrichment for GQ
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Region Class GQ" |sed "s/ /\t/g" >functional/'${sp}'/enrichment_GQ_corr_fullgenome.tsv
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        echo "Looking at $class" 
        tmp="genome "$class
        gcfile=`echo $bedfile |sed "s/.bed/.withGC.bed/"`
        gc=`awk '"'"'{sum+=$7+$8}END{print sum}'"'"' ref/annotation/'$sp'/$gcfile`
        cat densities/'${sp}'_pri_GCcorr_nonB_genome_wide.txt |grep "GQ" | while read -r non_b tot dens;
        do
            echo "looking at $non_b"
            d=`intersectBed -a <(cut -f1-4 ref/annotation/human/$bedfile) -b <(cat nonB_annotation/'${sp}'_pri/genome_${non_b}.bed) -wo |awk -v l=$gc -v dtot=$dens '"'"'{sum+=$8}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
            tmp=`echo $tmp" "$d`
            echo $tmp >tmp
        done
        cat tmp |sed "s/ /\t/g" >>functional/'${sp}'/enrichment_GQ_corr_fullgenome.tsv
    done
    ' | sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done


#Calculate AT corrected enrichment for APR
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Region Class APR" |sed "s/ /\t/g" >functional/'${sp}'/enrichment_APR_corr_fullgenome.tsv
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        echo "Looking at $class" 
        tmp="genome "$class
        gcfile=`echo $bedfile |sed "s/.bed/.withGC.bed/"`
        at=`awk '"'"'{sum+=$6+$9}END{print sum}'"'"' ref/annotation/'$sp'/$gcfile`
        cat densities/'${sp}'_pri_ATcorr_nonB_genome_wide.txt |grep "APR" | while read -r non_b tot dens;
        do
            echo "looking at $non_b"
            d=`intersectBed -a <(cut -f1-4 ref/annotation/human/$bedfile) -b <(cat nonB_annotation/'${sp}'_pri/genome_${non_b}.bed) -wo |awk -v l=$at -v dtot=$dens '"'"'{sum+=$8}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
            tmp=`echo $tmp" "$d`
            echo $tmp >tmp
        done
        cat tmp |sed "s/ /\t/g" >>functional/'${sp}'/enrichment_APR_corr_fullgenome.tsv
    done
    ' | sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done

# Need to add empty columns for the uncorrected types 
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    echo 
    module load bedtools/2.31.0
    echo "Region Class DR IR MR STR Z" |sed "s/ /\t/g" >functional/${sp}/enrichment_un_corr_fullgenome.tsv
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        echo "genome "$class" 0 0 0 0 0" >>functional/${sp}/enrichment_un_corr_fullgenome.tsv
    done 
done 

### NEW METHOD, JUST DIVIDE FOLD-ENRICHMENT ON (GC_CONT_REGION / GC_CONT_GENOME)
# FIRST GET GC-CONT and AT-CONT
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
# AND APR
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    at_genome=`awk  '{at=$2; tot=$1+$2}END{atcont=at/tot; print atcont}' ref/assemblies/$filename.nuc`
    paste <(cut -f1,2,3 functional/$sp/enrichment_fullgenome.tsv) ref/annotation/$sp/all_classes_GC_AT.txt | \
    awk -v atg=$at_genome -v OFS="\t" '{if(NR==1){print $1,$2,$3}else{new=$3/($9/atg); print $1,$2,new}}'  >functional/$sp/enrichment_APR_MYcorr_fullgenome.tsv
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



 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STATISTICAL TESTS  

# PERMUTATION 
mkdir functional/human/permutation

# Shuffle the regions (outside of the original regions), calculate GC and 
# density (for now, not GC corrected)
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    for i in {1..1}
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        cat T2T_primate_nonB/helpfiles/functional_classes.txt |grep 'Repeats' |while read -r class bedfile
        do 
       #     bedtools shuffle -i ref/annotation/human/$bedfile -excl ref/annotation/human/$bedfile -g <(cut -f1,2 ref/'$filename'.fai |grep -v chrM) >functional/'$sp'/permutation/$class.perm_'$i'.bed 
      #      bedtools nuc -fi ref/assemblies/'$filename' -bed <(cut -f1-3 functional/'$sp'/permutation/$class.perm_'$i'.bed) >functional/'$sp'/permutation/$class.perm_'$i'.withGC.bed
            echo "Looking at $class" 
            tmp="genome "$class
            len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' ref/annotation/human/$bedfile`
            cat densities/'${sp}'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
            do
                echo "looking at $non_b"
                d=`intersectBed -a <(cut -f1-4 functional/'$sp'/permutation/$class.perm_'$i'.bed) -b <(cat nonB_annotation/'${sp}'_pri/genome_${non_b}.bed) -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$8}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
                tmp=`echo $tmp" "$d`
                echo $tmp |sed "s/ /\t/g" >functional/'$sp'/permutation/enrichment.$class.perm_'$i'.tsv
            done
        done
        ' | sbatch -J $i --ntasks=1 --cpus-per-task=2  --mem-per-cpu=8G --time=10:00:00 --partition=open 
    done
done





