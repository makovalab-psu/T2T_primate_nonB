################## ENRICHMENT OF NON-B IN FUNCTIONAL ANNOTATION #################
# written by Linnéa Smeds, October 2024
# Code for analysing non-B DNA enrichment in funtional regions of the genome. 


# Regions annotated by Sasway Mohanty, taken from the preprint 
# https://www.biorxiv.org/content/10.1101/2024.11.05.621973v1.full
# Code used for generating the annotation files are found in 
# https://github.com/makovalab-psu/GreatApeT2T-G4s/tree/main/datasets/dataForAnalysis/geneAnnotations
# And the actual files can be downloaded from: 
# https://github.com/makovalab-psu/GreatApeT2T-G4s/tree/main/datasets/functionalOutputs/Homo_sapiens

# Abbreviations for functional classes and file names used inhouse are found in 
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

# Check size of annotated regions 
cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
do 
    len1=`awk '{sum+=$3-$2}END{print sum}' ref/annotation/human/$bedfile`
    len2=`awk '{sum+=$3-$2}END{print sum}' <(mergeBed -i ref/annotation/human/$bedfile)`

    echo $class":" $len1" "$len2
done 

# INTERSECT FUNCTIONAL ANNOTATION WITH NON-B
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "human" |while read -r sp latin filename;
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Region Class APR DR GQ IR MR TRI STR Z" |sed "s/ /\t/g" >functional/'${sp}'/enrichment_fullgenome.tsv
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        echo "Looking at $class" 
        tmp="genome "$class
        len=`mergeBed -i ref/annotation/'$sp'/$bedfile |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
        echo "Length of $class is $len"
        cat densities/'${sp}'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
        do
            echo "looking at $non_b"
            d=`intersectBed -a <(mergeBed -i ref/annotation/human/$bedfile) -b <(cat nonB_annotation/'${sp}'_pri/genome_${non_b}.bed) -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
            tmp=`echo $tmp" "$d`
            echo $tmp >tmp
        done
        cat tmp |sed "s/ /\t/g" >>functional/'${sp}'/enrichment_fullgenome.tsv
    done
    ' | sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --partition=open
done

# Plot with R/plot_fig5_functional_enrichment.R


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CHECK SIGNIFICANCE OF THE BARS IN FIGURE 5. 

# Subsample the functional regions to get a distribution of the enrichment
# values. If distribution overlaps with 1 (genome average), we deem this 
# comparison not to be significally different from the genome average.

# Subsample X sets 
sp="human"
mkdir -p functional/$sp/resample/
perc=50

cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
do 
    num=`awk -v p=$perc '{}END{num=int(NR*(1/p)); print num}' ref/annotation/$sp/$bedfile`
    # Subsample the sequences 10 times 
    echo '#!/bin/bash
    module load bedtools/2.31.0
    for i in {11..20}
    do 
        shuf -n '$num' ref/annotation/'$sp'/'$bedfile' |sort -k1,1 -k2,2n |mergeBed -i -  >functional/'$sp'/resample/Resamp.'$perc'perc.'$class'.$i.bed
    done
    ' |sbatch -J $class -o slurm-"$perc"perc-$class.%j.out --requeue --ntasks=1 --mem-per-cpu=1G --cpus-per-task=1 --time=1:00:00 --partition=open
done 

# Calculate the enrichment for each of the subsamples
sp="human"
for subset in "50perc" # "10perc" # 
do
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        echo '#!/bin/bash
        module load bedtools/2.31.0
        for i in {11..20}
        do 
            echo "Class NonB Enrichment" |sed "s/ /\t/g" >functional/'$sp'/resample/'$subset'.'$class'.$i.summary.txt
            len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' functional/'$sp'/resample/Resamp.'$subset'.'$class'.$i.bed`
            cat densities/'${sp}'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
            do
                d=`intersectBed -a functional/'$sp'/resample/Resamp.'$subset'.'$class'.$i.bed -b <(cat nonB_annotation/'${sp}'_pri/genome_${non_b}.bed) -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
                echo '$class' $non_b $d >>functional/'$sp'/resample/'$subset'.'$class'.$i.summary.txt
            done
        done 
        ' |sbatch -J $class.$subset -o slurm-$subset-$class.%j.out --requeue --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --partition=open
    done 
done

# Merge the results
# This is done separately from above since the jobs are run in parallell
sp="human"
rep=20
for subset in "50perc" "10perc" #
do
    echo "Rep Class NonB Enrichment" |sed "s/ /\t/g" >functional/$sp/$subset.${rep}rep.summary.txt
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do  
        for i in $( seq 1 $rep )
        do 
         awk '(NR>1){print "'$i'",$0}' functional/$sp/resample/$subset.$class.$i.summary.txt |sed "s/ /\t/g" >>functional/$sp/$subset.${rep}rep.summary.txt
       done 
    done 
done

# Get the min and the max values for each class and non-B type
sp="human"
rep=20
for subset in "10perc" "50perc" #
do
    echo "Class NonB Min Max" |sed "s/ /\t/g" >functional/$sp/$subset.${rep}rep.minmax.txt
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        for nonb in APR DR GQ IR MR TRI STR Z
        do 
            min=`grep $class functional/$sp/$subset.${rep}rep.summary.txt |awk -v nb=$nonb -v OFS="\t" '($3==nb){print $4}' |sort -n |head -n1`
            max=`grep $class functional/$sp/$subset.${rep}rep.summary.txt |awk -v nb=$nonb -v OFS="\t" '($3==nb){print $4}' |sort -n |tail -n1`
            echo $class $nonb $min $max |sed "s/ /\t/g" >>functional/$sp/$subset.${rep}rep.minmax.txt
        done 
    done 
done

# For the 100 runs, we actually want to remove the ~5% most extreme (meaning
# we remove the top and bottom two values, saving a 96% CI)
sp="human"
rep=100
for subset in "10perc"
do
    echo "Class NonB Min Max" |sed "s/ /\t/g" >functional/$sp/$subset.${rep}rep.96CI.txt
    cat T2T_primate_nonB/helpfiles/functional_classes.txt |while read -r class bedfile
    do 
        for nonb in APR DR GQ IR MR TRI STR Z
        do 
            min=`grep $class functional/$sp/$subset.${rep}rep.summary.txt |awk -v nb=$nonb -v OFS="\t" '($3==nb){print $4}' |sort -n |head -n3 |tail -n1`
            max=`grep $class functional/$sp/$subset.${rep}rep.summary.txt |awk -v nb=$nonb -v OFS="\t" '($3==nb){print $4}' |sort -n |tail -n3 |head -n1`
            echo $class $nonb $min $max |sed "s/ /\t/g" >>functional/$sp/$subset.${rep}rep.96CI.txt
        done 
    done 
done


# Plot with R/plot_fig10_functional_enrichment_errorBars.R


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
        bedtools nuc -fi ref/assemblies/'$filename' -bed <(cut -f1-3 ref/annotation/'$sp'/'$bedfile' |mergeBed -i -) >ref/annotation/'$sp'/$newfile
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

# Then go through the enrichment file and correct the GQ value (GQ is in the 5th column)
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
# There are no clear relationships with GC for the non-B types other than G4!

# Make a supplementary figure where only G4 is corrected for GC content
# Plot with R/plot_figS11_enrichment_functional_GCcorrected.R



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exploring enhancers and promoters, are there overlaps in the annotation?

# Number of promoters 
 wc ref/annotation/human/promoter_protCode_singleTrans.bed
# 19578 117468 851593 ref/annotation/human/promoter_protCode_singleTrans.bed
# Total length of promoters (mean length should be close to 1000 by def)
awk '{sum+=$3-$2}END{print sum}' ref/annotation/human/promoter_protCode_singleTrans.bed
# 19234458
# Mean length (slightly shorter than 1000)
awk '{sum+=$3-$2}END{print sum/NR}' ref/annotation/human/promoter_protCode_singleTrans.bed
982.453

# Number of enhancers
wc ref/annotation/human/enhancer_protCode_singleTrans.bed 
# 106998  641988 5894744 ref/annotation/human/enhancer_protCode_singleTrans.bed
#Total length of enhancers
awk '{sum+=$3-$2}END{print sum}' ref/annotation/human/enhancer_protCode_singleTrans.bed 
# 55359339
# Mean length 
awk '{sum+=$3-$2}END{print sum/NR}' ref/annotation/human/enhancer_protCode_singleTrans.bed
# 517

 
#Intersect: 
intersectBed -a ref/annotation/human/enhancer_protCode_singleTrans.bed -b ref/annotation/human/promoter_protCode_singleTrans.bed -wo |wc
# 11261  146393 1154808
# Size of intersection

intersectBed -a ref/annotation/human/enhancer_protCode_singleTrans.bed -b ref/annotation/human/promoter_protCode_singleTrans.bed |awk '{sum+=$3-$2}END{print sum}' 
# 3075721

# Not so large proportion of enhancers overlapping with promoters.

# Check within each of these categories how much overlap there is between GQ and Z-DNA. 
# First create bed files of each set of regions
intersectBed -a ref/annotation/human/promoter_protCode_singleTrans.bed -b nonB_annotation/human_pri/genome_GQ.bed |cut -f1-3 >functional/human/overlap/test.promoter.GQ.bed
intersectBed -a ref/annotation/human/promoter_protCode_singleTrans.bed -b nonB_annotation/human_pri/genome_Z.bed |cut -f1-3 >functional/human/overlap/test.promoter.Z.bed
intersectBed -a ref/annotation/human/enhancer_protCode_singleTrans.bed -b nonB_annotation/human_pri/genome_GQ.bed |cut -f1-3 >functional/human/overlap/test.enhancer.GQ.bed
intersectBed -a ref/annotation/human/enhancer_protCode_singleTrans.bed -b nonB_annotation/human_pri/genome_Z.bed |cut -f1-3 >functional/human/overlap/test.enhancer.Z.bed

# Checking sizes of GQ and Z in promoters 
awk '{sum+=$3-$2}END{print sum}' functional/human/overlap/test.promoter.GQ.bed
# 1012270
awk '{sum+=$3-$2}END{print sum}' functional/human/overlap/test.promoter.Z.bed 
# 86799
# Then intersect the regions - only a TINY fraction is overlapping 
 intersectBed -a functional/human/overlap/test.promoter.GQ.bed -b functional/human/overlap/test.promoter.Z.bed |  awk '{sum+=$3-$2}END{print sum}' 
# 1379

# Checking sizes of GQ and Z in enhancers 
awk '{sum+=$3-$2}END{print sum}' functional/human/overlap/test.enhancer.GQ.bed
#  1876727
awk '{sum+=$3-$2}END{print sum}' functional/human/overlap/test.enhancer.Z.bed
#  203501 
# Intersect the regions - even less, in terms of fraction.
intersectBed -a functional/human/overlap/test.enhancer.GQ.bed -b functional/human/overlap/test.enhancer.Z.bed |  awk '{sum+=$3-$2}END{print sum}'
# 1545

# Checking number of enhancers and promoters that has both GQ and Z-DNA 
# (even if not overlapping)
paste functional/human/overlap/promoter_protCode_singleTrans.GQ.bed functional/human/overlap/promoter_protCode_singleTrans.Z.bed |awk '($5>0 && $11>0){print}' |wc
# 2805   33660  271785
paste functional/human/overlap/enhancer_protCode_singleTrans.GQ.bed functional/human/overlap/enhancer_protCode_singleTrans.Z.bed |awk '($5>0 && $11>0){print}' |wc
# 4360   52320  494837
# This corresponds to 4.1% and 14.4% respectively


# Checking number of promoters that has both GQ and Z DNA annotated for the same bp
intersectBed -a functional/human/overlap/test.promoter.GQ.bed -b functional/human/overlap/test.promoter.Z.bed | uniq |intersectBed -a ref/annotation/human/promoter_protCode_singleTrans.bed -b - |cut -f4 |uniq |wc -l
# 205
# And enhancers 
intersectBed -a functional/human/overlap/test.enhancer.GQ.bed -b functional/human/overlap/test.enhancer.Z.bed | uniq |intersectBed -a ref/annotation/human/enhancer_protCode_singleTrans.bed -b - |cut -f4 |uniq |wc -l
# 221
# This is 1.0% and 0.2% of the total number of promoters and enhancers respectively.

# Number of overlapping bp between GQ and Z-DNA in promoters and enhancers
intersectBed -a functional/human/overlap/test.promoter.GQ.bed -b functional/human/overlap/test.promoter.Z.bed |uniq | awk '{sum+=$3-$2}END{print sum}'
#1152
intersectBed -a functional/human/overlap/test.enhancer.GQ.bed -b functional/human/overlap/test.enhancer.Z.bed |uniq | awk '{sum+=$3-$2}END{print sum}'
#1382
# We want to check how much this is in terms of GQ and Z-DNA regions annotated, 
# dividing on the number given above:
# 1152/1012270=0.0011 Of GQ in promoters 
# 1382/1876727=0.0007 Of GQ in enhancers 
# 1152/86799=0.0133 Of Z in promoters
# 1382/203501=0.0068 Of Z in enhancers

# As we also mentioned Origins of replication in the paper, we should also check
# how much overlap there is with these.
intersectBed -a ref/annotation/human/core_stochastic_origins_chm13v2.0_hglft.ucsc.bed  -b nonB_annotation/human_pri/genome_GQ.bed |cut -f1-3 >functional/human/overlap/test.origin.GQ.bed
intersectBed -a ref/annotation/human/core_stochastic_origins_chm13v2.0_hglft.ucsc.bed  -b nonB_annotation/human_pri/genome_Z.bed |cut -f1-3 >functional/human/overlap/test.origin.Z.bed
# Checking sizes of GQ and Z in origins of replication 
awk '{sum+=$3-$2}END{print sum}' functional/human/overlap/test.origin.GQ.bed
#  4212032
awk '{sum+=$3-$2}END{print sum}' functional/human/overlap/test.origin.Z.bed
#  649601
# Intersect the regions:
intersectBed -a functional/human/overlap/test.origin.GQ.bed -b functional/human/overlap/test.origin.Z.bed |  awk '{sum+=$3-$2}END{print sum}'
# 6422
# Almost nothing in terms of fraction! 