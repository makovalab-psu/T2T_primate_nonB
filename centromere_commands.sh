################################################################################
# CODE FOR CALCULATING ENRICHMENT IN CENTROMERES AND CENTROMERE FLANKS. ALSO 
# LOOKS INTO CENP-B PRESENCE, AND A DEEPER LOOK INTO THE BORNEAN ORANGUTAN 
# CHROMOSOME 10 WHERE ONE HAPLOTYPE LACKS A CENTROMERE. 


####################### ENRICHMENT IN CENTROMERES ##############################
# written by Linnéa Smeds, October 2024
# Enrichment and depletion of non-B DNA in centromeres 
# At the point of writing, Siamang had no centromeres annotated in their 
# GenomeFeature file, and were therefore excluded from this analysis. 

# Genome feature files downloaded from GenomeArk. 
mkdir ref/centromeres
mkdir ref
mkdir -p centromeres/human
mkdir -p centromeres/chimp
mkdir -p centromeres/gorilla
mkdir -p centromeres/bonobo
mkdir -p centromeres/borang
mkdir -p centromeres/sorang

# Convert from bigbed
bigBedToBed ref/centromeres/mPanPan1_v2.0.GenomeFeature_v1.0.bb ref/centromeres/mPanPan1_v2.0.GenomeFeature_v1.0.bed
bigBedToBed ref/centromeres/mPanTro3_v2.0.GenomeFeature_v0.9.bb ref/centromeres/mPanTro3_v2.0.GenomeFeature_v0.9.bed
bigBedToBed ref/centromeres/mGorGor1_v2.0.GenomeFeature_v0.9.bb ref/centromeres/mGorGor1_v2.0.GenomeFeature_v0.9.bed
bigBedToBed ref/centromeres/mPonAbe1_v2.0.GenomeFeature_v1.0.bb ref/centromeres/mPonAbe1_v2.0.GenomeFeature_v1.0.bed
bigBedToBed ref/centromeres/mPonPyg2_v2.0.GenomeFeature_v0.9.bb ref/centromeres/mPonPyg2_v2.0.GenomeFeature_v0.9.bed
bigBedToBed ref/centromeres/mSymSyn1_v2.0.GenomeFeature_v0.9.bb ref/centromeres/mSymSyn1_v2.0.GenomeFeature_v0.9.bed

# In the human file Centromere is annotated with "CEN", in the others with "Cen"
# The siamang file does not have any centromeres annotated.
# Divide into haplotypes 
for hap in "pri" "alt" #
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        mkdir -p centromeres/${sp}_$hap/
        cenname=`echo $filename |cut -f1 -d"."`
        rm -f centromeres/${sp}_$hap/cen.bed
        cenfile=$(ls ref/centromeres/$cenname*GenomeFeature*.bed)
        ls $cenfile
        for chr in $(cut -f1 ref/$filename.fai) 
        do 
            awk -v OFS="\t" -v c=$chr '($1==c && ($4 ~ /Cen/ || $4 ~ /CEN/)){print $1,$2,$3,$4}' $cenfile >>centromeres/${sp}_$hap/cen.bed
        done
    done
 done

# MIGHT NOT USE THIS 
# Calculate the lengths of the centromeres for each chromosome type
for hap in "pri" "alt" #
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        grep "chrX" centromeres/${sp}_$hap/cen.bed | awk -v sum=0 '{sum+=$3-$2}END{print sum}' >centromeres/${sp}_$hap/chrX.cen_len.txt
        grep "chrY" centromeres/${sp}_$hap/cen.bed | awk -v sum=0 '{sum+=$3-$2}END{print sum}' >centromeres/${sp}_$hap/chrY.cen_len.txt
        grep -v "chrY" centromeres/${sp}_$hap/cen.bed |grep -v "chrX" | awk -v sum=0 '{sum+=$3-$2}END{print sum}' >centromeres/${sp}_$hap/autosomes.cen_len.txt
    done
done 

# Enrichment per chromosome
module load bedtools/2.31.0
for hap in "pri"  "alt" 
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        echo "Chr APR DR GQ IR MR TRI STR Z" |sed "s/ /\t/g" >centromeres/'$sp'_'$hap'/enrichment_per_chrom.tsv
        for chr in $(cut -f1 centromeres/'$sp'_'$hap'/cen.bed |uniq)
        do
            c_len=`awk -v c=$chr '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' centromeres/'$sp'_'$hap'/cen.bed`
            chralias=`echo $chr |cut -f1 -d"_"`
            tmp="$chralias"
            cat densities/'$sp'_'$hap'_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
            do
                d=`intersectBed -a <(awk -v c=$chr '"'"'($1==c){print}'"'"' centromeres/'$sp'_'$hap'/cen.bed) -b nonB_annotation/'$sp'_'$hap'/genome_${non_b}.bed -wo |awk -v l=$c_len -v dtot=$dens '"'"'{sum+=$8}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
                tmp=`echo $tmp" "$d`
                echo $tmp >tmp.$SLURM_JOB_ID
            done
            cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/'$sp'_'$hap'/enrichment_per_chrom.tsv
        done
        '| sbatch -J $sp.$hap --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
     done
done


# Merge the species for plotting (use only "chr*", not chr*_hap*_hsa*)
for hap in "pri"  "alt" 
do
    echo "Species Chr APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g'  >centromeres/${hap}_6sp_enrichment_merged.tsv
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        awk -v sp=$sp '(NR>1){split($1,s,"_"); $1=s[1]; print sp,$0}' centromeres/${sp}_${hap}/enrichment_per_chrom.tsv |sed 's/^ //' |sed 's/ /\t/g' >>centromeres/${hap}_6sp_enrichment_merged.tsv
    done
done

# Check GC content for centromeres 
for hap in "pri"  "alt" 
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep bonobo |grep -v siamang |while read -r sp latin filename;
    do
        echo "Looking at species $sp, hap $hap"
        ls ref/assemblies/$filename 
        bedtools nuc -fi ref/assemblies/$filename -bed centromeres/${sp}_$hap/cen.bed >centromeres/${sp}_$hap/cen_with_GC.bed
    done 
done 
# And merge in a similar manner (only one value per chromosome = merge split centromeres)
for hap in "pri"  "alt" 
do
    echo "Species Chr GCcont" |sed 's/ /\t/g'  >centromeres/${hap}_6sp_GCcont_merged.tsv
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        awk -v sp=$sp '(NR>1){if(NR==2){chr=$1; gc=$8+$9; len=$13}else{if(chr==$1){gc+=$8+$9; len+=$13}else{gccont=gc/len; print sp, chr, gccont; chr=$1; gc=$8+$9; len=$13}}}END{gccont=gc/len; print sp, chr, gccont}' centromeres/${sp}_$hap/cen_with_GC.bed |awk '{split($2,s,"_"); $2=s[1]; print $0}' |sed 's/^ //' |sed 's/ /\t/g' >>centromeres/${hap}_6sp_GCcont_merged.tsv
    done
done

# Combined file with centromere locations 
for hap in "pri"  "alt" 
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        awk -v sp=$sp -v OFS='\t' '{print sp,$1,$2,$3}' centromeres/${sp}_${hap}/cen.merged.bed >>centromeres/${hap}_6sp_centromeres_merged.tsv
        awk -v sp=$sp -v OFS='\t' '{print sp,$1,$2,$3}' centromeres/${sp}_${hap}/cen.bed >>centromeres/${hap}_6sp_centromeres_unmerged.tsv
    done 
done 


##### FOR THE PAPER, CALCULATE HOW MANY CHR THAT HAS AT LEAST ONE non-B ENRICHED
# I do this in the R scripts that generate the centromere figures!


# ~~~~~~~~~~~~~~~~~ CENTROMERE DENSITIES INSTEAD OF ENRICHMENT ~~~~~~~~~~~~~~~~~
# This results were not presented in the final paper

# Density per chromosome
module load bedtools/2.31.0
for hap in "pri" "alt" 
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        echo "Chr APR DR GQ IR MR TRI STR Z" |sed "s/ /\t/g" >centromeres/'$sp'_'$hap'/density_per_chrom.tsv
        for chr in $(cut -f1 centromeres/'$sp'_'$hap'/cen.bed |uniq)
        do
            c_len=`awk -v c=$chr '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' centromeres/'$sp'_'$hap'/cen.bed`
            chralias=`echo $chr |cut -f1 -d"_"`
            tmp="$chralias"
            cat densities/'$sp'_'$hap'_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
            do
                d=`intersectBed -a <(awk -v c=$chr '"'"'($1==c){print}'"'"' centromeres/'$sp'_'$hap'/cen.bed) -b nonB_annotation/'$sp'_'$hap'/genome_${non_b}.bed -wo |awk -v l=$c_len '"'"'{sum+=$8}END{d=sum/l; print d}'"'"'`
                tmp=`echo $tmp" "$d`
                echo $tmp >tmp.$SLURM_JOB_ID
            done
            cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/'$sp'_'$hap'/density_per_chrom.tsv
        done
        '| sbatch -J $sp.$hap --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
     done
done
# Merge the species for plotting (use only "chr*", not chr*_hap*_hsa*)
for hap in "pri"  "alt" 
do
    echo "Species Chr APR DR GQ IR MR STR Z" |sed 's/ /\t/g'  >centromeres/${hap}_6sp_density_merged.tsv
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        awk -v sp=$sp '(NR>1){split($1,s,"_"); $1=s[1]; print sp,$0}' centromeres/${sp}_${hap}/density_per_chrom.tsv |sed 's/^ //' |sed 's/ /\t/g' >>centromeres/${hap}_6sp_density_merged.tsv
    done
done



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STATISTICAL TEST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Divide the rest of the chromosome into similar sized regions, remove the real 
# centromere (entire region if separated into two or more regions), and 
# calculate the non-B density for each region. 

module load bedtools/2.31.0

# We don't want to generate anything in the real centromere, or - if there are
# more than one annotated centromere - between them. Therefore, make merged 
# versions of the cen files first. 

for hap in "pri"  "alt" 
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do 
      awk -v OFS="\t" '{if(NR==1){chr=$1;st=$2;e=$3}else{if($1==chr){e=$3}else{print chr,st,e; ;chr=$1;st=$2;e=$3}}}END{print chr,st,e}' centromeres/${sp}_$hap/cen.bed >centromeres/${sp}_$hap/cen.merged.bed
    done 
done 

# Generating 100 windows of the same size as the centromere, not overlapping with
# the real centromere. Randomly selected if there are more than 100 windows,
# overlapping if there are less. For very large centromeres, removing overlap
# with the real centromere and non-complete windows at the end of the chromosome
# can result in less than 100 windows. To make sure there always are more than 
# 100, I make the overlap so there are 120 windows to choose from, and allow end
# windows to be as small as 1M (if this is smaller than the centromere size). 
for hap in "alt" "pri" 
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        mkdir -p centromeres/${sp}_$hap/windows/
        
        for chr in $(cut -f1 centromeres/${sp}_$hap/cen.bed |uniq)
        do
            echo '#!/bin/bash
            echo "Looking at '$sp' and '$chr'"
            module load bedtools/2.31.0
            c_len=`awk -v c='$chr' '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' centromeres/'${sp}'_'$hap'/cen.bed`
            chrlen=`awk -v c='$chr' '"'"'($1==c){print $2}'"'"' ref/'$filename'.fai`
            poswin=`echo $chrlen"/"$c_len |bc`
            minsize=1000000
            if [ "$poswin" -lt "120" ]
            then
              echo "For '$chr', there will be less than 120 nonovl windows"
              slide=`echo "("$chrlen"-"$c_len")/120"|bc`
              bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$filename'.fai) -s $slide -w $c_len |subtractBed -a - -b centromeres/'${sp}'_'$hap'/cen.merged.bed |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/'${sp}'_'$hap'/windows/'$chr'.exclCen.bed
            else 
              echo "For '$chr', there will be more than 120 nonovl windows"
              bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$filename'.fai) -w $c_len |subtractBed -a - -b centromeres/'${sp}'_'$hap'/cen.merged.bed |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/'${sp}'_'$hap'/windows/'$chr'.exclCen.bed
            fi
            bedtools nuc -fi ref/assemblies/'$filename' -bed centromeres/'${sp}'_'$hap'/windows/'$chr'.exclCen.bed >centromeres/'${sp}'_'$hap'/windows/'$chr'.exclCen.with_GC.bed
            rm -f centromeres/'${sp}'_'$hap'/windows/'$chr'.100rand.enrichment.tsv
            tmp="'$chr'"
            while read line
            do 
              cat densities/'${sp}'_'${hap}'_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
              do
                  echo $line |sed "s/ /\t/g" >tmp.$SLURM_JOB_ID.bed
                  d=`intersectBed -a tmp.$SLURM_JOB_ID.bed -b nonB_annotation/'$sp'_'$hap'/genome_${non_b}.bed -wo |awk -v l=$c_len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
                  tmp=`echo $tmp" "$d`
                  echo $tmp >tmp.$SLURM_JOB_ID
              done
              cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/'${sp}'_'$hap'/windows/'$chr'.100rand.enrichment.tsv
            done <centromeres/'${sp}'_'$hap'/windows/'$chr'.exclCen.bed
            '| sbatch -J $sp.$hap --ntasks=1 --cpus-per-task=1 --partition=open --time=5:00:00
        done 
    done 
done 


# Merge background densities and GC 
for hap in "pri"  "alt" 
do
  echo "Window Species Chr APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g'  >centromeres/${hap}_6sp_background_enrichment_merged.tsv
  echo "Window Species Chr GCcont" |sed 's/ /\t/g'  >centromeres/${hap}_6sp_background_GCcont_merged.tsv
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
  do
    echo "Window Chr APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g' >centromeres/${sp}_$hap/AllChr.100rand.enrichment.tsv
    for file in $(ls centromeres/${sp}_$hap/windows/chr*.100rand.enrichment.tsv)
    do
      awk '{print NR,$0}' $file |sed 's/ /\t/g' >>centromeres/${sp}_$hap/AllChr.100rand.enrichment.tsv
      awk -v sp=$sp '{split($1,s,"_"); $1=s[1]; print NR,sp,$0}' $file |sed 's/ /\t/g' >>centromeres/${hap}_6sp_background_enrichment_merged.tsv
    done 
    for file in $(ls centromeres/${sp}_$hap/windows/chr*.exclCen.with_GC.bed)
    do 
      awk -v sp=$sp '(NR>1){split($1,s,"_"); $1=s[1]; rownum=NR-1; print rownum,sp,$1,$5}' $file |sed 's/ /\t/g' >>centromeres/${hap}_6sp_background_GCcont_merged.tsv
    done 
  done 
done 

# Make an ordered chromosome file to be used for plotting 
# (chromosomes are correctly ordered in the lengths file)
for hap in "pri"  "alt" 
do
  grep "_$hap" T2T_primate_nonB/helpfiles/all_lengths.txt |grep -v siamang |grep -v "random" |sed "s/_$hap//" |awk '{split($2,s,"_"); $2=s[1]; print $0}' >centromeres/${hap}_6sp_order.tsv
done



# ~~~~~~~~~~~~~~~~~~~~ A LOOK INTO BORNEAN ORANGUTAN CHR10 ~~~~~~~~~~~~~~~~~~~~~
# Bornean orangutan has no annotated centromere in the primary chr 10. 
# Extract the two chromosome 10 from the primary and alternative assemblies 
module load samtools/1.19.2
echo "chr10_hap1_hsa12" >bornean_pri.temp 
echo "chr10_hap2_hsa12" >bornean_alt.temp 
samtools faidx ref/assemblies/mPonPyg2.pri.cur.20231122.fasta -r bornean_pri.temp -o centromeres/borang_pri/chr10.fa
samtools faidx ref/assemblies/mPonPyg2.alt.cur.20231122.fasta -r bornean_alt.temp -o centromeres/borang_alt/chr10.fa

# Align the two haplotypes to each other using winnowmap  
echo '#!/bin/bash
~/software/meryl-1.4.1/bin/meryl count k=19 output centromeres/borang_alt/chr10.meryldb centromeres/borang_alt/chr10.fa
~/software/meryl-1.4.1/bin/meryl print greater-than distinct=0.9998 centromeres/borang_alt/chr10.meryldb > centromeres/borang_alt/chr10.repeats
~/software/Winnowmap/bin/winnowmap -x asm20 -c --eqx -t 4 -W centromeres/borang_alt/chr10.repeats centromeres/borang_alt/chr10.fa centromeres/borang_pri/chr10.fa >centromeres/borang_pri/chr10_to_alt.winnowmap.paf
' | sbatch -J winnow --ntasks=1 --cpus-per-task=4 --time=5:00:00 --partition=open

# Convert to readable format 
cut -f1,3,4,6,8,9 centromeres/borang_pri/chr10_to_alt.winnowmap.paf |sort -k1,1 -k2,2n >centromeres/borang_pri/chr10_to_alt.aligned.bed
#Take only blocks that overlap with region of interest
awk '($5>28000000 && $6<40000000){print}' centromeres/borang_pri/chr10_to_alt.aligned.bed >centromeres/borang_pri/chr10_to_alt.aligned.centromere.txt


# From looking at the alignment blocks, we see that the centromere present in the
# alternative assembly is not aligning to anything in the primary assembly, this
# region seem to be entirely missing. 

#Extract non-B density in these regions 
grep chr10 densities/borang_pri_comb_100kb.bed |awk '($2>29000000 && $3<39000000){print}' >densities/borang_pri_centro_100kb.bed
grep chr10 densities/borang_alt_comb_100kb.bed |awk '($2>29000000 && $3<39000000){print}' >densities/borang_alt_centro_100kb.bed

# Plot this in R using 
# plot_figS15_Bornean_orangutan.R


# ~~~~~~~~~~~~~~~~~~~~ COMPARING CENTROMERES WITH AND WITHOUT CENP-B ~~~~~~~~~~~~~~~~~~~~~
# 
cd ref/centromeres
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/cen/mPanPan1_v2.0.cenpb_sites_v1.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/cen/mPanTro3_v2.0.cenpb_sites_v1.0.bb 
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/repeats/cen/mGorGor1_v2.0.cenpb_sites_v1.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/cen/mPonAbe1_v2.0.cenpb_sites_v1.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/repeats/cen/mPonPyg2_v2.0.cenpb_sites_v1.0.bb 
# Human file is a supplementary datafile to Altermose et al 2022 
https://www.science.org/doi/suppl/10.1126/science.abl4178/suppl_file/science.abl4178_databases_s1_to_s21.zip

# Extract big bed files 
~/software/bigBedToBed mPanPan1_v2.0.cenpb_sites_v1.0.bb mPanPan1_v2.0.cenpb_sites_v1.0.bed
~/software/bigBedToBed mPanTro3_v2.0.cenpb_sites_v1.0.bb mPanTro3_v2.0.cenpb_sites_v1.0.bed
~/software/bigBedToBed mGorGor1_v2.0.cenpb_sites_v1.0.bb mGorGor1_v2.0.cenpb_sites_v1.0.bed
~/software/bigBedToBed mPonAbe1_v2.0.cenpb_sites_v1.0.bb mPonAbe1_v2.0.cenpb_sites_v1.0.bed
~/software/bigBedToBed mPonPyg2_v2.0.cenpb_sites_v1.0.bb mPonPyg2_v2.0.cenpb_sites_v1.0.bed
cd ../.. 

# Check overlap between centromeres and CENP-B
module load bedtools/2.31.0

# Human (different filename format)
intersectBed -a centromeres/human_pri/cen.bed -b <(cut -f1-4 ref/centromeres/databaseS15_cenpB_pJalpha.bed)  -wo |cut -f1,2,3 |uniq >centromeres/human_pri/cen_with_cenpb.bed
intersectBed -v -a centromeres/human_pri/cen.bed -b centromeres/human_pri/cen_with_cenpb.bed >centromeres/human_pri/cen_without_cenpb.bed

# Non-humans
for hap in "pri" "alt"
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v human |grep -v siamang |while read -r sp latin filename;
    do 
        echo $filename 
        prefix=`echo $filename |cut -f1 -d"."`
        echo $prefix 
        intersectBed -a centromeres/${sp}_${hap}/cen.bed -b <(cut -f1-4 ref/centromeres/${prefix}_v2.0.cenpb_sites_v1.0.bed) -wo |cut -f1,2,3 |uniq >centromeres/${sp}_${hap}/cen_with_cenpb.bed
        intersectBed -v -a centromeres/${sp}_${hap}/cen.bed -b centromeres/${sp}_${hap}/cen_with_cenpb.bed >centromeres/${sp}_${hap}/cen_without_cenpb.bed
    done 
done 

#Merge files 
for hap in "pri" "alt"
do
    echo "Species Chr CENPB" | sed 's/ /\t/g' >centromeres/${hap}_6sp_cenpb.tsv
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do 
        awk -v OFS="\t" -v sp=$sp '{split($1,s,"_"); print sp, s[1], "Yes"}' centromeres/${sp}_${hap}/cen_with_cenpb.bed >>centromeres/${hap}_6sp_cenpb.tsv
        awk -v OFS="\t" -v sp=$sp '{split($1,s,"_"); print sp, s[1], "No"}' centromeres/${sp}_${hap}/cen_without_cenpb.bed >>centromeres/${hap}_6sp_cenpb.tsv
    done
done 

# SAME, BUT INCLUDE 1MB FLANKS 
# Humans
intersectBed -a <(awk -v OFS="\t" '{$2=$2-1000000; $3=$3+1000000; print}' centromeres/human_pri/cen.bed) -b <(cut -f1-4 ref/centromeres/databaseS15_cenpB_pJalpha.bed) -wo |cut -f1,2,3 |uniq >centromeres/human_pri/cen1MbFlank_with_cenpb.bed
intersectBed -v -a centromeres/human_pri/cen.bed -b centromeres/human_pri/cen1MbFlank_with_cenpb.bed >centromeres/human_pri/cen1MbFlank_without_cenpb.bed
# Non-humans
for hap in "pri" "alt"
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v human |grep -v siamang |while read -r sp latin filename;
    do 
        prefix=`echo $filename |cut -f1 -d"."`
        intersectBed -a <(awk -v OFS="\t" '{$2=$2-1000000; $3=$3+1000000; print}' centromeres/${sp}_${hap}/cen.bed) -b <(cut -f1-4 ref/centromeres/${prefix}_v2.0.cenpb_sites_v1.0.bed) -wo |cut -f1,2,3 |uniq >centromeres/${sp}_${hap}/cen1MbFlank_with_cenpb.bed
        intersectBed -v -a centromeres/${sp}_${hap}/cen.bed -b centromeres/${sp}_${hap}/cen1MbFlank_with_cenpb.bed >centromeres/${sp}_${hap}/cen1MbFlank_without_cenpb.bed
    done 
done 
#Merge 
for hap in "pri" "alt"
do
    echo "Species Chr CENPB" | sed 's/ /\t/g' >centromeres/${hap}_6sp_cenpb_incl1MbFlanks.tsv
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do 
        awk -v OFS="\t" -v sp=$sp '{split($1,s,"_"); print sp, s[1], "Yes"}' centromeres/${sp}_${hap}/cen1MbFlank_with_cenpb.bed >>centromeres/${hap}_6sp_cenpb_incl1MbFlanks.tsv
        awk -v OFS="\t" -v sp=$sp '{split($1,s,"_"); print sp, s[1], "No"}' centromeres/${sp}_${hap}/cen1MbFlank_without_cenpb.bed >>centromeres/${hap}_6sp_cenpb_incl1MbFlanks.tsv
    done
done 



# Plot and calculate statistics with plot_figS16_enrichment_and_CENPB.R



# ~~~~~~~~~~~~~~~~~~~~~~~~ LOOKING AT CENTROMERE FLANKS ~~~~~~~~~~~~~~~~~~~~~~~~
# Only 1Mb flanks were used for the paper
module load bedtools/2.31.0
for hap in "pri" "alt" 
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        awk -v OFS="\t" '{start=$2-1000000; print $1,start,$2}' centromeres/${sp}_${hap}/cen.merged.bed >centromeres/${sp}_${hap}/cen.upstream1Mb.bed
  #      awk -v OFS="\t" '{start=$2-3000000; print $1,start,$2}' centromeres/${sp}_${hap}/cen.merged.bed >centromeres/${sp}_${hap}/cen.upstream3Mb.bed
        awk -v OFS="\t" '{end=$3+1000000; print $1,$3,end}' centromeres/${sp}_${hap}/cen.merged.bed >centromeres/${sp}_${hap}/cen.downstream1Mb.bed
  #      awk -v OFS="\t" '{end=$3+3000000; print $1,$3,end}' centromeres/${sp}_${hap}/cen.merged.bed >centromeres/${sp}_${hap}/cen.downstream3Mb.bed
        # Also make files with both flanks+centromeres to subtract from background
        awk -v OFS="\t" '{start=$2-1000000; end=$3+1000000; print $1,start,end}' centromeres/${sp}_${hap}/cen.merged.bed >centromeres/${sp}_${hap}/cen.merged.with1Mbflanks.bed
  #      awk -v OFS="\t" '{start=$2-3000000; end=$3+3000000; print $1,start,end}' centromeres/${sp}_${hap}/cen.merged.bed >centromeres/${sp}_${hap}/cen.merged.with3Mbflanks.bed
         join centromeres/${sp}_${hap}/cen.merged.bed <(grep -v "random" ref/$filename.fai) |awk -v OFS='\t' '{s=$3+1000000; print $1,s,$4}' >centromeres/${sp}_${hap}/cen.q_excl1Mbflank.bed
         join centromeres/${sp}_${hap}/cen.merged.bed <(grep -v "random" ref/$filename.fai) |awk -v OFS='\t' '{e=$2-1000000; print $1,0,e}' >centromeres/${sp}_${hap}/cen.p_excl1Mbflank.bed
   
    done 
done 

# Enrichment per chromosome
module load bedtools/2.31.0
for flank in  "upstream1Mb" "downstream1Mb" #"q_excl1Mbflank" "p_excl1Mbflank" #"upstream3Mb"  "downstream3Mb"
do 
    for hap in "pri"  "alt" 
    do
        cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
        do
            echo '#!/bin/bash
            module load bedtools/2.31.0
            echo "Chr APR DR GQ IR MR TRI STR Z" |sed "s/ /\t/g" >centromeres/'$sp'_'$hap'/'$flank'.enrichment_per_chrom.tsv
            for chr in $(cut -f1 centromeres/'$sp'_'$hap'/cen.'$flank'.bed |uniq)
            do
                c_len=`awk -v c=$chr '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' centromeres/'$sp'_'$hap'/cen.'$flank'.bed`
                chralias=`echo $chr |cut -f1 -d"_"`
                tmp="$chralias"
                cat densities/'$sp'_'$hap'_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
                do
                    d=`intersectBed -a <(awk -v c=$chr '"'"'($1==c){print}'"'"' centromeres/'$sp'_'$hap'/cen.'$flank'.bed) -b nonB_annotation/'$sp'_'$hap'/genome_${non_b}.bed -wo |awk -v l=$c_len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
                    tmp=`echo $tmp" "$d`
                    echo $tmp >tmp.$SLURM_JOB_ID
                done
                cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/'$sp'_'$hap'/'$flank'.enrichment_per_chrom.tsv
            done
            '| sbatch -J $sp.$hap --ntasks=1 --cpus-per-task=1 --partition=open
        done
    done
done 

# Merge the species for plotting (use only "chr*", not chr*_hap*_hsa*)
for flank in "upstream1Mb" "downstream1Mb"  #"q_excl1Mbflank" "p_excl1Mbflank" # "upstream3Mb"  "downstream3Mb"
do 
    for hap in "pri"  "alt" 
    do
        echo "Species Chr APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g'  >centromeres/${hap}_6sp_${flank}_enrichment_merged.tsv
        cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
        do
            awk -v sp=$sp '(NR>1){split($1,s,"_"); $1=s[1]; print sp,$0}' centromeres/${sp}_${hap}/${flank}.enrichment_per_chrom.tsv |sed 's/^ //' |sed 's/ /\t/g' >>centromeres/${hap}_6sp_${flank}_enrichment_merged.tsv
        done
    done
done 



# ~~~~~~~~~~~~~~~~ ADD SATELLITE ANNOTATION TO CENTROMERE PLOTS ~~~~~~~~~~~~~~~~
# In the CenSat annotation from the human genome, the categories are as follows:
# (Tabel from genome browser)
# Active αSat HOR (hor ... L)	red
# Inactive αSat HOR (hor)	orange
# Divergent αSat HOR (dhor)	dark red
# Monomeric αSat (mon)	peach/yellow
# Classical Human Satellite 1A (hsat1A)	light green
# Classical Human Satellite 1B (hsat1B)	dark green
# Classical Human Satellite 2 (hsat2)	light blue
# Classical Human Satellite 3 (hsat3)	blue
# Beta Satellite (bsat)	pink
# Gamma Satellite (gsat)	purple
# Other centromeric satellites (censat)	teal
# Centromeric transition regions (ct)	grey

# Downloading the censat annotation
cd ref/repeats/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/repeats/mGorGor1_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/repeats/mPonPyg2_v2.0_CenSat_v2.0.bb
wget https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/repeats/mSymSyn1_v2.1_CenSat_v2.0.bb 
~/software/bigBedToBed mGorGor1_v2.0_CenSat_v2.0.bb mGorGor1_v2.0_CenSat_v2.0.bed 
~/software/bigBedToBed mPanPan1_v2.0_CenSat_v2.0.bb mPanPan1_v2.0_CenSat_v2.0.bed
~/software/bigBedToBed mPanTro3_v2.0_CenSat_v2.0.bb mPanTro3_v2.0_CenSat_v2.0.bed
~/software/bigBedToBed mPonAbe1_v2.0_CenSat_v2.0.bb mPonAbe1_v2.0_CenSat_v2.0.bed
~/software/bigBedToBed mPonPyg2_v2.0_CenSat_v2.0.bb mPonPyg2_v2.0_CenSat_v2.0.bed
~/software/bigBedToBed mSymSyn1_v2.1_CenSat_v2.0.bb mSymSyn1_v2.1_CenSat_v2.0.bed
cd ../..

#Human 
awk -v OFS="\t" '(NR>1){if($4~/^hor/){if($4~/L)$/){print $1,$2,$3,"active_asat"}else{print $1,$2,$3,"inactive_asat"}}else{split($4,s,"_"); print $1,$2,$3,tolower(s[1])}}' ref/repeats/chm13v2.0_censat_v2.0.bed >centromeres/human_pri/cenSat_shortnames.bed

# Others
for hap in "pri" "alt" #
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |while read -r sp latin filename;
    do
        name=`echo $filename |cut -f1 -d"."`
        cut -f1 ref/$filename.fai |grep -v random |grep -f - ref/repeats/$name*CenSat*.bed |cut -f1 -d "(" |cut -f1,2,3,4 |awk -v OFS="\t" '{if($4=="active_hor"){print $1,$2,$3,"active_asat"}else if($4=="hor"){print $1,$2,$3,"inactive_asat"}else{split($4,s,"_"); split(s[1],t,","); print $1,$2,$3,tolower(t[1])}}' >centromeres/${sp}_$hap/cenSat_shortnames.bed
    done 
done 

#Merge the different species per haplotype
for hap in "pri"  "alt" 
do
    echo "Species Chr Sat_start Sat_stop Sat" |sed 's/ /\t/g' >centromeres/${hap}_6sp_cenSat_merged.tsv
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v siamang |while read -r sp latin filename;
    do
        awk -v OFS="\t" -v sp=$sp '{split($1,s,"_"); $1=s[1]; print sp,$0}' centromeres/${sp}_${hap}/cenSat_shortnames.bed >>centromeres/${hap}_6sp_cenSat_merged.tsv
    done 
done 



# ~~~~~~~~~~~~~~~~ CENTROMERE SUPRACHROMOSOMAL FAMILIES ~~~~~~~~~~~~~~~~
# Check the active centromere superfamilies 
sed 's/(/\t/g' chm13v2.0_censat_v2.0.bed |sed 's/)/\t/g' | awk -v OFS="\t" '(NR>1){if($4~/hor/){if($5~/L$/){print $1,$2,$3,$4,$5}}}' |less

# Non-humans
for hap in "pri" "alt" 
do
    echo "Species LongChr Chr SF" |sed 's/ /\t/g' >centromeres/${hap}_6sp_SF_merged.tsv
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v human |grep -v siamang |while read -r sp latin filename;
    do
        echo "Looking at $sp"
        name=`echo $filename |cut -f1 -d"."`
        cut -f1 ref/$filename.fai |grep -v random |grep -f - ref/repeats/$name*CenSat*.bed | grep "active_hor" |cut -f1,4 |sed 's/(/\t/' |sed 's/)/\t/' |awk -v OFS="\t" -v sp=$sp '{split($1,s,"_"); print sp,$1,s[1],$3}' >>centromeres/${hap}_6sp_SF_merged.tsv
    done
done 
# Human file is different from the other files
awk -v OFS="\t" '(NR>1){if($4~/^hor/ && $4~/L)$/){print}}' ref/repeats/chm13v2.0_censat_v2.0.bed  |cut -f1,4 |sed 's/(/\t/' |sed 's/\//\t/' |sed 's/C/\t/' |cut -f1,3 |sed 's/S/SF/g' |awk '{print "human\t"$1,"\t"$1,"\t",$2}' >>centromeres/pri_6sp_SF_merged.tsv

# Combine SF annotation with non-B enrichment in R:
# plot_figS14_enrichment_and_SF.R