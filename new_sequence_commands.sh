################################################################################
# CODE FOR COMPARING THE T2T PRIMATE NON-B DNA ANNOTATION WITH OLD ASSEMBLIES



################### ENRICHMENT OF NON-B DNA IN NEW SEQUENCE ####################
# written by Linnéa Smeds, September 2024
# This code was also used for the primates autosome paper (Yoo et al 2024) 
# Code for winnowmap alignments was adapted from Bob Harris (from Makova et 
# al. Nature 2024)

# Old and new assemblies divided into separate files for each chromosomes,
# named after the new chromosome number, example: "new_sequence/chimp/T2T.chr2.fa" 
# and "new_sequence/chimp/old.chr2.fa" (the latter contains chr3 from the old 
# chimpanzee assembly)

# DOWLOAD ASSEMBLIES 
mkdir ref/assemblies/
cd ref/assemblies/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fa
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.alt.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.alt.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.pri.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.alt.cur.20231122.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.pri.cur.20231205.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.alt.cur.20231205.fasta.gz
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.pri.cur.20231122.fasta.gz 
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.alt.cur.20231122.fasta.gz 
wget https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/mSymSyn1.pri.cur.20240514.fasta.gz 
wget https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/mSymSyn1.alt.cur.20240514.fasta.gz 
cd ../.. 
# Already had the indexes from before 
cp ref/*fai ref/assemblies/

#OLD ASSEMBLIES 
mkdir ref/assemblies/old_versions
cd ref/assemblies/old_versions
wget https://hgdownload.soe.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.fa.gz
wget https://hgdownload2.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.fa.gz
wget https://hgdownload2.soe.ucsc.edu/goldenPath/panPan3/bigZips/panPan3.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/ponAbe3/bigZips/ponAbe3.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/015/021/865/GCA_015021865.1_gorGor.msY.makovalab.ver3/GCA_015021865.1_gorGor.msY.makovalab.ver3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/015/021/855/GCA_015021855.1_panPan.msY.makovalab.ver1/GCA_015021855.1_panPan.msY.makovalab.ver1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/015/021/835/GCA_015021835.1_ponAbe.msY.makovalab.ver3/GCA_015021835.1_ponAbe.msY.makovalab.ver3_genomic.fna.gz


# Unzip and index old 
for $file in $(ls ref/assemblies/old_versions/*.f*a.gz)
do
    echo '#!/bin/bash
    module load samtools/1.19.2
    name=`basename '$file' |sed "s/.gz//"`
    gunzip -c '$file' >ref/assemblies/old_versions/$name
    samtools faidx ref/assemblies/old_versions/$name
    ' |sbatch -J $file --ntasks=1 --cpus-per-task=4 --time=1:00:00 --partition=open
done 
cd ../../..

# Unzip and index new 
for file in $(ls ref/assemblies/*.fasta.gz)
do 
    echo '#!/bin/bash
        module load samtools/1.19.2
        name=`basename '$file' |sed "s/.gz//"`
        echo $name 
        gunzip -c '$file' >ref/assemblies/$name 
        samtools faidx ref/assemblies/$name 
     ' |sbatch -J $file --ntasks=1 --cpus-per-task=4 --time=1:00:00 --partition=open
done 


# Extract fasta sequence 
for hap in "pri" "alt" 
do 
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep "human" |grep -v "borang" |grep -v siamang | while read -r sp latin filename;
    do 
        mkdir -p new_sequence/${sp}_$hap
        mkdir -p new_sequence/${sp}_old
        old_ass=`grep $sp T2T_primate_nonB/helpfiles/old_assemblies.txt |cut -f2 -d" "`
        echo $sp $old_ass
        # Go through all chromosome pairs (old might consist of more than one chr)
        for chr in $(cut -f1 T2T_primate_nonB/helpfiles/chr_translation/${sp}_${hap}.txt |grep -v "#" |uniq)
        do 
            echo $chr
            # Get new chromosome
            awk -v c=$chr '($1==c){print $2}' T2T_primate_nonB/helpfiles/chr_translation/${sp}_${hap}.txt |uniq >$sp.$hap.temp
            samtools faidx ref/assemblies/$filename -r $sp.$hap.temp -o new_sequence/${sp}_$hap/t2t.$chr.fa
            # get old chromosome(s)
            awk -v c=$chr '($1==c){print $3}' T2T_primate_nonB/helpfiles/chr_translation/${sp}_${hap}.txt |uniq >$sp.$hap.temp
            samtools faidx ref/assemblies/old_versions/$old_ass -r $sp.$hap.temp -o new_sequence/${sp}_old/old.$chr.fa
        done 
    done 
done 
# Gorilla, Bonobo and Sumatran orangutan does not have Y chromosomes in these
# assembly versions, we use the best available non-T2T sequences. 
# Bonobo: GCA_015021855.1_panPan.msY.makovalab.ver1_genomic.fna
# Gorilla: GCA_015021865.1_gorGor.msY.makovalab.ver3_genomic.fna
# Sumatran orangutang: GCA_015021835.1_ponAbe.msY.makovalab.ver3_genomic.fna
awk '{if(/>/){print $1}else{print}}' ref/assemblies/old_versions/GCA_015021855.1_panPan.msY.makovalab.ver1_genomic.fna >new_sequence/bonobo_pri/old.chrY.fa 
awk '{if(/>/){print $1}else{print}}' ref/assemblies/old_versions/GCA_015021865.1_gorGor.msY.makovalab.ver3_genomic.fna >new_sequence/gorilla_pri/old.chrY.fa 
awk '{if(/>/){print $1}else{print}}' ref/assemblies/old_versions/GCA_015021835.1_ponAbe.msY.makovalab.ver3_genomic.fna >new_sequence/sorang_pri/old.chrY.fa 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ALIGN SEQUENCE: Run meryl and winnowmap 
for sp in "human" "chimp" "bonobo" "gorilla" "sorang"  
do
  for i in {1..23} "X" "Y" 
  do
    job=$(echo '#!/bin/bash
    ~/software/meryl-1.4.1/bin/meryl count k=19 output new_sequence/'$sp'_old/old.chr'$i'.meryldb new_sequence/'$sp'_old/old.chr'$i'.fa
    ~/software/meryl-1.4.1/bin/meryl print greater-than distinct=0.9998 new_sequence/'$sp'_old/old.chr'$i'.meryldb > new_sequence/'$sp'_old/old.chr'$i'.repeats
    ' | sbatch -J meryl.$sp.chr$i --ntasks=1 --cpus-per-task=1 --time=05:00 --partition=open |cut -f4 -d" ")
    for hap in "pri" "alt" 
    do 
      echo '#!/bin/bash
      ~/software/Winnowmap/bin/winnowmap -x asm20 -c --eqx -t 4 -W new_sequence/'$sp'_old/old.chr'$i'.repeats new_sequence/'$sp'_old/old.chr'$i'.fa new_sequence/'$sp'_'$hap'/t2t.chr'$i'.fa >new_sequence/'$sp'_'$hap'/chr'$i'.winnowmap.paf
      ' | sbatch -J $i.$sp.$hap.winnowmap --ntasks=1 --cpus-per-task=4 --time=5:00:00 --partition=open -d afterok:$job
    done
  done
done
# NOTE! This analysis extracts sequences that does not align at all to the 
# older assembly versions. In case of former collapsed repeats that are still
# represented in one copy in the old assembly, several repeat copies can align 
# many-to-one. Using a one-to-one approach would result in even more new seq.

# Convert to bed
for hap in "pri" "alt" 
do 
   cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "borang" |grep -v siamang | while read -r sp latin filename;
  do
    echo "Looking at ${sp}_$hap"
    rm -f new_sequence/${sp}_$hap/merged_aligned.bed
    for i in {1..23} "X" "Y" 
    do
      if [ -e new_sequence/${sp}_$hap/chr${i}.winnowmap.paf ]
      then
        echo Looking at chr$i
        cut -f1,3,4 new_sequence/${sp}_$hap/chr${i}.winnowmap.paf |sort -k1,1 -k2,2n |grep -v "random" |mergeBed -i - |awk -v c=$i '{print $1"\t"$2"\t"$3}' >new_sequence/${sp}_$hap/chr${i}.aligned.bed
        cat new_sequence/${sp}_$hap/chr${i}.aligned.bed >>new_sequence/${sp}_$hap/merged_aligned.bed
      fi
    done
    echo "get complement"
    cut -f1,2 ref/$filename.fai |grep -v "random" |grep -v "chrM" |sort -k1,1 | complementBed -i <(sort -k1,1 -k2,2n new_sequence/${sp}_$hap/merged_aligned.bed) -g - >new_sequence/${sp}_$hap/merged_unaligned.bed
  done
done

# Separate merged bed files into autosomes and sex chromosomes
# (to simplify downstream analysis) 
for hap in "pri" "alt" 
do 
   cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "borang" |grep -v siamang | while read -r sp latin filename;
  do
    grep -v chrX new_sequence/${sp}_$hap/merged_unaligned.bed |grep -v chrY |sort -k1,1 -k2,2n >new_sequence/${sp}_$hap/autosomes.unaligned.bed 
    grep -v chrX new_sequence/${sp}_$hap/merged_aligned.bed |grep -v chrY >new_sequence/${sp}_$hap/autosomes.aligned.bed 
    grep chrX new_sequence/${sp}_$hap/merged_unaligned.bed >new_sequence/${sp}_$hap/chrX.unaligned.bed 
    grep chrY new_sequence/${sp}_$hap/merged_unaligned.bed >new_sequence/${sp}_$hap/chrY.unaligned.bed 
  done 
done 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DENSITY AND STATS CALCULATIONS

# Count 'new' (=unaligned) basepairs for human:
# Autosomes 
grep -v chrX new_sequence/human_pri/merged_unaligned.bed |grep -v chrY |awk '{sum+=$3-$2}END{print sum}'
# 88575866
grep chrX new_sequence/human_pri/merged_unaligned.bed |awk '{sum+=$3-$2}END{print sum}'
# 136784
grep chrY new_sequence/human_pri/merged_unaligned.bed |awk '{sum+=$3-$2}END{print sum}'
# 15105017

# Calculate non-B density inside and outside new sequence
for hap in "pri" "alt" 
do 
   cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "borang" |grep -v siamang | while read -r sp latin filename;
  do
    echo "looking at $sp"
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "#Chr NonB NewTotBp NewNonB NewDens OldTotBp OldNonB OldDens" >new_sequence/'$sp'_'$hap'/summary.txt
    for chr in "autosomes" "chrX" "chrY"
      do
      for nb in "APR" "DR" "STR" "IR" "MR" "TRI" "GQ" "Z" "all"
      do
        new=`intersectBed -a new_sequence/'$sp'_'$hap'/$chr.unaligned.bed -b nonB_annotation/'$sp'_'$hap'/${chr}_${nb}.bed  -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk -v nb=$nb '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print nb,sum_l,sum_nb,d}'"'"'`
        old=`intersectBed -a new_sequence/'$sp'_'$hap'/$chr.aligned.bed -b nonB_annotation/'$sp'_'$hap'/${chr}_${nb}.bed -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print sum_l,sum_nb,d}'"'"'`
        echo $chr" "$new" "$old >>new_sequence/'$sp'_'$hap'/summary.txt
      done 
    done
     ' |sbatch -J $hap.$sp --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
  done 
done


# Calculate enrichment for autosomes (as in Yoo et al)
for hap in "pri" "alt" 
do 
  echo "Species Non-B Enrichment" | sed 's/ /\t/g' >new_sequence/${hap}_autosomal_enrichment.tsv 
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "borang" |grep -v siamang | while read -r sp latin filename;
  do
    cat new_sequence/${sp}_$hap/summary.txt |grep "autosomes" |awk -v OFS="\t" -v sp=$sp '{fold=$5/$8; print sp,$2,fold}' >>new_sequence/${hap}_autosomal_enrichment.tsv  
  done 
done


# Make a table with base pairs in and outside for statistic tests
for hap in "pri" "alt" 
do 
   cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt | grep -v "borang" |grep -v siamang | while read -r sp latin filename;
  do
    echo "Chr nonB Region bp_outside_nonB bp_inside_nonB" | sed 's/ /\t/' >new_sequence/${sp}_$hap/stat_table.tsv
    awk '(NR>1){out_new=$3-$4; out_old=$6-$7; print $1,$2,"New",out_new,$4,"\n",$1,$2,"Old",out_old,$7}' new_sequence/${sp}_$hap/summary.txt |sed 's/^ //' |sed 's/ /\t/g' >>new_sequence/${sp}_$hap/stat_table.tsv
  done 
done 

# Subsample the "new and old" (take ~half the regions from each file?) and recalculate statistics 
for hap in "alt" #"pri" #"alt" 
do 
   cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v human |grep -v "borang" |grep -v siamang | while read -r sp latin filename;
  do
    echo "looking at $sp, $hap"
    for chr in  "chrX" "chrY" #"autosomes"
    do
      unum=`awk '{}END{num=int(NR/2); print num}' new_sequence/${sp}_${hap}/$chr.unaligned.bed`
      anum=`awk '{}END{num=int(NR/2); print num}' new_sequence/${sp}_${hap}/$chr.aligned.bed`
      # Subsample half the sequences 10 times 
      for i in {1..10}
      do
        shuf -n $unum new_sequence/${sp}_${hap}/$chr.unaligned.bed |sort -k1,1 -k2,2n >new_sequence/${sp}_${hap}/resample.$i.$chr.unaligned.bed
        shuf -n $anum new_sequence/${sp}_${hap}/$chr.aligned.bed |sort -k1,1 -k2,2n >new_sequence/${sp}_${hap}/resample.$i.$chr.aligned.bed 
      done
    done
  done 
done 

# Calculate enrichment for each of these 
for hap in "pri" "alt" 
do 
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "borang" |grep -v siamang | while read -r sp latin filename;
  do
    for i in {1..10}
    do 
      echo '#!/bin/bash
      module load bedtools/2.31.0
     # echo "#Chr NonB NewTotBp NewNonB NewDens OldTotBp OldNonB OldDens" >new_sequence/'$sp'_'$hap'/resample.'$i'.summary.txt
      for chr in "autosomes" "chrX" "chrY"
      do
        for nb in "TRI" #"APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z" "all"
        do
          new=`intersectBed -a new_sequence/'$sp'_'$hap'/resample.'$i'.$chr.unaligned.bed -b nonB_annotation/'$sp'_'$hap'/${chr}_${nb}.bed  -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk -v nb=$nb '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print nb,sum_l,sum_nb,d}'"'"'`
          old=`intersectBed -a new_sequence/'$sp'_'$hap'/resample.'$i'.$chr.aligned.bed -b nonB_annotation/'$sp'_'$hap'/${chr}_${nb}.bed -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print sum_l,sum_nb,d}'"'"'`
          echo $chr" "$new" "$old >>new_sequence/'$sp'_'$hap'/resample.'$i'.summary.txt
        done 
      done 
      ' |sbatch -J $hap.$sp --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
    done
  done 
done

# Table with base pairs from the subsampled regions for the statistics test 
for hap in "pri" "alt" 
do 
   cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "borang" |grep -v siamang | while read -r sp latin filename;
  do
    echo "Resample Chr nonB Region bp_outside_nonB bp_inside_nonB" | sed 's/ /\t/' >new_sequence/${sp}_$hap/resample.stat_table.tsv
    for i in {1..10}
    do 
      awk -v i=$i '(NR>1){out_new=$3-$4; out_old=$6-$7; print i,$1,$2,"New",out_new,$4,"\n",i,$1,$2,"Old",out_old,$7}' new_sequence/${sp}_$hap/resample.$i.summary.txt |sed 's/^ //' |sed 's/ /\t/g' >>new_sequence/${sp}_$hap/resample.stat_table.tsv
    done
  done 
done 

# Numbers for Table 1
echo "Chr APR DR IR MR TRI G4 STR Z all" |sed 's/ /\t/g'
for chr in "autosomes" "chrX" "chrY"
do
  tmp="$chr"
  for nb in "APR" "DR"  "STR" "IR" "MR" "TRI" "GQ" "Z" "all"
  do
    fold=`grep $chr new_sequence/human_pri/summary.txt |grep $nb |awk '{fold=$5/$8; print fold}'`
    tmp="$tmp $fold"
  done 
  echo $tmp |sed 's/ /\t/g'
done 

# Stats for Table 1 is calculated in R, using the script    
# T2T_primate_nonB/R/chisquare_test_for_table1.R
   