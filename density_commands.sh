
################################## DENSITIES ###################################
# written by Linnéa Smeds, August 2024
# Summarizes the non-B DNA content along the genome in 100kb windows

# Create genomic windows and make length file
module load bedtools/2.31.0
mkdir -p ref/windows
#rm -f T2T_primate_nonB/helpfiles/all_lengths.txt
for hap in  "alt" "pri"
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep siamang |while read -r sp latin filename;
  do
    mkdir -p ref/windows/${sp}_$hap
    bedtools makewindows -g ref/$filename.fai -w 100000 >ref/windows/${sp}_$hap/all_100k_windows.bed;
    grep "chrX" ref/windows/${sp}_$hap/all_100k_windows.bed >ref/windows/${sp}_$hap/chrX_100k_windows.bed;
    grep "chrY" ref/windows/${sp}_$hap/all_100k_windows.bed >ref/windows/${sp}_$hap/chrY_100k_windows.bed;
    grep -v "chrY" ref/windows/${sp}_$hap/all_100k_windows.bed |grep -v "chrX" |grep -v "chrM" >ref/windows/${sp}_$hap/autosomes_100k_windows.bed;
  cut -f1,2 ref/$filename.fai |awk -v sp=${sp}_$hap '{print sp"\t"$1"\t"$2}' >>T2T_primate_nonB/helpfiles/all_lengths.txt
  done
done 

# Try longer windows for testing correlations 
module load bedtools/2.31.0
for hap in "pri" "alt" 
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep human |while read -r sp latin filename;
  do
    bedtools makewindows -g ref/$filename.fai -w 1000000 >ref/windows/${sp}_$hap/all_1Mb_windows.bed;
    bedtools makewindows -g ref/$filename.fai -w 5000000 >ref/windows/${sp}_$hap/all_5Mb_windows.bed;
  done
done 

# Fill windows with number of covered bases 
# (run as one batch job per hap/sp/chrtype/nonB)
mkdir densities
for hap in  "alt" "pri"
do
  for chr in "chrX" "chrY" "autosomes" 
  do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |while read -r sp latin filename;
    do
      for n in "TRI" #"APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
      do      
        # One bed per sp/chrtype/nonBtype
        echo '#!/bin/bash  
          module load bedtools/2.31.0
        rm -f densities/'${sp}'_'${hap}'/'${chr}'_'${n}'_100kb.bed
        intersectBed -wao -a ref/windows/'${sp}'_'${hap}'/'${chr}'_100k_windows.bed -b nonB_annotation/'${sp}'_'${hap}'/'${chr}'_'${n}'.bed | cut -f1,2,3,7 |awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" >densities/'${sp}'_'${hap}'/'${chr}'_'${n}'_100kb.bed 
        '| sbatch -J $sp.X --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G
      done 
    done 
  done 
done 
# Note that chrX and chrY are empty for the alternative haplotype. 

# Also combine into one file per species/haplo 
for hap in "alt" # "pri" 
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep siamang |while read -r sp latin filename;
  do
    rm densities/${sp}_${hap}_comb_100kb.bed
    for chr in "autosomes" "chrX" "chrY" #
    do
      for n in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
      do
        awk -v OFS="\t" -v nb=$n '{print $0,nb}' densities/${sp}_${hap}/${chr}_${n}_100kb.bed >>densities/${sp}_${hap}_comb_100kb.bed
      done
    done 
  done 
done
       
# For plotting 6 species into one single plot, merge the data:
for hap in "pri" "alt" #
do
  rm densities/${hap}_6sp_lengths.txt
  rm densities/${hap}_6sp_centromeres.txt
  rm densities/${hap}_6sp_merged.txt
  grep "_$hap" T2T_primate_nonB/helpfiles/all_lengths.txt |sed "s/_$hap//" |grep -v "siamang" \
  |grep -v "random" |awk -v OFS="\t" '{split($2,s,"_"); print $0,s[1]}' >>densities/${hap}_6sp_lengths.txt 
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "siamang" |while read -r sp latin filename;
  do
    echo $sp_${hap}
    cut -f1 ref/$filename.fai |sort -k1,1 |join - <(sort -k1,1 -k2,2n centromeres/$sp/cen.bed) \
    | awk -v OFS="\t" -v s=$sp '{print s,$1,$2,$3}' \
    |awk -v OFS="\t" '{if(NR==1){sp=$1;chr=$2;st=$3;e=$4}else{if($2==chr){e=$4}else{print sp,chr,st,e; sp=$1;chr=$2;st=$3;e=$4}}}END{print sp,chr,st,e}' >>densities/${hap}_6sp_centromeres.txt
    awk -v OFS="\t" -v s=$sp '{print s,$0}'  densities/${sp}_${hap}_comb_100kb.bed >>densities/${hap}_6sp_merged.txt
  done
done 

# Plot haplotypes with circos, see circos_commands.sh an example code in circos/

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate GC content per window 

module load bedtools/2.31.0
for hap in "pri" "alt" 
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep "human" |while read -r sp latin filename;
  do
    bedtools nuc -fi ref/assemblies/$filename -bed ref/windows/${sp}_$hap/all_100k_windows.bed >ref/windows/${sp}_$hap/all_100k_windows.withGC.bed
    bedtools nuc -fi ref/assemblies/$filename -bed ref/windows/${sp}_$hap/all_1Mb_windows.bed >ref/windows/${sp}_$hap/all_1Mb_windows.withGC.bed
    bedtools nuc -fi ref/assemblies/$filename -bed ref/windows/${sp}_$hap/all_5Mb_windows.bed >ref/windows/${sp}_$hap/all_5Mb_windows.withGC.bed
  done 
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Instead of intersecting the new window sizes with nonB and calculate density
# again, I just intersect the 1Mb window file with the 100kb window file and 
# sum up the nonB count for each window 

for hap in "pri" "alt" 
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep "human" |while read -r sp latin filename;
  do
    for size in "1Mb" "5Mb"
    do
      echo "Chr Start Stop GCcont nonB_count nonB" |sed 's/ /\t/g' >densities/${sp}_${hap}_nonB_GC.$size.bed
      for n in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
        do
          intersectBed -a <(cut -f1,2,3,5 ref/windows/${sp}_$hap/all_${size}_windows.withGC.bed) -b <(cat densities/${sp}_$hap/autosomes_${n}_100kb.bed densities/human_pri/chrX_${n}_100kb.bed densities/${sp}_$hap/chrY_${n}_100kb.bed) -wo |awk -v OFS="\t" -v nb=$n '{if(NR==1){chr=$1; start=$2; end=$3; gc=$4; sum=$8}else{if(chr==$1 && start==$2){sum+=$8}else{print chr,start,end,gc,sum,nb; chr=$1; start=$2; end=$3; gc=$4; sum=$8}}}END{print chr,start,end,gc,sum,nb}' >densities/${sp}_$hap/genome_${n}_${size}.withGG.bed
          cat densities/${sp}_$hap/genome_${n}_${size}.withGG.bed >>densities/${sp}_${hap}_nonB_GC.$size.bed
       done 
    done 
  done 
done 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For later enrichment analysis, we need genomewide density of non-B
# (also make files with autosomes,chrX and chrY separately)
for hap in "pri" "alt" #
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v human |while read -r sp latin filename;
  do
    echo ${sp}"_"$hap
    rm -f densities/${sp}_${hap}_nonB_genome_wide.txt
    rm -f densities/${sp}_${hap}_nonB_autosomes.txt
    rm -f densities/${sp}_${hap}_nonB_chrX.txt
    rm -f densities/${sp}_${hap}_nonB_chrY.txt
    totlen=`cat ref/$filename.fai | awk '{sum+=$2}END{print sum}'`
    autolen=`grep -v "chrX" ref/$filename.fai |grep -v "chrY" |grep -v "chrM"| awk '{sum+=$2}END{print sum}'`
    xlen=`grep "chrX" ref/$filename.fai | awk '{sum+=$2}END{print sum}'`
    ylen=`grep "chrY" ref/$filename.fai | awk '{sum+=$2}END{print sum}'`
    for non_b in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
    do
      cat nonB_annotation/${sp}_$hap/autosomes_${non_b}.bed \
      | awk -v n=$non_b -v tot=$autolen '{sum+=$3-$2}END{d=sum/tot; print n, sum, d}' >>densities/${sp}_${hap}_nonB_autosomes.txt
      cat nonB_annotation/${sp}_$hap/chrX_${non_b}.bed \
      | awk -v n=$non_b -v tot=$xlen '{sum+=$3-$2}END{d=sum/tot; print n, sum, d}' >>densities/${sp}_${hap}_nonB_chrX.txt
      cat nonB_annotation/${sp}_$hap/chrY_${non_b}.bed \
      | awk -v n=$non_b -v tot=$ylen '{sum+=$3-$2}END{d=sum/tot; print n, sum, d}' >>densities/${sp}_${hap}_nonB_chrY.txt
      cat nonB_annotation/${sp}_$hap/autosomes_${non_b}.bed nonB_annotation/${sp}_$hap/chrX_${non_b}.bed nonB_annotation/${sp}_$hap/chrY_${non_b}.bed \
      | awk -v n=$non_b -v tot=$totlen '{sum+=$3-$2}END{d=sum/tot; print n, sum, d}' >>densities/${sp}_${hap}_nonB_genome_wide.txt
    done
  done
done

# Per chromosome:
for hap in "pri" "alt"
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt  |while read -r sp latin filename;
    do
      rm -f densities/${sp}_${hap}_nonB_per_chrom.txt
      cat ref/$filename.fai |cut -f1,2 |grep -v "chrM" |while read -r chr len;
      do
          for non_b in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
          do
              cat nonB_annotation/${sp}_$hap/*_${non_b}.bed | awk -v n=$non_b -v tot=$len -v c=$chr '($1==c){sum+=$3-$2}END{d=sum/tot; print c, n, sum, d}' >>densities/${sp}_${hap}_nonB_per_chrom.txt
        done
      done
  done
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GC Corrected 'densities' (Divide by GC-bp instead of full length)
# THIS PART WAS NOT USED FOR THE FINAL PAPER. 

# Find GC and AT of full genomes 
for hap in "pri" "alt" #
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep siamang |while read -r sp latin filename;
  do
    echo '#!/bin/bash
    cat ref/assemblies/'$filename' | awk '"'"'(!/^>/){gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"")}END{print gc"\t"at}'"'"' >ref/assemblies/'$filename'.nuc
  '| sbatch -J $ind --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
  done 
done 

# GC CORRECTED 
module load bedtools/2.31.0
for hap in "pri" "alt" #
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep human |while read -r sp latin filename;
  do
    echo ${sp}"_"$hap
    rm -f densities/${sp}_${hap}_GCcorr_nonB_genome_wide.txt
    totGC=`cut -f1 ref/assemblies/$filename.nuc`
    echo $totGC
    for non_b in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
    do
      cat nonB_annotation/${sp}_$hap/autosomes_${non_b}.bed nonB_annotation/${sp}_$hap/chrX_${non_b}.bed nonB_annotation/${sp}_$hap/chrY_${non_b}.bed \
      | awk -v n=$non_b -v tot=$totGC '{sum+=$3-$2}END{d=sum/tot; print n, sum, d}' >>densities/${sp}_${hap}_GCcorr_nonB_genome_wide.txt
    done
  done
done

# AT CORRECTED 
module load bedtools/2.31.0
for hap in "pri" "alt" #
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep human |while read -r sp latin filename;
  do
    echo ${sp}"_"$hap
    rm -f densities/${sp}_${hap}_ATcorr_nonB_genome_wide.txt
    totAT=`cut -f2 ref/assemblies/$filename.nuc`
    echo $totAT
    for non_b in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
    do
      awk -v n=$non_b -v tot=$totAT '{sum+=$3-$2}END{d=sum/tot; print n, sum, d}' nonB_annotation/${sp}_$hap/genome_${non_b}.bed  >>densities/${sp}_${hap}_ATcorr_nonB_genome_wide.txt
    done
  done
done


