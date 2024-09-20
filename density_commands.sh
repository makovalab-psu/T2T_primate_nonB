
################################## DENSITIES ###################################
# written by Linnéa Smeds, August 2024
# Summarizes the non-B DNA content along the genome in 100kb windows

# Create genomic windows and make length file
module load bedtools/2.31.0
mkdir -p ref/windows
rm -f T2T_primate_nonB/helpfiles/all_lengths.txt
for hap in "pri" "alt" 
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |while read -r trivial latin filename;
  do
    mkdir -p ref/windows/${trivial}_$hap
    bedtools makewindows -g ref/$filename.fai -w 100000 >ref/windows/${trivial}_$hap/all_100k_windows.bed;
    grep "chrX" ref/windows/${trivial}_$hap/all_100k_windows.bed >ref/windows/${trivial}_$hap/chrX_100k_windows.bed;
    grep "chrY" ref/windows/${trivial}_$hap/all_100k_windows.bed >ref/windows/${trivial}_$hap/chrY_100k_windows.bed;
    grep -v "chrY" ref/windows/${trivial}_$hap/all_100k_windows.bed |grep -v "chrX" |grep -v "chrM" >ref/windows/${trivial}_$hap/autosomes_100k_windows.bed;
  cut -f1,2 ref/$filename.fai |awk -v sp=${trivial}_$hap '{print sp"\t"$1"\t"$2}' >>T2T_primate_nonB/helpfiles/all_lengths.txt
  done
done 

# Fill windows with number of covered bases 
# (run as one batch job per hap/sp/chrtype/nonB)
mkdir densities
for hap in  "alt" #"pri"
do
  for chr in "chrX" "chrY" #"autosomes" 
  do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |while read -r trivial latin filename;
    do
      for n in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
      do      
        # One bed per sp/chrtype/nonBtype
        echo '#!/bin/bash     
        rm -f densities/'${trivial}'_'${hap}'/'${chr}'_'${n}'_100kb.bed
        intersectBed -wao -a ref/windows/'${trivial}'_'${hap}'/'${chr}'_100k_windows.bed -b nonB_annotation/'${trivial}'_'${hap}'/'${chr}'_'${n}'.bed | cut -f1,2,3,7 |awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" >densities/'${trivial}'_'${hap}'/'${chr}'_'${n}'_100kb.bed 
        '| sbatch -J $trivial.X --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G
      done 
    done 
  done 
done 
# Note that chrX and chrY are empty for the alternative haplotype. 

# Also combine into one file per species/haplo 
for hap in "pri" "alt" 
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |while read -r trivial latin filename;
  do
    rm densities/${trivial}_${hap}_comb_100kb.bed
    for chr in "autosomes" "chrX" "chrY" #
    do
      for n in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
      do
        awk -v OFS="\t" -v nb=$n '{print $0,nb}' densities/${trivial}_${hap}/${chr}_${n}_100kb.bed >>densities/${trivial}_${hap}_comb_100kb.bed
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
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "siamang" |while read -r trivial latin filename;
  do
    echo $trivial_${hap}
    cut -f1 ref/$filename.fai |sort -k1,1 |join - <(sort -k1,1 -k2,2n centromeres/$trivial/cen.bed) \
    | awk -v OFS="\t" -v s=$trivial '{print s,$1,$2,$3}' \
    |awk -v OFS="\t" '{if(NR==1){sp=$1;chr=$2;st=$3;e=$4}else{if($2==chr){e=$4}else{print sp,chr,st,e; sp=$1;chr=$2;st=$3;e=$4}}}END{print sp,chr,st,e}' >>densities/${hap}_6sp_centromeres.txt
    awk -v OFS="\t" -v s=$trivial '{print s,$0}'  densities/${trivial}_${hap}_comb_100kb.bed >>densities/${hap}_6sp_merged.txt
  done
done 

# Plot primary haplotypes with 
Rscript T2T_primate_nonB/R/plot_fig3_density_merged.R
# And secondary haploype with
Rscript T2T_primate_nonB/R/plot_figSX_density_alt_merged.R


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For later enrichment analysis, we need genomewide density of non-B
# (also make one with only autosomes)
for hap in "pri" "alt" #
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep siamang |while read -r trivial latin filename;
    do
        rm -f densities/${trivial}_${hap}_nonB_genome_wide.txt
        rm -f densities/${trivial}_${hap}_nonB_autosomes.txt
        totlen=`cat ref/$filename.fai | awk '{sum+=$2}END{print sum}'`
        autolen=`grep -v "chrX" ref/$filename.fai |grep -v "chrY" |grep -v "chrM"| awk '{sum+=$2}END{print sum}'`
        for non_b in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
        do
            cat nonB_annotation/${trivial}_$hap/autosomes_${non_b}.bed \
            | awk -v n=$non_b -v tot=$autolen '{sum+=$3-$2}END{d=sum/tot; print n, sum, d}' >>densities/${trivial}_${hap}_nonB_autosomes.txt
            cat nonB_annotation/${trivial}_$hap/autosomes_${non_b}.bed nonB_annotation/${trivial}_$hap/chrX_${non_b}.bed nonB_annotation/${trivial}_$hap/chrY_${non_b}.bed \
            | awk -v n=$non_b -v tot=$totlen '{sum+=$3-$2}END{d=sum/tot; print n, sum, d}' >>densities/${trivial}_${hap}_nonB_genome_wide.txt
        done
    done
done

# Per chromosome:
for hap in "pri" "alt"
do
    cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt  |while read -r trivial latin filename;
    do
      rm -f densities/${trivial}_${hap}_nonB_per_chrom.txt
      cat ref/$filename.fai |cut -f1,2 |grep -v "chrM" |while read -r chr len;
      do
          for non_b in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
          do
              cat nonB_annotation/${trivial}_$hap/*_${non_b}.bed | awk -v n=$non_b -v tot=$len -v c=$chr '($1==c){sum+=$3-$2}END{d=sum/tot; print c, n, sum, d}' >>densities/${trivial}_${hap}_nonB_per_chrom.txt
        done
      done
  done
done

