########################## ENRICHMENT IN REPEATS ###############################
# written by Linnéa Smeds, August 2024
# Code for analysing enrichment of non-B DNA motifs in repetitive sequences. 
# This analysis is done on the primary haplotype only.

# Dowload Repeat masker output and curated repeat files
mkdir ref/repeats/
cd ref/repeats/ 
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_composite-repeats_2022DEC.bed
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_new-satellites_2022DEC.bed
wget https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out
wget https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out
wget https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/repeats/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out
wget https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out
wget https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/repeats/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out
wget https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/repeats/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out
cd ../..

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Make directories for all species
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
  mkdir -p repeats/$sp
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert the downloaded files to bed (adapt the names slightly)
# Human (also composite and new repeats)
awk -v OFS="\t" '{if(/Satellite/){name=$11"_"$10}else{name=$11}; s=$6-1; print $5,s,$7,name}' ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |sort -k1,1 -k2,2n >repeats/human/RepeatMasker.bed
cp ref/repeats/chm13v2.0_new-satellites_2022DEC.bed repeats/human/new_satellites.bed
cp ref/repeats/chm13v2.0_composite-repeats_2022DEC.bed repeats/human/composite_repeats.bed
# Chimp
awk -v OFS="\t" '{if(/Satellite/){name=$11"_"$10}else{name=$11}; s=$6-1; print $5,s,$7,name}' ref/repeats/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out |sort -k1,1 -k2,2n >repeats/chimp/RepeatMasker.bed
# Bonobo
awk -v OFS="\t" '{if(/Satellite/){name=$11"_"$10}else{name=$11}; s=$6-1; print $5,s,$7,name}' ref/repeats/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out |sort -k1,1 -k2,2n >repeats/bonobo/RepeatMasker.bed
# Gorilla
awk -v OFS="\t" '{if(/Satellite/){name=$11"_"$10}else{name=$11}; s=$6-1; print $5,s,$7,name}' ref/repeats/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out |sort -k1,1 -k2,2n >repeats/gorilla/RepeatMasker.bed
# B_orangutan
awk -v OFS="\t" '{if(/Satellite/){name=$11"_"$10}else{name=$11}; s=$6-1; print $5,s,$7,name}' ref/repeats/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out |sort -k1,1 -k2,2n >repeats/borang/RepeatMasker.bed
# S_orangutan
awk -v OFS="\t" '{if(/Satellite/){name=$11"_"$10}else{name=$11}; s=$6-1; print $5,s,$7,name}' ref/repeats/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out |sort -k1,1 -k2,2n >repeats/sorang/RepeatMasker.bed
# Siamang - note, this repeatmasker is run on the older siamang, but the only
# difference is chr12 and chr19 that are swapped => so swap them manually!
awk -v OFS="\t" '{if(/Satellite/){name=$11"_"$10}else{name=$11}; s=$6-1; print $5,s,$7,name}' ref/repeats/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out |awk -v OFS="\t" '{if($1=="chr19_hap1"){$1="chr12_hap1"}else if($1=="chr12_hap1"){$1="chr19_hap1"}; print}' |sort -k1,1 -k2,2n >repeats/siamang/RepeatMasker.bed


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge repeats with non-B DNA
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
  mkdir -p repeats/$sp/overlap
  for chr in "autosomes" "chrX" "chrY"
  do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    for non_b in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
    do
        intersectBed -a repeats/'$sp'/RepeatMasker.bed -b nonB_annotation/'${sp}'_pri/'${chr}'_${non_b}.bed >repeats/'$sp'/overlap/'${chr}'_${non_b}_repmask.bed
        if [[ "'$sp'" == "human" ]]
        then
          echo "only for '$sp'"
          intersectBed -a repeats/'$sp'/new_satellites.bed -b nonB_annotation/'${sp}'_pri/'${chr}'_${non_b}.bed >repeats/'$sp'/overlap/'${chr}'_${non_b}_newsat.bed
          intersectBed -a repeats/'$sp'/composite_repeats.bed -b nonB_annotation/'${sp}'_pri/'${chr}'_${non_b}.bed >repeats/'$sp'/overlap/'${chr}'_${non_b}_compos.bed
        fi
    done
    '| sbatch -J $chr.$sp --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
  done
done

# AND FULL GENOME COMBINED!
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |grep "gorilla" |while read -r sp latin filename;
do
  mkdir -p repeats/$sp/overlap
  echo '#!/bin/bash
  module load bedtools/2.31.0
  for non_b in "APR" "DR" "GQ" "IR" "MR" "STR" "Z" "all"
  do
      intersectBed -a repeats/'$sp'/RepeatMasker.bed -b <(cat nonB_annotation/'${sp}'_pri/autosomes_${non_b}.bed nonB_annotation/'${sp}'_pri/chrX_${non_b}.bed nonB_annotation/'${sp}'_pri/chrY_${non_b}.bed) >repeats/'${sp}'/overlap/genome_${non_b}_repmask.bed
      if [[ "'$sp'" == "human" ]]
      then
        echo "only for '$sp'"
        intersectBed -a repeats/'$sp'/new_satellites.bed -b <(cat nonB_annotation/'${sp}'_pri/autosomes_${non_b}.bed nonB_annotation/'${sp}'_pri/chrX_${non_b}.bed nonB_annotation/'${sp}'_pri/chrY_${non_b}.bed) >repeats/'$sp'/overlap/genome_${non_b}_newsat.bed
        intersectBed -a repeats/'$sp'/composite_repeats.bed -b <(cat nonB_annotation/'${sp}'_pri/autosomes_${non_b}.bed nonB_annotation/'${sp}'_pri/chrX_${non_b}.bed nonB_annotation/'${sp}'_pri/chrY_${non_b}.bed) >repeats'/$sp'/overlap/genome_${non_b}_compos.bed
      fi
  done
  '| sbatch -J $chr.$sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=1:00:00 --partition=open
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a list with all repeat types and their total lengths for each chr type
module load python/3.11.2
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
  grep "chrX" repeats/$sp/RepeatMasker.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/chrX.repeat_lengths.txt
  grep "chrY" repeats/$sp/RepeatMasker.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/chrY.repeat_lengths.txt
  grep -v "chrX" repeats/$sp/RepeatMasker.bed|grep -v "chrY" |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/autosomes.repeat_lengths.txt
  cat repeats/$sp/RepeatMasker.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/genome_repeat_lengths.txt
done

# And for human satellites/composite
sp="human"
grep "chrX" repeats/$sp/composite_repeats.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/chrX.compos_lengths.txt
grep "chrY" repeats/$sp/composite_repeats.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/chrY.compos_lengths.txt
grep -v "chrX" repeats/$sp/composite_repeats.bed|grep -v "chrY" |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/autosomes.compos_lengths.txt
cat repeats/$sp/composite_repeats.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/genome_compos_lengths.txt
grep "chrX" repeats/$sp/new_satellites.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/chrX.newsat_lengths.txt
grep "chrY" repeats/$sp/new_satellites.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/chrY.newsat_lengths.txt
grep -v "chrX" repeats/$sp/new_satellites.bed|grep -v "chrY" |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/autosomes.newsat_lengths.txt
cat repeats/$sp/new_satellites.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/genome_newsat_lengths.txt



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summarize the enrichment of each non-B motif in all repeat classes
# (Print NA if there are no such repeats in the genome)
cat T2T_primate_nonB/helpfiles/pri_species_list.txt  |grep "human" |while read -r sp latin filename;
do
  for chr in  "autosomes" "chrX" "chrY" 
  do
    echo "Repeat APR DR GQ IR MR STR Z" |sed 's/ /\t/g' >repeats/${sp}_${chr}_repeat_enrichment.tsv
    cat repeats/$sp/${chr}.repeat_lengths.txt | while read -r rep replen;
    do
      #repname=`echo $rep |sed 's/\//_/g'`
      tmp=`echo $rep`
      #Go through non-B
      cat densities/${sp}_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
      do
        d=`cat repeats/$sp/overlap/${chr}_${non_b}_repmask.bed | awk -v r=$rep -v l=$replen -v dtot=$dens '($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print frac}}'`
        tmp=`echo $tmp" "$d`
        echo $tmp >tmp
      done
      cat tmp |sed 's/ /\t/g' >>repeats/${sp}_${chr}_repeat_enrichment.tsv
    done
  done
done

# FULL GENOME 
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r trivial latin filename;
do
  echo '#!/bin/bash
  echo "Repeat APR DR GQ IR MR STR Z" |sed "s/ /\t/g" >repeats/'$trivial'_genome_repeat_enrichment.tsv
  cat repeats/'$trivial'/genome_repeat_lengths.txt | while read -r rep replen;
  do
    tmp=`echo $rep`
    #Go through non-B
    cat densities/'$trivial'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
    do
      d=`cat repeats/'$trivial'/overlap/genome_${non_b}_repmask.bed | awk -v r=$rep -v l=$replen -v dtot=$dens '"'"'($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print frac}}'"'"'`
      tmp=`echo $tmp" "$d`
      echo $tmp >tmp
    done
    cat tmp |sed "s/ /\t/g" >>repeats/'$trivial'_genome_repeat_enrichment.tsv
  done
  '| sbatch -J $trivial --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done

# And for human satellites / composite repeats
# Per chrom type # CHECK THIS PART, UNFINISHED OUTPUT FILES 
trivial="human"
for type in "compos" "newsat"
do
  echo '#!/bin/bash
  for chr in "chrY" "autosomes" "chrX"
  do
    echo "Repeat APR DR GQ IR MR STR Z" |sed "s/ /\t/g" >repeats/'$trivial'_${chr}_'$type'_enrichment.tsv
    cat repeats/'$trivial'/${chr}.'$type'_lengths.txt | while read -r rep replen;
    do
      tmp=`echo $rep`
      cat densities/'$trivial'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
      do
        d=`cat repeats/'$trivial'/overlap/${chr}_${non_b}_'$type'.bed | awk -v r=$rep -v l=$replen -v dtot=$dens '"'"'($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print frac}}'"'"'`
        tmp=`echo $tmp" "$d`
        echo $tmp >tmp
      done
      cat tmp |sed "s/ /\t/g" >>repeats/'$trivial'_${chr}_'$type'_enrichment.tsv
    done
  done
  '| sbatch -J $type --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done
# All chromosomes together 
for type in "compos" "newsat"
do
  echo '#!/bin/bash
  echo "Repeat APR DR GQ IR MR STR Z" |sed "s/ /\t/g" >repeats/'$trivial'_genome_'$type'_enrichment.tsv
  cat repeats/'$trivial'/genome_'$type'_lengths.txt | while read -r rep replen;
  do
    tmp=`echo $rep`
    cat densities/'$trivial'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
    do
      d=`cat repeats/'$trivial'/overlap/genome_${non_b}_'$type'.bed | awk -v r=$rep -v l=$replen -v dtot=$dens '"'"'($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print frac}}'"'"'`
      tmp=`echo $tmp" "$d`
      echo $tmp >tmp
    done
    cat tmp |sed "s/ /\t/g" >>repeats/'$trivial'_genome_'$type'_enrichment.tsv
  done
  '| sbatch -J $type --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Before plotting, shorten the longest names 
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r trivial latin filename;
do
 cat repeats/${trivial}_genome_repeat_enrichment.tsv |sed 's/-subunit_LSAU-BSAT_rnd-/*/' |sed 's/_family-/*/' |sed 's/-VAR_rnd-/*/' | sed 's/-chromosome_ALRY-/*/' >repeats/${trivial}_genome_repeat_enrichment_renamed.tsv
 cat repeats/${trivial}/genome_repeat_lengths.txt |sed 's/-subunit_LSAU-BSAT_rnd-/*/' |sed 's/_family-/*/' |sed 's/-VAR_rnd-/*/' |sed 's/-chromosome_ALRY-/*/' >repeats/${trivial}/genome_repeat_lengths_renamed.txt
done

# For non-humans, we want to merge the files
echo "Species Repeat APR DR GQ IR MR STR Z" |sed 's/ /\t/g' >repeats/6sp_genome_repeat_enrichment_renamed.tsv
echo "Species Repeat Len" |sed 's/ /\t/g' >repeats/6sp_genome_repeat_lengths_renamed.txt
cat species_list.txt |grep -v "human" |while read -r trivial latin filename;
do
   awk -v sp=$trivial '(NR>1){print sp,$0}' repeats/${trivial}_genome_repeat_enrichment_renamed.tsv |sed 's/^ //' |sed 's/ /\t/g' >>repeats/6sp_genome_repeat_enrichment_renamed.tsv
   awk -v sp=$trivial '{print sp,$0}' repeats/${trivial}/genome_repeat_lengths_renamed.txt |sed 's/^ //' |sed 's/ /\t/g' >>repeats/6sp_genome_repeat_lengths_renamed.txt
done

# And merge human composite repeats and new satellites 
echo "Type Repeat APR DR GQ IR MR STR Z" |sed 's/ /\t/g' >repeats/human_genome_compsat_enrichment.tsv
echo "Type Repeat Len" |sed 's/ /\t/g' >repeats/human_genome_compsat_lengths.txt
for type in "newsat" "compos"
do
   awk -v sp=$type '(NR>1){print sp,$0}' repeats/human_genome_${type}_enrichment.tsv |sed 's/^ //' |sed 's/ /\t/g' >>repeats/human_genome_compsat_enrichment.tsv
   awk -v sp=$type '{print sp,$0}' repeats/human/genome_${type}_lengths.txt |sed 's/^ //' |sed 's/ /\t/g' >>repeats/human_genome_compsat_lengths.txt
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot main figure (human) 
Rscript T2T_primate_nonB/R/plot_fig5_repeats.R 
# Plot supplementary figure (other species)
Rscript T2T_primate_nonB/R/plot_figS2_repeats.R
# Plot supplementary figure (new satellites and composite repeats)
Rscript T2T_primate_nonB/R/plot_figS3_repeats_compsat.R



