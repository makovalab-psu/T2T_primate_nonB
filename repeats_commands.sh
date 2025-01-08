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
  for chr in "genome" "autosomes" "chrX" "chrY"
  do
    echo '#!/bin/bash
    echo "Looking at '$sp' '$chr'"
    module load bedtools/2.31.0
    for non_b in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z" "all"
    do
        intersectBed -a repeats/'$sp'/RepeatMasker.bed -b nonB_annotation/'${sp}'_pri/'${chr}'_${non_b}.bed >repeats/'$sp'/overlap/'${chr}'_${non_b}_repmask.bed
        if [[ "'$sp'" == "human" ]]
        then
          echo "only for '$sp'"
          intersectBed -a repeats/'$sp'/new_satellites.bed -b nonB_annotation/'${sp}'_pri/'${chr}'_${non_b}.bed >repeats/'$sp'/overlap/'${chr}'_${non_b}_newsat.bed
          intersectBed -a repeats/'$sp'/composite_repeats.bed -b nonB_annotation/'${sp}'_pri/'${chr}'_${non_b}.bed >repeats/'$sp'/overlap/'${chr}'_${non_b}_compos.bed
        fi
    done
    '| sbatch -J $chr.$sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=1:00:00 --partition=open
  done
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a list with all repeat types and their total lengths
module load python/3.11.2
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
   cat repeats/$sp/RepeatMasker.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/genome_repeat_lengths.txt
done

# And for human satellites/composite
sp="human"
cat repeats/$sp/composite_repeats.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/genome_compos_lengths.txt
cat repeats/$sp/new_satellites.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/genome_newsat_lengths.txt
# Note:Composite repeat 'Charlie 5' was manually edited to Charlie 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summarize the enrichment of each non-B motif in all repeat classes
# (Print NA if there are no such repeats in the genome)
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
  echo '#!/bin/bash
  echo "Looking at repeats in '$sp'"
  echo "Repeat APR DR GQ IR MR TRI STR Z" |sed "s/ /\t/g" >repeats/'$sp'_genome_repeat_enrichment.tsv
  cat repeats/'$sp'/genome_repeat_lengths.txt | while read -r rep replen;
  do
    tmp=`echo $rep`
    #Go through non-B
    cat densities/'$sp'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
    do
      d=`cat repeats/'$sp'/overlap/genome_${non_b}_repmask.bed | awk -v r=$rep -v l=$replen -v dtot=$dens '"'"'($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print frac}}'"'"'`
      tmp=`echo $tmp" "$d`
      echo $tmp >tmp.'$sp'
    done
    cat tmp.'$sp' |sed "s/ /\t/g" >>repeats/'$sp'_genome_repeat_enrichment.tsv
  done
  '| sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=1:00:00 --partition=open
done

# COMPOSITE AND NEW SATELLITES, HUMAN
sp="human"
for type in "compos" "newsat"
do
  echo '#!/bin/bash
  echo "Looking at repeats in "'$sp' '$type'
  echo "Repeat APR DR GQ IR MR TRI STR Z" |sed "s/ /\t/g" >repeats/'$sp'_genome_'$type'_enrichment.tsv
  cat repeats/'$sp'/genome_'$type'_lengths.txt | while read -r rep replen;
  do
    tmp=`echo $rep`
    cat densities/'$sp'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
    do
      d=`cat repeats/'$sp'/overlap/genome_${non_b}_'$type'.bed | awk -v r=$rep -v l=$replen -v dtot=$dens '"'"'($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print frac}}'"'"'`
      tmp=`echo $tmp" "$d`
      echo $tmp >tmp.'$type'
    done
    cat tmp.'$type' |sed "s/ /\t/g" >>repeats/'$sp'_genome_'$type'_enrichment.tsv
  done
  '| sbatch -J $type --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Before plotting, shorten the longest names 
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
 cat repeats/${sp}_genome_repeat_enrichment.tsv |sed 's/-subunit_LSAU-BSAT_rnd-/*/' |sed 's/_family-/*/' |sed 's/-VAR_rnd-/*/' | sed 's/-chromosome_ALRY-/*/' >repeats/${sp}_genome_repeat_enrichment_renamed.tsv
 cat repeats/${sp}/genome_repeat_lengths.txt |sed 's/-subunit_LSAU-BSAT_rnd-/*/' |sed 's/_family-/*/' |sed 's/-VAR_rnd-/*/' |sed 's/-chromosome_ALRY-/*/' >repeats/${sp}/genome_repeat_lengths_renamed.txt
done

# For non-humans, we want to merge the files
echo "Species Repeat APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g' >repeats/6sp_genome_repeat_enrichment_renamed.tsv
echo "Species Repeat Len" |sed 's/ /\t/g' >repeats/6sp_genome_repeat_lengths_renamed.txt
cat species_list.txt |grep -v "human" |while read -r sp latin filename;
do
   awk -v sp=$sp '(NR>1){print sp,$0}' repeats/${sp}_genome_repeat_enrichment_renamed.tsv |sed 's/^ //' |sed 's/ /\t/g' >>repeats/6sp_genome_repeat_enrichment_renamed.tsv
   awk -v sp=$sp '{print sp,$0}' repeats/${sp}/genome_repeat_lengths_renamed.txt |sed 's/^ //' |sed 's/ /\t/g' >>repeats/6sp_genome_repeat_lengths_renamed.txt
done

# And merge human composite repeats and new satellites 
echo "Type Repeat APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g' >repeats/human_genome_compsat_enrichment.tsv
echo "Type Repeat Len" |sed 's/ /\t/g' >repeats/human_genome_compsat_lengths.txt
for type in "newsat" "compos"
do
   awk -v sp=$type '(NR>1){print sp,$0}' repeats/human_genome_${type}_enrichment.tsv |sed 's/^ //' |sed 's/ /\t/g' >>repeats/human_genome_compsat_enrichment.tsv
   awk -v sp=$type '{print sp,$0}' repeats/human/genome_${type}_lengths.txt |sed 's/^ //' |sed 's/ /\t/g' >>repeats/human_genome_compsat_lengths.txt
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot main figure (human, need to first generate methylation data)
Rscript R/plot_fig6ABC_repeats_methylation.R
# Plot supplementary figure (other species)
Rscript R/plot_figS9_enrichment_repeats.R
# Plot supplementary figure (new satellites and composite repeats)
Rscript R/plot_figS8_enrichment_compsat_human.R



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Repeating analysis of human WALUSAT repeat on chr14:260778-634253
# using G4Hunter online predictor (https://bioinformatics.ibp.cz/#/analyse/quadruplex)
# with default settings (window=25 and score=1.2).
# Results were downloaded as a csv and processed as follows. 

# Find number of predicted G4s
sed 's/"//g' repeats/human/6f162de4-0e77-4b52-86d9-d6b236c81392.csv |tail -n+2 |wc
#   22701  113505 1366880

# Number of different sequences 
sed 's/"//g' repeats/human/6f162de4-0e77-4b52-86d9-d6b236c81392.csv |tail -n+2 |cut -f5 |sort |uniq -c |wc
#    531    1062   16461
# With more than 50 occurences:
sed 's/"//g' repeats/human/6f162de4-0e77-4b52-86d9-d6b236c81392.csv |tail -n+2 |cut -f5 |sort |uniq -c |awk '($1>+50){print}' |wc
#     25      50     775

# With more than 1000 occurences:
sed 's/"//g' repeats/human/6f162de4-0e77-4b52-86d9-d6b236c81392.csv |tail -n+2 |cut -f5 |sort |uniq -c |awk '($1>1000){print}' 
#4898 AGGCGGGGTCAGAGGAATAGAAAGG
#4909 GCGGGGTCAGAGGAATAGAAAGGGA
#4902 GGCGGGGTCAGAGGAATAGAAAGGG
#2622 GGGGTCAGAGGAATAGAAAGGGATA
# (all these share a 21 bp motif)
# MOST COMMON SCORES OF THESE 
# AGGCGGGGTCAGAGGAATAGAAAGG 1.2000000476837158
# GCGGGGTCAGAGGAATAGAAAGGGA 1.2400000095367432
# GGCGGGGTCAGAGGAATAGAAAGGG 1.3200000524520874
# GGGGTCAGAGGAATAGAAAGGGATA 1.2000000476837158

# Split string every 64bp and count occcurences:
grep -v ">" repeats/human/walusat.chr14.fa |tr "\n" " " |sed 's/ //g' | fold -w 64 |sort |uniq -c |sort -k1,1nr |less
# This doesn't work as there are some indels that make the pattern shift.
# Instead, split string on a common pattern, and then split on 64mer 
# (since some units have mismatches in the pattern:
grep -v ">" repeats/human/walusat.chr14.fa |tr "\n" " " |sed 's/ //g' | sed 's/ggggtca/\nggggtca/g'|awk '{for(i=1; i<=length; i+=64) print substr($0,i,64)}'|sort |uniq -c |sort -k1,1nr |head 
#1186 ggggtcagaggaatagaaagggacagggctgaagaacacaggtcgctgcatttagaaaggaggc
# 779 ggggtcagaggaatagaaagggatagggctgaagaacacaggtcgctgcatttagaaaggaggc
# 478 ggggtcagaggaatagaaagggacagggctgaagaacagaggtcgctgcatttagaaaggaggc
# 442 ggggtcagaggaatagaaagggatagggctgaagaacagaggtcgctgcatttagaaaggaggc
# 265 ggggtcagaggaatagaaagggataggcctgaagaacacaggtcgctgcatttagaaaggaggc
# 214 ggggtcagaggaatagaaagggataggcctgaagaacagaggtcgctgcatttagaaaggaggc
# 168 ggggtcagaggaatagaaagggataggactgaagaacacaggtcgctgcatttagaaaggaggc
#  97 ggggtcagaggaatagaaagggatagggatgaagaacacaggtcgctgcatttagaaaggaggc
#  86 ggggttagaggaatagaaagggacagggctgaagaacacaggtcgctgcatttagaaaggaggc
#  85 ggggtcagaggaatagaaagggatggggctgaagaacagaggtcgctgcatttagaaaggaggc
# Number of 64 mers in total:
grep -v ">" repeats/human/walusat.chr14.fa |tr "\n" " " |sed 's/ //g' | sed 's/ggggtca/\nggggtca/g'|awk '{for(i=1; i<=length; i+=64) print substr($0,i,64)}'|wc
#    5849    5849  379325



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A DEEPER LOOK INTO SST1 
# SST1 is thought to be the breakpoint in the formation of Robertsonian 
# chromosomes by chr 13, 14 and 21 in human. In our general analysis, SST1 is 
# enriched for G4s but no other non-B motifs. In the recent paper by Gomes deLima
# et al, they show that the SST1 repeat array is rather different on 13/14/21 
# compared to other chromosomes (chr4/17/19 is clustered together as a different 
# type, and chrY has a third type). I want to see if the non-B DNA content also
# differs between the chromosomes. 

# check space in between
grep SST1 ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |cut -f5-7 |sort -k1,1 -k2,2n |awk '{if(NR==1){chr=$1; end=$3}else{if(chr==$1){diff=$2-end; print $0"\t"diff; end=$3}else{chr=$1; end=$3}}}' |less

# Look at occurences in the array on chr13,
# both in total..
awk '($1=="chr13" && $2>12301367 && $3<12440010){print}' nonB_annotation/human_pri/autosomes_GQ.bed |wc
     99     297    2376
awk '($1=="chr13" && $2>12301367 && $3<12440010){print}' nonB_annotation/human_pri/autosomes_Z.bed |wc
    194     582    4656
# ..and only in annotated repeats
awk '($1=="chr13" && $2>12301367 && $3<12440010){print}' repeats/human/overlap/autosomes_GQ_repmask.bed |wc
     99     396    4467
awk '($1=="chr13" && $2>12301367 && $3<12440010){print}' repeats/human/overlap/autosomes_Z_repmask.bed |grep "SST1" |wc
      6      24     270
# This means there are almost no Z-DNA motifs in SST1, but many in the space between 
# the repeat units


# Make files with SST1 on each of the chromosomes of interest
mkdir repeats/human/SST1
for c in "chr13" "chr14" "chr21" "chr4" "chr17" "chr19" "chrY"
do
  grep SST1 ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |cut -f5,6,7 |awk -v c=$c '($1==c){print}' |sort -k2,2n >repeats/human/SST1/$c.exactRep.bed
  grep SST1 ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |cut -f5,6,7 |awk -v c=$c '($1==c){print}' |sort -k2,2n | awk -v OFS="\t" '{if(NR==1){start=$2; end=$3}else{if($2-end<1000){end=$3}else{print $1,start,end; start=$2;end=$3}}}END{print $1,start,end}' >repeats/human/SST1/$c.merge1kb.bed
  grep SST1 ref/repeats/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out |cut -f5,6,7 |awk -v c=$c '($1==c){print}' |sort -k2,2n | awk -v OFS="\t" '{if(NR>1){if($2-end<1000){print $1,end,$2}}; end=$3}' >repeats/human/SST1/$c.intermediate1kb.bed
done

# Go through the bed files and calculate density and fold enrichment compared 
# to genomewide density
for type in "intermediate1kb" "exactRep" #"merge1kb"
do
  echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Chr NonB Sum Density FoldEnr" |sed "s/ /\t/g" >repeats/human/SST1/'$type'.enrichment.tsv
    for c in "chr13" "chr14" "chr21" "chr4" "chr17" "chr19" "chrY"
    do

      len=`awk -v sum=0 '"'"'{sum+=$3-$2}END{print sum}'"'"' repeats/human/SST1/$c.'$type'.bed`
      echo $len
      cat densities/human_pri_nonB_genome_wide.txt |grep -v "all" | while read -r nb tot dens;
      do
        intersectBed -a repeats/human/SST1/$c.'$type'.bed -b nonB_annotation/human_pri/genome_$nb.bed -wo |awk -v l=$len -v gwd=$dens -v nb=$nb -v c=$c -v OFS="\t" '"'"'{sum+=$7}END{dens=sum/l; fold=dens/gwd; print c,nb,sum,dens,fold}'"'"' >>repeats/human/SST1/'$type'.enrichment.tsv
      done 
    done 
  ' | sbatch -J $type --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done 

# G4 density is higher in all SST (both acrocentric and other), while Z density
# in 'spacers' is higher compared to genome-wide in all except on chrY. 

# Make three groups 
for type in "intermediate1kb" #"exactRep" "merge1kb"
do
  cat repeats/human/SST1/chr13.$type.bed repeats/human/SST1/chr14.$type.bed repeats/human/SST1/chr21.$type.bed >repeats/human/SST1/sf1.$type.bed
  cat repeats/human/SST1/chr4.$type.bed repeats/human/SST1/chr17.$type.bed repeats/human/SST1/chr19.$type.bed >repeats/human/SST1/sf2.$type.bed
  cat repeats/human/SST1/chrY.$type.bed >repeats/human/SST1/sf3.$type.bed
done

# Repeat Fold enrichment analysis 
for type in "intermediate1kb" "exactRep" #"merge1kb"
do
  echo '#!/bin/bash
    module load bedtools/2.31.0
    echo "Chr NonB Sum Density FoldEnr" |sed "s/ /\t/g" >repeats/human/SST1/'$type'.sf.enrichment.tsv
    for c in "sf1" "sf2" "sf3"
    do

      len=`awk -v sum=0 '"'"'{sum+=$3-$2}END{print sum}'"'"' repeats/human/SST1/$c.'$type'.bed`
      echo $len
      cat densities/human_pri_nonB_genome_wide.txt |grep -v "all" | while read -r nb tot dens;
      do
        intersectBed -a repeats/human/SST1/$c.'$type'.bed -b nonB_annotation/human_pri/genome_$nb.bed -wo |awk -v l=$len -v gwd=$dens -v nb=$nb -v c=$c -v OFS="\t" '"'"'{sum+=$7}END{dens=sum/l; fold=dens/gwd; print c,nb,sum,dens,fold}'"'"' >>repeats/human/SST1/'$type'.sf.enrichment.tsv
      done 
    done 
  ' | sbatch -J $type --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done 

for nb in "DR" "STR" "MR" "TRI" "GQ" "Z"
do
  awk '($1=="chr13" && $2>12301367 && $3<12440010){print}' nonB_annotation/Homo_sapiens/autosomes_$nb.bed >repeats/human/SST1/chr13.$nb.bed
done 


################################################################################ 
# Enrichment in CenSat (possibly overlapping with the other repeat annotation)

# Use the "clean" bed files I extracted for the centromere analyis.
module load python/3.11.2
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
  cp centromeres/${sp}_pri/cenSat_shortnames.bed repeats/$sp/cenSat.bed
  cat repeats/$sp/cenSat.bed |python3 T2T_primate_nonB/python/repeat_summary.py >repeats/$sp/genome_cenSat_lengths.txt
done

# Get overlap 
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  for non_b in "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z" "all"
  do
      intersectBed -a repeats/'$sp'/cenSat.bed -b nonB_annotation/'${sp}'_pri/genome_${non_b}.bed >repeats/'${sp}'/overlap/genome_${non_b}_cenSat.bed
  done
  '| sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=1:00:00 --partition=open
done

# Calculate enrichment for the full genome
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
  echo '#!/bin/bash
  echo "Repeat APR DR GQ IR MR TRI STR Z" |sed "s/ /\t/g" >repeats/'$sp'_genome_censat_enrichment.tsv
  cat repeats/'$sp'/genome_cenSat_lengths.txt | while read -r rep replen;
  do
    tmp=`echo $rep`
    #Go through non-B
    cat densities/'$sp'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
    do
      d=`cat repeats/'$sp'/overlap/genome_${non_b}_cenSat.bed | awk -v r=$rep -v l=$replen -v dtot=$dens '"'"'($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print frac}}'"'"'`
      tmp=`echo $tmp" "$d`
      echo $tmp >tmp.'$sp'
    done
    cat tmp.'$sp' |sed "s/ /\t/g" >>repeats/'$sp'_genome_censat_enrichment.tsv
  done
  '| sbatch -J $sp --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
done

# Merge
echo "Species Repeat APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g' >repeats/7sp_genome_censat_enrichment.tsv
echo "Species Repeat Len" |sed 's/ /\t/g' >repeats/7sp_genome_censat_lengths.txt
cat T2T_primate_nonB/helpfiles/pri_species_list.txt |while read -r sp latin filename;
do
   awk -v sp=$sp '(NR>1){print sp,$0}' repeats/${sp}_genome_censat_enrichment.tsv |sed 's/^ //' |sed 's/ /\t/g' >>repeats/7sp_genome_censat_enrichment.tsv
   awk -v sp=$sp '{print sp,$0}' repeats/${sp}/genome_cenSat_lengths.txt |sed 's/^ //' |sed 's/ /\t/g' >>repeats/7sp_genome_censat_lengths.txt
done

# Plot censat enrichment with R/plot_figS6_enrichment_censat.R


################################################################################
# Check enrichment in repeats vs not repeats (all non-B combined) for main text

# ONLY REPEATMASKER FILES:
# Calculate non-B density inside and outside repeats sequence
#(only pri, no repeatmasker files for alternative haplotypes)
echo "Species NonB RepTotBp RepNonB RepDens NRTotBp NRNonB NRDens" |sed 's/ /\t/g' >repeats/rep_vs_nonrep_summary.txt
cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
do
  echo "looking at $sp"
  echo '#!/bin/bash
  module load bedtools/2.31.0
  #cut -f1,2 ref/'$filename'.fai |sort -k1,1 | complementBed -i <(sort -k1,1 -k2,2n repeats/'$sp'/RepeatMasker.bed) -g -  |grep -v "random" |grep -v "chrM"  >repeats/'$sp'/RepeatMasker_Complement.bed
  replen=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' repeats/'$sp'/RepeatMasker.bed`
  nrlen=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' repeats/'$sp'/RepeatMasker_Complement.bed`
  for nb in "all" "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z" 
  do
    rep=`intersectBed -a repeats/'$sp'/RepeatMasker.bed -b nonB_annotation/'$sp'_pri/genome_${nb}.bed  -wo | awk -v OFS="\t" -v rlen=$replen '"'"'{sum+=$8}END{dens=sum/rlen; print rlen,sum,dens}'"'"'`
    nr=`intersectBed -a repeats/'$sp'/RepeatMasker_Complement.bed -b nonB_annotation/'$sp'_pri/genome_${nb}.bed -wo | awk -v OFS="\t" -v nrlen=$nrlen '"'"'{sum+=$7}END{dens=sum/nrlen; print nrlen,sum,dens}'"'"'`
    echo '$sp'" "$nb" "$rep" "$nr |sed "s/ /\t/g" >>repeats/rep_vs_nonrep_summary.txt
  done
    ' |sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=16G --time=1:00:00 
done 

# COMBINED (this is what we used):
echo "Species NonB RepTotBp RepNonB RepDens NRTotBp NRNonB NRDens" |sed 's/ /\t/g' >repeats/allrep_vs_nonrep_summary.txt
cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
do
  echo "looking at $sp"
  echo '#!/bin/bash
  echo "looking at '$sp'"
  module load bedtools/2.31.0
  cat repeats/'$sp'/RepeatMasker.bed repeats/'$sp'/cenSat.bed repeats/'$sp'/composite_repeats.bed repeats/'$sp'/new_satellites.bed |cut -f1,2,3|sort -k1,1 -k2,2n |mergeBed -i -  >repeats/'$sp'/AllRepeats.bed
  cut -f1,2 ref/'$filename'.fai |sort -k1,1 | complementBed -i <(sort -k1,1 -k2,2n repeats/'$sp'/AllRepeats.bed) -g -  |grep -v "random" |grep -v "chrM"  >repeats/'$sp'/AllRepeats_Complement.bed
  replen=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' repeats/'$sp'/AllRepeats.bed`
  nrlen=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' repeats/'$sp'/AllRepeats_Complement.bed`
  for nb in "all" "APR" "DR" "GQ" "IR" "MR" "TRI" "STR" "Z" 
  do
    rep=`intersectBed -a repeats/'$sp'/AllRepeats.bed -b nonB_annotation/'$sp'_pri/genome_${nb}.bed  -wo | awk -v OFS="\t" -v rlen=$replen '"'"'{sum+=$7}END{dens=sum/rlen; print rlen,sum,dens}'"'"'`
    nr=`intersectBed -a repeats/'$sp'/AllRepeats_Complement.bed -b nonB_annotation/'$sp'_pri/genome_${nb}.bed -wo | awk -v OFS="\t" -v nrlen=$nrlen '"'"'{sum+=$7}END{dens=sum/nrlen; print nrlen,sum,dens}'"'"'`
    echo '$sp'" "$nb" "$rep" "$nr |sed "s/ /\t/g" >>repeats/allrep_vs_nonrep_summary.txt
  done
    ' |sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=16G --time=5:00:00
done 


# Calculate Fold enrichment and rearrange for supplementary table
echo "Species APR DR STR IR MR TRI G4 Z all" |sed 's/ /\t/g' >repeats/7sp_repeat_vs_nonrepeat_enrichment.txt
cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
do
  tmp="$sp"
  for nb in "APR" "DR" "STR" "IR" "MR" "TRI" "GQ" "Z" "all" 
  do
    f=`awk -v s=$sp -v n=$nb '($1==s && $2==n){fold=$5/$8; print fold}' repeats/allrep_vs_nonrep_summary.txt`
    tmp=$tmp" "$f
  done 
  echo $tmp >>repeats/7sp_repeat_vs_nonrepeat_enrichment.txt
done 

# Do Enrichment analysis for the different TYPES of repeats: 
# TE, Satellite, RNA, Other

# Prepare new bedfiles 
module load bedtools/2.31.0
cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
do
  mkdir -p repeats/$sp/GroupAnnotation
  awk -v OFS="\t" '{if($4~/^DNA/ || $4~/^LTR/ || $4~/^LINE/ ||$4~/^SINE/ ||$4~/^RC/ ||$4~/^Retroposon/){$4="TEs"}else if($4~/^Satellite/){$4="Satellites"}else if($4~/RNA/){$4="RNA"}else{$4="Other"}; print $0}' repeats/$sp/RepeatMasker.bed >repeats/$sp/RepeatMasker.groups.bed
  for type in "TEs" "RNA" "Other"
  do
    grep $type repeats/$sp/RepeatMasker.groups.bed |mergeBed -i - >repeats/$sp/GroupAnnotation/$type.bed
  done 
  grep "Satellites" repeats/$sp/RepeatMasker.groups.bed | cat - repeats/$sp/new_satellites.bed  <(grep -v "gap" repeats/$sp/cenSat.bed ) |sort -k1,1 -k2,2n |cut -f1-3 |mergeBed -i - >repeats/$sp/GroupAnnotation/Satellites.bed
done 

cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
do
  echo "Repeat APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g' >repeats/${sp}_group_enrichment.tsv
  echo '#!/bin/bash
  module load bedtools/2.31.0
  for rep in "TEs" "Satellites" "RNA" "Other"
  do
    replen=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' repeats/'$sp'/GroupAnnotation/$rep.bed`
    tmp=`echo $rep`
    cat densities/'$sp'_pri_nonB_genome_wide.txt |grep -v "all" | while read -r nb tot dens;
    do
      fold=`intersectBed -a repeats/'$sp'/GroupAnnotation/$rep.bed -b nonB_annotation/'$sp'_pri/genome_${nb}.bed  -wo | awk -v OFS="\t" -v rlen=$replen -v d=$dens '"'"'{sum+=$7}END{dens=sum/rlen; fold=dens/d; print fold}'"'"'`
      tmp=`echo $tmp" "$fold`
      echo $tmp >tmp.'$sp'
    done
    cat tmp.'$sp' |sed "s/ /\t/g" >>repeats/'$sp'_group_enrichment.tsv
  done 
  ' |sbatch -J $sp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=1:00:00
done 

# Merge before plotting 
echo "Species Repeat APR DR GQ IR MR TRI STR Z" |sed 's/ /\t/g' >repeats/7sp_group_enrichment.tsv
cat T2T_primate_nonB/helpfiles/pri_species_list.txt | while read -r sp latin filename;
do
  awk -v OFS="\t" -v sp=$sp '{print sp,$0}' <(tail -n+2 repeats/${sp}_group_enrichment.tsv) >>repeats/7sp_group_enrichment.tsv
done 

# Plot with R/plot_fig6ABC_repeats_methylation.R (part B)
