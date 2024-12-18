# CODE TO CALCULATE THE AMOUNT OF NEW SEQUENCE IN THE GENOME

mkdir -p circos/human/

# Create karyotype file
awk '{gsub(/chr/,""); print "chr - hs"$1,$1,"0",$2,"chr"$1}' ref/chm13v2.0.fa.fai |grep -v chrM >circos/human/karyotype.human.ucscCol.txt
# or in grey
awk '{gsub(/chr/,""); print "chr - hs"$1,$1,"0",$2,"vvlgrey"}' ref/chm13v2.0.fa.fai |grep -v hsM >circos/human/karyotype.human.txt
# With centromeres highlighted
awk '(/CEN/){gsub(/chr/,""); print "hs"$1,$2,$3}' ref/centromeres/chm13v2.0_GenomeFeature_v1.0.bed   >>circos/human/data/highlight_CEN.txt

# Create data files
for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
do
  cat densities/human_*_${nb}_100kb.bed |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "hs"); $3=$3-1; print $1,$2,$3,$4}' >circos/human/data/${nb}.100kb.txt
done

# Run circos
cd circos/human/
~/software/circos-0.69-9/bin/circos -conf circos_halfInnerCircle.conf

awk '(NR>1){if($4~/^hor/){if($4~/L)$/){print $1,$2,$3,"active_asat"}else{print $1,$2,$3,"active_asat"}}}else{split($4,s,"_"); print $1,$2,$3,s[1]} }' ref/repeats/chm13v2.0_censat_v2.0.bed |less



################################################################################
### BONOBO
mkdir -p bonobo_pri/data
grep -v random ref/mPanPan1.pri.cur.20231122.fasta.fai | awk '{gsub(/chr/,""); ; split($1,s,"_"); print "chr - pp"s[1],s[1],"0",$2,"vvlgrey"}'  >circos/bonobo_pri/karyotype.txt

# With centromeres highlighted
cut -f1  ref/mPanPan1.pri.cur.20231122.fasta.fai |grep -v random |grep -f - ref/centromeres/mPanPan1_v2.0.GenomeFeature_v1.0.bed | awk '(/Cen/){gsub(/chr/,""); split($1,s,"_"); print "pp"s[1],$2,$3}' >>circos/bonobo_pri/data/highlight_CEN.txt

for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
do
  cat densities/bonobo_pri/*_${nb}_100kb.bed |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "pp"); split($1,s,"_"); $3=$3-1; print s[1],$2,$3,$4}' >circos/bonobo_pri/data/${nb}.100kb.txt
done

################################################################################
### ALL OTHER SPECIES
#Naming chromosomes "nn", doesn't matter if all sp have the sama names
for hap in "alt" #"pri" #
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep "siamang" |grep -v "human" |while read -r sp latin filename;
  do
#    cp -r circos/template circos/${sp}_${hap}/
#    grep -v random ref/$filename.fai | awk '{gsub(/chr/,""); ; split($1,s,"_"); print "chr - nn"s[1],s[1],"0",$2,"vvlgrey"}'  >circos/${sp}_${hap}/karyotype.txt
    for nb in "APR" "DR" "GQ" "IR" "MR" "STR" "Z"
    do
      cat densities/${sp}_${hap}/*_${nb}_100kb.bed |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "nn"); split($1,s,"_"); $3=$3-1; print s[1],$2,$3,$4}' >circos/${sp}_${hap}/data/${nb}.100kb.txt
    done
#    spshort=`echo $filename |cut -f1 -d"."`
#    cut -f1  ref/$filename.fai |grep -v random |grep -f - ref/centromeres/$spshort*GenomeFeature*.bed | awk '(/Cen/){gsub(/chr/,""); split($1,s,"_"); print "nn"s[1],$2,$3}' >circos/${sp}_${hap}/data/highlight_CEN.txt
  done
done

# Run Circos for each of them
for hap in "alt" #"pri" #
do
  cat T2T_primate_nonB/helpfiles/${hap}_species_list.txt |grep -v "human" |while read -r sp latin filename;
  do
    echo " Looking at $sp, $hap"
    cd circos/${sp}_${hap}/
    ~/software/circos-0.69-9/bin/circos -conf circos_halfInnerCircle.conf
    mv circos.png ${sp}_${hap}_circos.png
    cd ../..
  done
done
