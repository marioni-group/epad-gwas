# Download Hapmap reference data and merge with EPAD QC'd data for ancestry PCA ###
# Protocol taken from https://cran.r-project.org/web/packages/plinkQC/vignettes/HapMap.pdf

##########################################
##### Download reference data (Hapmap) ###
##########################################
cd /Cluster_Filespace/Marioni_Group/EPAD_GWAS/refdir

ftp=ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/
prefix=hapmap3_r2_b36_fwd.consensus.qc.poly

wget $ftp/$prefix.map.bz2
bunzip2 $prefix.map.bz2

wget $ftp/$prefix.ped.bz2
bunzip2 $prefix.ped.bz2

wget $ftp/relationships_w_pops_121708.txt

# Make bed files
plink19 --file ./$prefix \
      --make-bed \
      --out ./HapMapIII_NCBI36
mv ./HapMapIII_NCBI36.log ./log

# Update annotation to hg19 (download hg18 to hg19 liftover file from ucsc)
# http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/
gunzip hg18ToHg19.over.chain.gz

# Convert chrX/Y to 23/24
awk '{print "chr" $1, $4 -1, $4, $2 }' ./HapMapIII_NCBI36.bim | \
    sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > \
    ./HapMapIII_NCBI36.tolift

# Liftover
liftOver ./HapMapIII_NCBI36.tolift ./hg18ToHg19.over.chain \
    ./HapMapIII_CGRCh37 ./HapMapIII_NCBI36.unMapped    

# extract mapped variants
awk '{print $4}' ./HapMapIII_CGRCh37 > ./HapMapIII_CGRCh37.snps
# ectract updated positions
awk '{print $4, $3}' ./HapMapIII_CGRCh37 > ./HapMapIII_CGRCh37.pos

# Update the reference data
plink19 --bfile ./HapMapIII_NCBI36 \
    --extract ./HapMapIII_CGRCh37.snps \
    --update-map ./HapMapIII_CGRCh37.pos \
    --make-bed \
    --out ./HapMapIII_CGRCh37

	mv ./HapMapIII_CGRCh37.log ./log




##########################################
##########################################
##### Merge data and run PCA   ###########
##########################################

cd /Cluster_Filespace/Marioni_Group/EPAD_GWAS
mkdir qcdir
mkdir qcdir/plink_log
cp ./PLINK_COMBINED/EPAD.QC.bim ./qcdir/
cp ./PLINK_COMBINED/EPAD.QC.bed ./qcdir/
cp ./PLINK_COMBINED/EPAD.QC.fam ./qcdir/


# Merge datasets
qcdir='./qcdir'
refdir='./refdir'
name='EPAD.QC'
refname='HapMapIII_CGRCh37'



# Filter reference and study data for non A-T or G-C SNPs
##########################################################
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $qcdir/$name.bim  > \
    $qcdir/$name.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $refdir/$refname.bim  > \
    $qcdir/$refname.ac_gt_snps
   
plink19 --bfile  $refdir/$refname \
      --exclude $qcdir/$refname.ac_gt_snps \
      --make-bed \
      --out $qcdir/$refname.no_ac_gt_snps

mv  $qcdir/$refname.no_ac_gt_snps.log $qcdir/plink_log/$refname.no_ac_gt_snps.log

plink19 --bfile  $qcdir/$name \
      --exclude $qcdir/$name.ac_gt_snps \
      --make-bed \
      --out $qcdir/$name.no_ac_gt_snps
mv  $qcdir/$name.no_ac_gt_snps.log $qcdir/plink_log/$name.no_ac_gt_snps.log

#####################################################################

#### Prune Data ###
plink19 --bfile  $qcdir/$name.no_ac_gt_snps \
      --exclude range /home/dmccartn/R/x86_64-suse-linux-gnu-library/3.6/plinkQC/extdata/high-LD-regions-hg19-GRCh37.txt \
      --indep-pairwise 50 5 0.2 \
      --out $qcdir/$name.no_ac_gt_snps

mv  $qcdir/$name.no_ac_gt_snps.log $qcdir/plink_log/$name.no_ac_gt_snps.pruned.log 

plink19 --bfile  $qcdir/$name.no_ac_gt_snps \
      --extract $qcdir/$name.no_ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/$name.pruned
mv  $qcdir/$name.pruned.log $qcdir/plink_log/$name.pruned.log

#########################################################

### Filter reference data for the same SNP set as in study ###
plink19 --bfile  $refdir/$refname \
      --extract $qcdir/$name.no_ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/$refname.pruned
mv  $qcdir/$refname.pruned.log $qcdir/plink_log/$refname.pruned.log

#########################################################

### Check and correct chromosome mismatch ###

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim | \
    sed -n '/^[XY]/!p' > $qcdir/$refname.toUpdateChr

plink19 --bfile $qcdir/$refname.pruned \
      --update-chr $qcdir/$refname.toUpdateChr 1 2 \
      --make-bed \
      --out $qcdir/$refname.updateChr
mv $qcdir/$refname.updateChr.log $qcdir/plink_log/$refname.updateChr.log

#################################################################

### Position mismatch ### 
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/${refname}.toUpdatePos

#################################################################

### Possible allele flips ### 
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/$refname.toFlip

#################################################################

### Update positions and flip alleles ### 

plink19 --bfile $qcdir/$refname.updateChr \
      --update-map $qcdir/$refname.toUpdatePos 1 2 \
      --flip $qcdir/$refname.toFlip \
      --make-bed \
      --out $qcdir/$refname.flipped
mv $qcdir/$refname.flipped.log $qcdir/plink_log/$refname.flipped.log

#################################################################

### Remove mismatches ### 

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.flipped.bim > \
    $qcdir/$refname.mismatch

plink19 --bfile $qcdir/$refname.flipped \
      --exclude $qcdir/$refname.mismatch \
      --make-bed \
      --out $qcdir/$refname.clean
mv $qcdir/$refname.clean.log $qcdir/plink_log/$refname.clean.log

#################################################################

### Merge study genotype and reference data ###
plink19 --bfile $qcdir/$name.pruned  \
      --bmerge $qcdir/$refname.clean.bed $qcdir/$refname.clean.bim \
         $qcdir/$refname.clean.fam  \
      --make-bed \
      --out $qcdir/$name.merge.$refname
mv $qcdir/$name.merge.$refname.log $qcdir/plink_log


#################################################################

### PCA on the merged data ###
plink19 --bfile $qcdir/$name.merge.$refname \
      --pca \
      --out $qcdir/$name.$refname
mv $qcdir/$name.$refname.log $qcdir/plink_log

#################################################################


##########################################
##### Ancestry check in R ################
##########################################


# Plink QC in R (Ancestry check)
require(plinkQC)
setwd("/Cluster_Filespace/Marioni_Group/EPAD_GWAS/")
package.dir <- find.package('plinkQC')
indir = "qcdir"
qcdir <- tempdir()
name <- 'EPAD.QC'
path2plink <- "/usr/local/bin/plink19"

# Ancestry check
indir <- "qcdir"
refname <- 'HapMapIII_CGRCh37'
prefixMergedDataset <- paste(name, ".", refname, sep="")
extdir = "/home/dmccartn/R/x86_64-pc-linux-gnu-library/4.0/plinkQC/extdata"

exclude_ancestry <-
    evaluate_check_ancestry(indir=indir, name=name, europeanTh=3, # EuropeanTh default=1.5 but eyeballing plot suggests this is too conservative
                            prefixMergedDataset=prefixMergedDataset,
                            refSamplesFile=paste(extdir, "/HapMap_ID2Pop.txt",
                                                 sep=""),
                            refColorsFile=paste(extdir, "/HapMap_PopColors.txt",
                                                sep=""),
                            interactive=TRUE)
 anc_fails = exclude_ancestry$fail_ancestry
 write.table(anc_fails, file="PLINK_COMBINED/ancestry_outliers.txt", sep=" ", quote=F, row.names=F)


