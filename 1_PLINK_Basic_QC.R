# Generate list of QC'd samples and SNPs:
# https://choishingwan.github.io/PRS-Tutorial/target/
# Caroline H suggested MAF could be decreased to 0.005 (now more commonly used in GWAS) from 0.01

plink19 \
    --bfile EPAD \
    --maf 0.005  \
    --hwe 1e-6 \
    --geno 0.01 \
    --chr 1-22, X, XY, Y \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EPAD.QC


# Prune data:
plink19 --bfile EPAD --keep EPAD.QC.fam --extract EPAD.QC.snplist --indep-pairwise 200 50 0.25 --out EPAD.QC  

# Compute heterozygosity rates
plink19 \
    --bfile EPAD \
    --extract EPAD.QC.prune.in \
    --keep EPAD.QC.fam \
    --het \
    --out EPAD.QC   


# Remove samples with high heterozygosity (in R)
dat <- read.table("EPAD.QC.het", header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], "EPAD.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
q() # exit R
###############################################################################################################

# # Check sex
plink19 --bfile EPAD --keep EPAD.valid.sample --check-sex --out EPAD.QC
    

# ### YCHR SEXCHECK FOR COMPARISON ### Not likely to use this
# plink19 --check-sex y-only --bfile EPAD --keep EPAD.valid.sample  --out checksex_y_test


# Read in file (in R)
valid <- read.table("EPAD.valid.sample", header=T)
dat <- read.table("EPAD.QC.sexcheck", header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], "EPAD.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) # Overwrite EPAD.QC.valid with sexcheck passes
q() # exit R
###############################################################################################################


# Generate map and ped files for --genome flag
plink19 --recode --bfile EPAD --out EPAD --extract EPAD.QC.prune.in
# Generate table of pairwise relationships
plink19 --genome --file EPAD --keep EPAD.QC.valid --out relatedness
	
# threshold 0.1875 represents the half-way point between 2nd and 3rd degree relatives and is a common cut-off to use

# Make QC-filtered bim/bed/bam files using SNPlist and valid samples file
plink19 --make-bed --bfile EPAD --extract EPAD.QC.snplist --keep EPAD.QC.valid --out EPAD.QC
