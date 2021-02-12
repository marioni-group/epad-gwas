# Generate list of related pairs and filter based on csf availability
setwd("/Cluster_Filespace/Marioni_Group/EPAD_GWAS/")
rel = read.table("PLINK_COMBINED/relatedness.genome", header=T) 

# Remove if already a QC failure/Ancestry failure
ancestry = read.table("PLINK_COMBINED/ancestry_outliers.txt", header=T)

# List of valid samples
# NA any invalid samples
valid = read.table("PLINK_COMBINED/EPAD.QC.valid", header=F)
valid = valid[-which(valid$V2 %in% ancestry$IID), ]

rel$IID1[which(!rel$IID1 %in% valid$V2)] = NA
rel$IID2[which(!rel$IID2 %in% valid$V2)] = NA

# threshold 0.1875 represents the half-way point between 2nd and 3rd degree relatives and is a common cut-off to use
rel = rel1 = rel[which(rel$PI_HAT > 0.1875), ]


csf = read.csv("Phenotypes/Phenotypes_Bert_02Nov2020/v_imi_epadlcs_csf.csv")

array = read.delim("Samples_16Dec2020.txt", header=T)
missing = csf[which(csf$abeta_1_42_result=="" & !is.na(csf$ID)),"ID"]

csf$ID = array[match(csf$csf_sample_id, array$Name2), "Sample.ID"]

csf = csf[which(csf$ptau_result!=""|csf$ttau_result!=""|csf$abeta_1_42_result!=""),]

rel$csf1 = rel$csf2 = NA
rel$csf1[which(rel$IID1 %in% csf$ID)] = "Yes"
rel$csf2[which(rel$IID2 %in% csf$ID)] = "Yes"

# Check for those with multiple relationships if CSF data is available for both IDs in each pair
counts = data.frame(table(c(rel$IID1, rel$IID2)))
counts = counts[order(counts$Freq, decreasing=T),]
multi = as.character(counts[which(counts$Freq>1), "Var1"])

rel[which(rel$IID1 %in% multi | rel$IID2 %in% multi),]

# No CSF pairs. Only 204390430081_R09C01  has CSF data

# Remove ID if no csf data
# Prioritise keeping males to shift the sex imbalance
ids = data.frame(table(c(rel$IID1, rel$IID2)))
ids$Var1 = as.character(ids$Var1)
ids$csf = "N"
ids$csf[ids$Var1 %in% csf$ID] = "Y"
fam = read.table("PLINK_COMBINED/EPAD.QC.fam")
ids$sex = fam[match(ids$Var1, fam$V2), "V5"]

# Aim: each row from rel should have at least one IID removed
rm_csf = ids$Var1[which(ids$csf=="N")] # remove those with no CSF data
rel$IID1[rel$IID1 %in% rm_csf] = NA
rel$IID2[rel$IID2 %in% rm_csf] = NA

# Remove females in male-female pairings (overall sample has more females than males)
rel$sex1 = fam[match(rel$IID1, fam$V2),"V5"]
rel$sex2 = fam[match(rel$IID2, fam$V2),"V5"]

rm_female = c(rel[which(rel$sex1==1 & rel$sex2==2), "IID2"], rel[which(rel$sex2==1 & rel$sex1==2), "IID1"])

rel$IID1[rel$IID1 %in% rm_female] = NA
rel$IID2[rel$IID2 %in% rm_female] = NA

# Take best call rate for remaining pairs (should be same-sex pairs both with CSF data)
rel$call1 = array[match(rel$IID1, array$Sample.ID), "Call.Rate"]
rel$call2 = array[match(rel$IID2, array$Sample.ID), "Call.Rate"]

tmp_rel = rel[which(!is.na(rel$IID1) & !is.na(rel$IID2)),]
mins = as.numeric(apply(tmp_rel[,c("call1","call2")],1 ,which.min))
# remove those with NA entries in call rate (i.e. already filtered rows)

rm_mins = c(tmp_rel$IID1[which(mins==1)], tmp_rel$IID2[which(mins==2)])

rel$IID1[rel$IID1 %in% rm_mins] = NA
rel$IID2[rel$IID2 %in% rm_mins] = NA
keep = na.omit(c(rel$IID1, rel$IID2))
keep = unique(keep)
remove = c(rel1$IID1, rel1$IID2)
remove = na.omit(unique(remove[-which(remove %in% keep)]))


# Make EPAD.GWAS.valid file
valid = valid[-which(valid$V2 %in% remove | valid$V2 %in% ancestry$IID), ]
write.table(valid, file="PLINK_COMBINED/EPAD.GWAS.valid", sep=' ', quote=F, row.names=F) # 1845 samples
####################################################################################################################

####################################################################################################################
# Make plink GWAS files with ancestry/relatedness filters
plink19 --make-bed --bfile EPAD.QC --keep EPAD.GWAS.valid --out EPAD.GWAS


## Make linker file (in R)
setwd("/Cluster_Filespace/Marioni_Group/EPAD_GWAS/PLINK_COMBINED")
info  = read.delim("../Samples_16Dec2020.txt")

info = info[which(!is.na(info$Name2)),]
link = cbind(info[,c("Name2", "Sample.ID")])
colnames(link) = c("CSF_ID", "GWAS_ID")

write.table(link, file="linkfile.txt", row.names=F, quote=F)