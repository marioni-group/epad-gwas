# csf GWAS phenotype - residualised for age, sex and 20PCs, study centre?
setwd("/Cluster_Filespace/Marioni_Group/EPAD_GWAS/qcdir")
pcs = read.table("EPAD.QC.eigenvec")

csf = read.csv("../Phenotypes/Phenotypes_Bert_02Nov2020/v_imi_epadlcs_csf.csv")
csf = csf[order(csf$patient_id, csf$visit),]

# Take first visit with a measurement for each id
csf_keep = list()
for(i in unique(csf$patient_id)){
	tmp = csf[which(csf$patient_id==i),]
	v1 = which(tmp$ptau_result!="")[1]
	csf_keep[[i]] = tmp[v1,]
}

baseline = do.call("rbind", csf_keep)

pt = read.csv("../Phenotypes/Phenotypes_Bert_02Nov2020/v_imi_epadlcs_sample_ids_and_socio_demographics.csv")
map = read.delim("../Samples_16Dec2020.txt")

pt$ID = map[match(pt$apoe_blood_sample_id, map$Name2), "Sample.ID"]
pt$age2 = pt$age_years + (pt$age_months/12)
pt = merge(pt, pcs[,2:22], by.x="ID", by.y="V2")

pt = merge(pt, csf, by.x="apoe_blood_sample_id", by.y="csf_sample_id")



# ptau
# Cap min and max values
 pt[grep(">|<", pt$ptau_result),"ptau_result"]
#  [1] "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8"
# [16] "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8"
# [31] "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8"
# [46] "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8" "<8"
 pt[grep(">|<", pt$ptau_result),"ptau_result"] = "8"
 pt$ptau_result = as.numeric(pt$ptau_result)

 # Trim outliers (5SD)
 pt$ptau_result[which(abs(scale(pt$ptau_result))>=5)] = NA

# ttau
# Cap min and max values
 pt[grep(">|<", pt$ttau_result),"ttau_result"]
# [1] "<80" "<80" "<80" "<80" "<80" "<80" "<80" "<80" "<80" "<80" "<80"

 pt[grep(">|<", pt$ttau_result),"ttau_result"] = "80"
 pt$ttau_result = as.numeric(pt$ttau_result)

 # Trim outliers (5SD)
 pt$ttau_result[which(abs(scale(pt$ttau_result))>=5)] = NA

 pt$ptau_ttau = pt$ptau_result/pt$ttau_result


 # Abeta capped at 200 and 1700
 # Comments present for those >1700 with 2nd measurement values - import these to phenotype
 pt[grep("<", pt$abeta_1_42_result),"abeta_1_42_result"] = "200"
  pt[grep(">", pt$abeta_1_42_result),"abeta_1_42_result"] = pt[grep(">", pt$abeta_1_42_result),"abeta_1_42_comments"]
pt$abeta_1_42_result = gsub(".*= | *PG.ML", "", pt$abeta_1_42_result)
 pt$abeta_1_42_result = as.numeric(pt$abeta_1_42_result)

 # Trim outliers (5SD)
 pt$abeta_1_42_result[which(abs(scale(pt$abeta_1_42_result))>=5)] = NA


# Remove QC fails
keep = read.table("../PLINK_COMBINED/EPAD.GWAS.valid")
pt = pt[which(pt$ID %in% keep[,2]), ]
rownames(pt) = pt$ID

fam = read.table("../PLINK_COMBINED/EPAD.QC.fam")
fam = fam[which(fam$V2 %in% pt$ID),]

test = read.table("/Cluster_Filespace/Marioni_Group/EPAD_GWAS/PLINK_COMBINED/imputation/HRC/test.fam")
test$V1 = sub("*._|*.._|*..._|*...._|", "", test$V1)
  # gwas phenotypes
abeta_gwas = resid(lm(abeta_1_42_result~age2 + as.factor(sex) + as.factor(site_id) + ptau_result + ttau_result +  V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + 
	                           V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22, data=pt))


ptau_gwas = resid(lm(ptau_result~age2 + as.factor(sex) + as.factor(site_id) + abeta_1_42_result + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + 
	                           V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22, data=pt))

ttau_gwas = resid(lm(ttau_result~age2 + as.factor(sex) + as.factor(site_id) +  abeta_1_42_result + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + 
	                           V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22, data=pt))

ptau_ttau_gwas = resid(lm(ptau_ttau~age2 + as.factor(sex) + as.factor(site_id) +  abeta_1_42_result + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + 
	                           V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22, data=pt))




abeta_phen = data.frame(fam[match(names(abeta_gwas), fam$V2),"V1"], names(abeta_gwas), abeta_gwas)
ptau_phen = data.frame(fam[match(names(ptau_gwas), fam$V2),"V1"], names(ptau_gwas), ptau_gwas)
ttau_phen = data.frame(fam[match(names(ttau_gwas), fam$V2),"V1"], names(ttau_gwas), ttau_gwas)
ptau_ttau_phen = data.frame(fam[match(names(ptau_ttau_gwas), fam$V2),"V1"], names(ptau_ttau_gwas), ptau_ttau_gwas)

# dir.create("../GWAS_Runs")
	write.table(abeta_phen, file="../GWAS_Runs/abeta_gwas.phen", sep=" ", col.names=F, row.names=F, quote=F)
	write.table(ttau_phen, file="../GWAS_Runs/ttau_gwas.phen", sep=" ", col.names=F, row.names=F, quote=F)
	write.table(ptau_phen, file="../GWAS_Runs/ptau_gwas.phen", sep=" ", col.names=F, row.names=F, quote=F)
	write.table(ptau_ttau_phen, file="../GWAS_Runs/ptau_ttau_gwas.phen", sep=" ", col.names=F, row.names=F, quote=F)

# Generate ped files for HRC GWAS
map$sex = map$Gender2
map$sex = gsub("Female", 2, map$sex)
map$sex = gsub("Male", 1, map$sex)
abeta_ped = data.frame(test[match(abeta_phen[,2], test$V1), "V2"], test[match(abeta_phen[,2], test$V1), "V2"],  abeta_phen[,3])
ttau_ped = data.frame(test[match(ttau_phen[,2], test$V1), "V2"], test[match(ttau_phen[,2], test$V1), "V2"],  ttau_phen[,3])
ptau_ped = data.frame(test[match(ptau_phen[,2], test$V1), "V2"], test[match(ptau_phen[,2], test$V1), "V2"], ptau_phen[,3])
ptau_ttau_ped = data.frame(test[match(ptau_ttau_phen[,2], test$V1), "V2"], test[match(ptau_ttau_phen[,2], test$V1), "V2"], ptau_ttau_phen[,3])

write.table(abeta_ped, file="/Cluster_Filespace/Marioni_Group/EPAD_GWAS/GWAS_Runs/HRC/ABeta_EPAD.ped", col.names=F, row.names=F, quote=F)
write.table(ttau_ped, file="/Cluster_Filespace/Marioni_Group/EPAD_GWAS/GWAS_Runs/HRC/TTau_EPAD.ped", col.names=F, row.names=F, quote=F)
write.table(ptau_ped, file="/Cluster_Filespace/Marioni_Group/EPAD_GWAS/GWAS_Runs/HRC/PTau_EPAD.ped", col.names=F, row.names=F, quote=F)
write.table(ptau_ttau_ped, file="/Cluster_Filespace/Marioni_Group/EPAD_GWAS/GWAS_Runs/HRC/Ratio_EPAD.ped", col.names=F, row.names=F, quote=F)
