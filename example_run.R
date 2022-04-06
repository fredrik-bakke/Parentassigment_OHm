####### Testing the example genotype file ########
inpgeno <- "example/example_001"
parentfile <- "example/parents.csv"
qc <- c(geno = 0.05, mind = 0, maf = 0.01, hwe = 1e-6, thin = 0.999, chrset = 30)
threshOMM <- 25
matchchecks <- NULL
outfile <- "example_001_OH"
outfolder <- "out"

source("OHmassign.R")
example_001_OH <- OHm(inpgeno, parentfile, qc, threshOMM, matchchecks, outfile, outfolder)

### import the orginal pedigree
origPED <- read.table(paste(inpgeno, ".fam", sep = ""), header = F, stringsAsFactors = F)
colnames(origPED) <- c("ID", "sireOrig", "damOrig")

#### merge the results from the OHm with orignal pedigree
checkPED <- merge(example_001_OH, origPED, by = "ID")

#####################################################################################
#### check mismatches for sires and dams (omit unassigned sires or dams)
siremismatch <- na.omit(checkPED[, c("ID", "sire", "sireOrig")])
siremismatch$check <- ifelse(siremismatch$sire == siremismatch$sireOrig, 1, 0)

## % mismatch for sires
cat("Mismath for sires:", 1 - sum(siremismatch$check) / nrow(siremismatch), "\n")

## % mismatch for sires over all data
cat("Mismath for sires over all data:", 1 - sum(siremismatch$check) / nrow(checkPED), "\n")
######################################################################################

######################################################################################
#### check mismatches for sires and dams (omit unassigned sires or dams)
dammismatch <- na.omit(checkPED[, c("ID", "dam", "damOrig")])
dammismatch$check <- ifelse(dammismatch$dam == dammismatch$damOrig, 1, 0)

## % mismatch for dams
cat("Mismath for dams:", 1 - sum(dammismatch$check) / nrow(dammismatch), "\n")

## % mismatch for dams over all data
cat("Mismath for dams over all data:", 1 - sum(dammismatch$check) / nrow(checkPED), "\n")
#####################################################################################

################# Specific checks on the mismatches  ###################
sirewrong <- siremismatch[which(siremismatch$check == 0), ]
sirewrong <- merge(sirewrong, example_001_OH, by.x = 1, by.y = 1)

damwrong <- dammismatch[which(dammismatch$check == 0), ]
damwrong <- merge(damwrong, example_001_OH, by.x = 1, by.y = 1)