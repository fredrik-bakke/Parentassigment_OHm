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
sirematch <- na.omit(checkPED[, c("ID", "sire", "sireOrig")])
sirematch$check <- ifelse(sirematch$sire == sirematch$sireOrig, 1, 0)

## % mismatch for sires
cat("Proportion of matches for sires:", 1 - sum(sirematch$check) / nrow(sirematch), "\n")

## % mismatch for sires over all data
cat("Proportion of matches for sires over all data:", 1 - sum(sirematch$check) / nrow(checkPED), "\n")
######################################################################################

######################################################################################
#### check mismatches for sires and dams (omit unassigned sires or dams)
dammatch <- na.omit(checkPED[, c("ID", "dam", "damOrig")])
dammatch$check <- ifelse(dammatch$dam == dammatch$damOrig, 1, 0)

## % mismatch for dams
cat("Proportion of matches for dams:", 1 - sum(dammatch$check) / nrow(dammatch), "\n")

## % mismatch for dams over all data
cat("Proportion of matches for dams over all data:", 1 - sum(dammatch$check) / nrow(checkPED), "\n")
#####################################################################################

################# Specific checks on the mismatches  ###################
sirewrong <- sirematch[which(sirematch$check == 0), ]
sirewrong <- merge(sirewrong, example_001_OH, by.x = 1, by.y = 1)

damwrong <- dammatch[which(dammatch$check == 0), ]
damwrong <- merge(damwrong, example_001_OH, by.x = 1, by.y = 1)