####### Testing the example genotype file ########
inpgeno='example_001';  
parentfile='parents.txt';
qc=c(geno=0.05,mind=0,maf=0.01,hwe=1e-6,thin=0.999,chrset=30);  
threshOMM=25;  
matchchecks=F;  
outfile='example_001_OH';  

source('OHmassign.R')  

example_001_OH <- OHm(inpgeno,parentfile,qc,threshOMM,matchchecks,outfile)  

### import the orginal pedigree  
origPED <- read.table('example_001.pedigree',header=F,stringsAsFactors = F)  
colnames(origPED) <- c('ID','sireOrig','damOrig')  

#### merge the results from the OHm with orignal pedigree  
checkPED <- merge(example_001_OH,origPED,by='ID')  

#####################################################################################
#### check mismatches for sires and dams (omit unassigned sires or dams)
siremismatch <- na.omit(checkPED[,c('ID','sire','sireOrig')])
siremismatch$check <- ifelse(siremismatch$sire==siremismatch$sireOrig,1,0)

## % mismatch for sires 
1-sum(siremismatch$check)/nrow(siremismatch)

## % mismatch for sires overall data
1-sum(siremismatch$check)/nrow(checkPED)
######################################################################################

######################################################################################
#### check mismatches for sires and dams (omit unassigned sires or dams)
dammismatch <- na.omit(checkPED[,c('ID','dam','damOrig')])
dammismatch$check <- ifelse(dammismatch$dam==dammismatch$damOrig,1,0)

## % mismatch for sires 
1-sum(dammismatch$check)/nrow(dammismatch)

## % mismatch for sires overall data
1-sum(dammismatch$check)/nrow(checkPED)
#####################################################################################

################# Specific checks on the mismatches  ###################
sirewrong <- siremismatch[which(siremismatch$check==0),]
sirewrong <- merge(sirewrong,example_001_OH,by.x=1,by.y=1)

damwrong <- dammismatch[which(dammismatch$check==0),]
damwrong <- merge(damwrong,example_001_OH,by.x=1,by.y=1)
