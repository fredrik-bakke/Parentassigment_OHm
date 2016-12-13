## Parentage assignment with opposite homozygote (OH) method
#### Part of the script is sourced from [Ferdosi M. and Boerner V. 2014](http://www.sciencedirect.com/science/article/pii/S1871141314002625) to undertake fast OH pairwise counts  
#### Plink v2 is also used to undertake quality control on the genotype data and also transform the data from binary/PED allele (AA/AB/BA/BB - 11/12/21/22) to genotype format (0/1/2).


### Requirement 
- Plink v2  
- Plink binary file format  
- comma separated (.csv) text file containing parental IDs and sex(M/F or 1/2). There should be a header  

### Argument  
- inpgeno='MHrecode'
- parentfile='parents.txt'
- qc=c(geno=0.05,mind=0.20,maf=0.05,hwe=1e-6,thin=0.9999,chrset=30);
- threshOMM=100
- matchchecks=F
- outfile='MHparassign'

##### inpgeno :: The prefix of the PLINK binary file  
##### parentfile :: comma separated file with header (ID, Sex) and genotyped parents.  See the example below
      eg ID,Sex  
        1,M  
        2,M  
        3,F  

##### qc :: Quality check to be undertaken on the genotype data. All qc are undertaken with PLINK, thus further information about the parameters could be read from there
- geno : SNP genotype call rate. please use a high value when no strong qc is required
- mind : Sample call rate. please use a high value when no strong qc is required
- maf  : Minor allele frequency 
- hwe  : Hardy Weinberg Equilibrium Fishers Exact Pvalue threshold
- thin : proportion of SNPs to keep for OH analyis, if marker data is too large thin can help reduce it
- chrset : The chromosome number


