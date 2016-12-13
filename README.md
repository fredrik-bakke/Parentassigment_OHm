# Parentage assignment with opposite homozygote (OH) method
### Part of the script is sourced from Ferdosi M. and Boerner V. 2014 to undertake fast OH pairwise counts  
### Plink v2 is also used to undertake quality control on the genotype data and also transform the data from binary/PED allele (AA/AB/BA/BB - 11/12/21/22) to genotype format (0/1/2).


### Requirement 
- Plink v2
- Plink binary file format
- comma separated (.csv) text file containing parental IDs and sex(M/F or 1/2). There should be a header
    eg ID,Sex
        1,M
        2,M
        3,F


