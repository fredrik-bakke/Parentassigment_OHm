## Parentage assignment with opposite homozygote (OH) method

#### Part of the script is sourced from [Ferdosi M. and Boerner V. 2014](http://www.sciencedirect.com/science/article/pii/S1871141314002625) to undertake fast OH pairwise counts  

#### Plink v2 is also used to undertake quality control on the genotype data and also transform the data from binary/PED allele (AA/AB/BA/BB - 11/12/21/22) to genotype format (0/1/2)

### Requirement

     - Plink v2  (PLINK windows version is need to be in the folder where the analysis is been undertaken)
     - Plink binary file format  for the genotype data
     - comma separated (.csv) text file containing parental IDs and sex(M/F or 1/2). There should be a header (ID, Sex)  

### Argument  (see the example_run.R file on how to run the example file)

      - inpgeno='example_001'  
      - parentfile='parents.csv'  
      - qc=c(geno=0.05,mind=0.20,maf=0.05,hwe=1e-6,thin=0.9999,chrset=30)  
      - threshOMM=100  
      - matchchecks=F  
      - outfile='example_001_OH'  

#### inpgeno :: The prefix of the PLINK binary file  

#### parentfile :: comma separated file with header (ID, Sex) and genotyped parents.  See the example below  

      example format of the parentfile  
       ID,Sex  
        1,M  
        2,M  
        3,F  

#### qc :: Quality check to be undertaken on the genotype data. All qc are undertaken with PLINK, thus further information about the parameters could be read from there  

       - geno : SNP genotype call rate. please use a high value when no strong qc is required  
       - mind : Sample call rate. please use a high value when no strong qc is required  
       - maf  : Minor allele frequency   
       - hwe  : Hardy Weinberg Equilibrium Fishers Exact Pvalue threshold  
       - thin : proportion of SNPs to keep for OH analyis, if marker data is too large thin can help reduce it  
       - chrset : The chromosome number  

#### threshOMM :: The number of SNPs that is expected to be genotyping error. It also means that, if parent-offspring comparison exceeds this number, they are not assigned.   It is important to look at the graph from the first run to inform you on the number of markers to allow as genotyping error. An example graph (graph.png) shows that, there is high chance of zero genotyping error, however, some OH combinations can be aroound 50, thus threshOMM can be set to 50 SNPs  

#### matchchecks :: This argument is used to signal for computing opposite homozygotes in the data of to check for sepcific parent offspring relationships. The parameter should be F (FALSE), when you want to compute OH and assign pedigree. Alternatively it should be a file containing the specific pedigree you want to check. That file should be a comma separated file with a header  

      matchchecks=F  
      matchchecks='pedigcheck.csv'   
      
      example format of the matchchecks file  
          ID,checkID
          10,1  
          11,5   
          12,5 
          13,2  

#### outfile :: output file name

### Output files  

- Pedigree file
- png plot of number of opposite homozygotes  
- Parent file after qc (always called 'parentsafterqc.csv')  

##### The generated pedigree with the following headers  

         - ID                    : The IDs of the offspring   
         - sire                  : The sire (from the opposite homozygote [OH] approach)  
         - OHsire                : Number of markers that are OH between animal (ID) and sire (pseudo-genotyping error)  
         - dam                   : The dam (from the opposite homozygote [OH] approach)  
         - OHdam                 : Number of markers that are OH between animal (ID) and dam (pseudo-genotyping error)  
         - sirepossib            : Alternative sires, this happens when another sire have low number of OH and below the threshold  
         - OHsirepossib          : Number of markers that are OH between animal (ID) and alternate sires (pseudo-genotyping error)  
         - dampossib             : Alternative dams, this happens when another dam have low number of OH and below the threshold  
         - OHdampossib           : Number of markers that are OH between animal (ID) and alternate dams (pseudo-genotyping error)  

### Disclaimer : You use the script at your own risk :)

#### The script was written in close discussions with

      Luqman Aslam (luqman.aslam@nofima.no)
      Matthew Baranski (matt.baranski@marineharvest.com)
