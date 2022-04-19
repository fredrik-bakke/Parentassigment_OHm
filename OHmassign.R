#!/usr/bin/Rscript

fastOH <- function(genotype) {
  #' opposite homozygous Ferdosi M. and Boerner V. (2014)
  #' Distinguished values in @param genotype are 0, and 2
  start_time <- Sys.time()
  cat("\n... Starting OH [ algorithm for fast OH by Ferdosi M. and Boerner V. (2014) ] assignment ...\n")

  # 2 maps to 1, all else maps to 0
  fpart <- + (genotype == 2)
  # 0 maps to 1, all else maps to 0
  tlpart <- + (genotype == 0)

  cat("... opposing homozygous loci counts started ...\n")
  result <- tcrossprod(fpart, tlpart)
  cat("... storing OH results in full matrix format ...\n")
  result <- result + t(result)
  cat("... OH pedigree calculated in", round(difftime(Sys.time(), start_time, units = "secs"), 2), "seconds ...\n")
  return(result)
}


generate_outfiles <- function(inpgeno, outfilestub, qc, plot = F) {
  # Setup
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc["chrset"], "--bfile", inpgeno, "--geno", qc["geno"], "--make-bed --out tmp"))
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc["chrset"], "--bfile tmp --mind", qc["mind"], "--make-bed --out tmp1"))
  system(paste("plink.exe --silent --allow-no-sex --nonfounders --chr-set", qc["chrset"], "--bfile tmp1 --maf", qc["maf"], "--hwe", qc["hwe"], "--make-bed --out tmp2"))
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc["chrset"], "--bfile tmp2 --het --out het"))

  # Calculate heterozygosity
  het <- read.table("het.het", header = T)
  het$pHET <- (het$N.NM. - het$O.HOM.) / het$N.NM.
  # # Plot histogram, found in Rplots.pdf
  if (plot) {
    hist(het$pHET, breaks = 100, col = "green", xlab = "Sample Heterozygosity", xlim = c(0.1, 0.60), main = "")
    dev.off()
  }
  # Heterozygosity limits
  less.H <- mean(het$pHET) - 6 * sd(het$pHET)
  high.H <- mean(het$pHET) + 6 * sd(het$pHET)
  # Filter out poor samples
  hetpoor <- het[which(het$pHET < less.H | het$pHET > high.H), ]
  # Write file with FID and IIDs of poor samples
  write.table(hetpoor[, 1:2], "het.poor", sep = "\t", quote = F, col.names = F, row.names = F)
  cat("...", nrow(hetpoor), "animals with poor heterozygosity ...\n")
  # Filter out samples with poor heterozygosity using file het.poor, and store files in output location
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc["chrset"], "--bfile tmp2 --remove het.poor --make-bed --out", outfilestub))
  unlink(c("het*", "nosex", "tmp*")) # Delete temporary files

  dat <- read.table(paste(outfilestub, ".bim", sep = ""))
  dat$sallele <- unique(dat[, 6])[2]
  write.table(dat[, c(2, 7)], "recodeallele.txt", quote = F, col.names = F, row.names = F)
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc["chrset"], "--bfile", outfilestub, "--thin", qc["thin"], "--recode A --recode-allele recodeallele.txt --out", outfilestub))
  unlink("recodeallele.txt")
}

assign_parent <- function(ids, parentIDs, offspringID, threshOMM, OH.PD) {
  # Lookup of only relevant IDs
  animparOH <- ids[ids$ID %in% c(parentIDs, offspringID), ]
  # Ordercodes probably point to rows in the OH.PD matrix
  offspringTested <- animparOH[animparOH$ID %in% offspringID, "ordercode"]
  # Accesses the OH.PD matrix at (x, o) where o is the current offspring and x ranges over sires and the current offspring
  OH.PDoffT <- data.frame(OPH = na.omit(OH.PD[animparOH$ordercode, offspringTested]), stringsAsFactors = F)
  OH.PDoffT$ID <- rownames(OH.PDoffT) # Insert IDs, although the data frame is already indexed by them.
  OH.PDoffT <- OH.PDoffT[which(OH.PDoffT$OPH <= threshOMM), ]

  chosen <- NA
  chosenOH <- NA
  possib <- NA
  possibOH <- NA
  if (nrow(OH.PDoffT) > 0) {
    OH.PDoffT <- OH.PDoffT[order(OH.PDoffT$OPH, decreasing = F), ]
    chosen <- OH.PDoffT$ID[1]
    chosenOH <- OH.PDoffT$OPH[1]
    if (nrow(OH.PDoffT) > 1) {
      possib <- paste(OH.PDoffT$ID[-1], collapse = "/")
      possibOH <- paste(OH.PDoffT$OPH[-1], collapse = "/")
    }
  }
  return(list(chosen = chosen, chosenOH = chosenOH, possib = possib, possibOH = possibOH))
}

plot_OH <- function(outfilestub, OH.PD) {
  OH.PDall <- OH.PD[lower.tri(OH.PD)]
  png(paste(outfilestub, ".png", sep = ""), width = 1000, height = 800)
  layout(mat = matrix(1:4, 2, ncol = 2, byrow = T))
  hist(OH.PDall, breaks = 250, xlab = "Number of opposite homozygote", main = "")
  hist(OH.PDall[which(OH.PDall < 500)], breaks = 250, xlab = "Number of opposite homozygote", main = "")
  dens <- density(OH.PDall, adjust = 2)
  plot(dens, lty = "dotted", col = "darkgreen", lwd = 2, main = "", xlab = "Number of opposite homozygote")
  dens <- density(OH.PDall[which(OH.PDall < 500)], adjust = 2)
  plot(dens, lty = "dotted", col = "darkgreen", lwd = 2, main = "", xlab = "Number of opposite homozygote")
  dev.off()
}

calculate_pedigree <- function(parentfile, outfilestub, outfolder, threshOMM) {
  dat <- read.table(paste(outfilestub, ".raw", sep = ""), skip = 1)
  cat("\n...", nrow(dat), "animals and", ncol(dat), "markers remaining for OH analysis ...\n\n")
  ids <- data.frame(ID = as.vector(dat[, 2]), ordercode = seq_len(nrow(dat)), stringsAsFactors = F)
  dat <- data.matrix(dat[, -1:-6])
  dat[is.na(dat)] <- 9
  rownames(dat) <- ids$ID

  cat("... Genotypes imported and edited ...\n")

  #### parental data with IDs and Sex of parent
  parents <- read.table(paste(parentfile, sep = ""), header = T, stringsAsFactors = F, sep = ",")
  colnames(parents) <- c("ID", "Sex")
  offspring <- ids[!ids$ID %in% parents$ID, ]

  cat("... Parental (ID, Sex) information imported ...\n")

  sires <- parents[which(parents[, 2] == "M" | parents[, 2] == "1"), "ID"]
  sires <- data.frame(ID = sires[sires %in% ids$ID], stringsAsFactors = F)
  dams <- parents[which(parents[, 2] == "F" | parents[, 2] == "2"), "ID"]
  dams <- data.frame(ID = dams[dams %in% ids$ID], stringsAsFactors = F)
  siredam <- rbind.data.frame(
    data.frame(ID = sires$ID, Sex = "M", stringsAsFactors = F),
    data.frame(ID = dams$ID, Sex = "F", stringsAsFactors = F),
    stringsAsFactors = F
  )
  write.table(siredam, paste(outfolder, "parentsafterqc.csv", sep = "/"), quote = F, row.names = F, col.names = T, sep = ",")
  cat("... total number of parents after QC", nrow(siredam), "...\n")
  cat("... with", nrow(sires), "sires and", nrow(dams), "dams ...\n")

  OH.PD <- fastOH(genotype = dat)
  diag(OH.PD) <- NA
  plot_OH(outfilestub, OH.PD)

  cat("... OH pedigree has been generated ...\n")
  pedigreconst <- data.frame(
    ID = character(),
    sire = character(), OHsire = numeric(),
    dam = character(), OHdam = numeric(),
    sirepossib = character(), OHsirepossib = character(),
    dampossib = character(), OHdampossib = character()
  )

  cat("... assigning", nrow(offspring), "offspring to", nrow(parents), "parents based on the OH counts ...\n")
  start_assign <- Sys.time()
  iterchecks.anim <- round(nrow(offspring) / 5, digits = 0)
  for (i in seq_len(nrow(offspring))) {
    father <- assign_parent(ids, sires$ID, offspring$ID[i], threshOMM, OH.PD)
    mother <- assign_parent(ids, dams$ID, offspring$ID[i], threshOMM, OH.PD)

    OHmdone <- cbind.data.frame(
      ID = offspring$ID[i],
      sire = father$chosen, OHsire = father$chosenOH,
      dam = mother$chosen, OHdam = mother$chosenOH,
      sirepossib = father$possib, OHsirepossib = father$possibOH,
      dampossib = mother$possib, OHdampossib = mother$possibOH
    )
    pedigreconst <- rbind.data.frame(pedigreconst, OHmdone, stringsAsFactors = F)

    if (i %% iterchecks.anim == 0) {
      cat("... offspring", i, "... out of", nrow(offspring), "... done\n")
    }
  }
  write.table(pedigreconst, paste(outfilestub, ".csv", sep = ""), quote = F, row.names = F, col.names = T, sep = ",")

  end_assign <- Sys.time()
  cat("... parent assignment took", round(end_assign - start_assign, 2), "seconds ...\n")
  return(pedigreconst)
}


calculate_matchchecks <- function(matchchecks, outfilestub) {
  if (!file.exists(paste(matchchecks))) {
    stop("... The file does not exist in the folder !! ...")
  }

  cat("... checking known matches !! ...\n")
  matchanims <- read.table(matchchecks, header = T, stringsAsFactors = F, sep = ",")
  colnames(matchanims) <- c("ID", "matchanims")
  matchanimsall <- data.frame(ID = as.vector(unique(c(matchanims$ID, matchanims$matchanims))), stringsAsFactors = F)

  dat <- read.table(paste(outfilestub, ".raw", sep = ""), skip = 1)
  dat <- merge(dat, matchanimsall, by.x = 2, by.y = 1)
  cat("...", nrow(dat), "animals and\n", ncol(dat), "markers remaining for OH analysis ...\n\n")
  ids <- data.frame(ID = as.vector(dat[, 1]), ordercode = seq_len(nrow(dat)), stringsAsFactors = F)
  dat <- data.matrix(dat[, -1:-6])
  dat[is.na(dat)] <- 9 # 9 symbolizes NA
  rownames(dat) <- ids$ID
  cat("... Genotypes imported and edited ...\n")

  matchanimsavail <- merge(matchanims, ids, by.x = 2, by.y = 1)
  OH.PD <- fastOH(genotype = dat)
  diag(OH.PD) <- NA

  matchacheckconst <- data.frame(ID = character(), check = character(), OHwithcheck = numeric())
  cat("... assigning offspring to parents based on the OH counts ...\n")
  iterchecks.anim <- round(nrow(matchanimsavail) / 5, digits = 0)
  for (i in seq_len(nrow(matchanimsavail))) {
    OHid <- data.frame(ID = c(matchanimsavail$matchanims[i], matchanimsavail$ID[i]), stringsAsFactors = F)
    animparOH <- ids[ids$ID %in% OHid$ID, 2]
    offsprinttested <- ids[ids$ID %in% tail(OHid$ID, 1), 2]
    OH.PDoffT <- data.frame(OPH = na.omit(OH.PD[animparOH, offsprinttested]), stringsAsFactors = F)
    OH.PDoffT$ID <- rownames(OH.PDoffT)

    OHmmatchcheck <- cbind.data.frame(
      ID = matchanimsavail$ID[i],
      check = matchanimsavail$matchanims[i],
      OHwithcheck = OH.PDoffT$OPH
    )
    matchacheckconst <- rbind.data.frame(matchacheckconst, OHmmatchcheck, stringsAsFactors = F)
    if (i %% iterchecks.anim == 0) {
      cat("...", i, "out of", nrow(matchanimsavail), "match checks ... done\n")
    }
  }
  write.table(matchacheckconst, paste(outfilestub, "match.csv", sep = ""), quote = F, row.names = F, col.names = T, sep = ",")
  return(matchacheckconst)
}


OHm <- function(inpgeno, parentfile, qc = c(geno = 0.05, mind = 0.10, maf = 0.01, hwe = 1e-6, thin = 1, chrset = 30), threshOMM = 25, matchchecks = NULL, outfilename, outfolder = ".") {
  strdate <- paste("started ...", date(), "...")

  cat("
  @*************************************************************************************@
  @              Parentage assignment based on opposing homozygotes (OH)                @
  @                                                                                     @
  @    Fast OH computation based on algorithm of Ferdosi M. and Boerner V. (2014)       @
  @                                                                    by: S.A. Boison  @
  @*************************************************************************************@
  \n\n")

  # Check for the presence of required files
  if (!file.exists("plink.exe"))
    stop("... Plink version 1.90 needed !! ...")
  if (!(file.exists(paste(inpgeno, ".bed", sep = "")) && file.exists(paste(inpgeno, ".bim", sep = "")) && file.exists(paste(inpgeno, ".fam", sep = ""))) &&
    !(file.exists(paste(inpgeno, ".map", sep = "")) && file.exists(paste(inpgeno, ".ped", sep = ""))))
    stop("... Plink genotype files needed\n The file you specified in not available !! ...")
  if (!file.exists(parentfile))
    stop("...The file you specified in not available !! ...")
  if (missing(matchchecks))
    stop("... Specify if you want to undertake match checks or not !! ...")

  outfilestub <- paste(outfolder, outfilename, sep = "/")
  generate_outfiles(inpgeno, outfilestub, qc)

  if (matchchecks == F || is.null(matchchecks) || is.na(matchchecks) || matchchecks == "") {
    resultconst <- calculate_pedigree(parentfile, outfilestub, outfolder, threshOMM)
  } else {
    resultconst <- calculate_matchchecks(matchchecks, outfilestub)
  }

  enddate <- paste("ended ...", date(), "...")
  cat("\n", strdate, "\n", enddate, "\n", sep = "")
  return(resultconst)
}