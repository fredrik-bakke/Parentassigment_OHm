#!/usr/bin/Rscript

OHm <- function(inpgeno, parentfile, qc = c(geno = 0.05, mind = 0.10, maf = 0.01, hwe = 1e-6, thin = 1, chrset = 30), threshOMM = 25, matchchecks = F, outfile) {
  cat("
  @*************************************************************************************@
  @              Parentage assignment based on opposing homozygotes (OH)                @
  @                                                                                     @
  @    Fast OH computation based on algorithm of Ferdosi M. and Boerner V. (2014)       @
  @                                                                    by: S.A. Boison  @
  @*************************************************************************************@
  \n")

  if (!file.exists("plink.exe")) {
    cat("... Plink version 1.90 needed !! ...")
    return()
  }
  if (!file.exists(paste(inpgeno, ".bim", sep = "")) | !file.exists(paste(inpgeno, ".bed", sep = "")) | !file.exists(paste(inpgeno, ".fam", sep = ""))) {
    cat("... Plink binary genotype file needed\n The file you specified in not available !! ...")
    return()
  }
  if (!file.exists(parentfile)) {
    cat("...The file you specified in not available !! ...")
    return()
  }
  if (missing(matchchecks)) {
    cat("... Specify if you want to undertake match checks or not !! ...")
    return()
  }
  cat("\n")
  strdate <- paste("started ...", date(), "...")
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc[6], "--bfile", inpgeno, "--geno", qc[1], "--make-bed --out tmp"))
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc[6], "--bfile tmp --mind", qc[2], "--make-bed --out tmp1"))
  system(paste("plink.exe --silent --allow-no-sex --nonfounders --chr-set", qc[6], "--bfile tmp1 --maf", qc[3], "--hwe", qc[4], "--make-bed --out tmp2"))
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc[6], "--bfile tmp2 --het --out het"))

  het <- read.table("het.het", header = T)
  het$pHET <- (het$N.NM. - het$O.HOM.) / het$N.NM.
  hist(het$pHET, breaks = 100, col = "green", xlab = "Sample Heterozygosity", xlim = c(0.1, 0.60), main = "")
  less.H <- (mean(het$pHET) - 6 * sd(het$pHET))
  high.H <- (mean(het$pHET) + 6 * sd(het$pHET))
  hetpoor <- het[which(het$pHET < less.H | het$pHET > high.H), ]
  write.table(hetpoor[, 1:2], "het.poor", sep = "\t", quote = F, col.names = F, row.names = F)
  cat(nrow(hetpoor), "... animals with poor heterozygosity ...\n")
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc[6], "--bfile tmp2 --remove het.poor --make-bed --out", outfile))
  unlink(c("het*", "nosex", "tmp*"))

  dat <- read.table(paste(outfile, ".bim", sep = ""))
  dat$sallele <- unique(dat[, 6])[2]
  write.table(dat[, c(2, 7)], "recodeallele.txt", quote = F, col.names = F, row.names = F)
  system(paste("plink.exe --silent --allow-no-sex --chr-set", qc[6], "--bfile", outfile, "--thin", qc[5], "--recode A --recode-allele recodeallele.txt --out", outfile))
  unlink(c("recodeallele.txt"))
  rm(het, less.H, high.H, hetpoor, dat)
  gc()

  ############ opposite homozygous Ferdosi M. and Boerner V. (2014)  ##############
  fastOH <- function(genotype) {
    cat("\n... Starting  OH [ algorithm for fast OH by Ferdosi M. and Boerner V. (2014) ] assignment ...\n")
    cm <- genotype - (floor(genotype / 9) * 8)
    fpart <- floor(cm / 2)
    lPart <- t(ceiling(((cm) - 2) / 2))
    cat("... opposing homozygous loci counts started ...\n")
    result <- (fpart %*% lPart) * (-1)
    cat("... storing OH results in full matrix format ...\n")
    result <- t(result) + result
    return(result)
  }
  ###################################################################################

  if (matchchecks == F) {
    dat <- read.table(paste(outfile, ".raw", sep = ""), skip = 1)
    cat("\n", nrow(dat), "... animals and\n", ncol(dat), " markers remaining for OH analysis ...\n\n")
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
    sires <- data.frame(sireID = sires[sires %in% ids$ID], stringsAsFactors = F)
    dams <- parents[which(parents[, 2] == "F" | parents[, 2] == "2"), "ID"]
    dams <- data.frame(damID = dams[dams %in% ids$ID], stringsAsFactors = F)
    siredam <- rbind.data.frame(data.frame(ID = sires$sireID, sex = "M", stringsAsFactors = F),
      data.frame(ID = dams$damID, sex = "F", stringsAsFactors = F),
      stringsAsFactors = F
    )
    write.table(siredam, paste("parentsafterqc.csv", sep = ""), quote = F, row.names = F, col.names = T, sep = ",")
    cat("... total number of parents after QC ", nrow(siredam), " ...\n")
    cat("... with ", nrow(sires), " sires and ", nrow(dams), " dams ...\n")

    OH.PD <- fastOH(genotype = dat)
    OH.PDall <- OH.PD[lower.tri(OH.PD)]
    diag(OH.PD) <- NA
    png(paste(outfile, ".png", sep = ""), width = 1000, height = 800)
    layout(mat = matrix(1:4, 2, ncol = 2, byrow = T))
    hist(OH.PDall, breaks = 250, xlab = "Number of opposite homozygote", main = "")
    hist(OH.PDall[which(OH.PDall < 500)], breaks = 250, xlab = "Number of opposite homozygote", main = "")
    dens <- density(OH.PDall, adjust = 2)
    plot(dens, lty = "dotted", col = "darkgreen", lwd = 2, main = "", xlab = "Number of opposite homozygote")
    dens <- density(OH.PDall[which(OH.PDall < 500)], adjust = 2)
    plot(dens, lty = "dotted", col = "darkgreen", lwd = 2, main = "", xlab = "Number of opposite homozygote")
    dev.off()

    ###
    rm(OH.PDall)
    gc()
    ###
    cat("... OH pedigree is been generated ...")
    pedigreconst <- data.frame(
      ID = character(), sire = character(), OHsire = numeric(),
      dam = character(), OHdam = numeric(), sirepossib = character(), OHsirepossib = character(),
      dampossib = character(), OHdampossib = character()
    )

    cat("\n... assigning offspring to parents based on the OH counts ...\n")
    iterchecks.anim <- round(nrow(offspring) / 5, digits = 0)
    for (i in seq_len(nrow(offspring))) {
      OHid <- data.frame(ID = c(sires$sireID, offspring$ID[i]), stringsAsFactors = F)
      animparOH <- ids[ids$ID %in% OHid$ID, ]
      offsprinttested <- animparOH[animparOH$ID %in% offspring$ID[i], "ordercode"]
      OH.PDoffT <- data.frame(OPH = na.omit(OH.PD[animparOH$ordercode, offsprinttested]), stringsAsFactors = F)
      OH.PDoffT$ID <- rownames(OH.PDoffT)
      OH.PDoffT <- OH.PDoffT[which(OH.PDoffT$OPH <= threshOMM), ]

      if (nrow(OH.PDoffT) == 0) {
        sirechosen <- NA
        sirechosenOH <- NA
        sirepossib <- NA
        sirepossibOH <- NA
      } else if (nrow(OH.PDoffT) == 1) {
        OH.PDoffT <- OH.PDoffT[order(OH.PDoffT$OPH, decreasing = F), ]
        sirechosen <- OH.PDoffT$ID[1]
        sirechosenOH <- OH.PDoffT$OPH[1]
        sirepossib <- NA
        sirepossibOH <- NA
      } else if (nrow(OH.PDoffT) > 1) {
        OH.PDoffT <- OH.PDoffT[order(OH.PDoffT$OPH, decreasing = F), ]
        sirechosen <- OH.PDoffT$ID[1]
        sirechosenOH <- OH.PDoffT$OPH[1]
        sirepossib <- paste(c(OH.PDoffT$ID[-1]), collapse = "/")
        sirepossibOH <- paste(c(OH.PDoffT$OPH[-1]), collapse = "/")
      }

      OHid <- data.frame(ID = c(dams$damID, offspring$ID[i]))
      animparOH <- ids[ids$ID %in% OHid$ID, ]
      offsprinttested <- animparOH[animparOH$ID %in% offspring$ID[i], "ordercode"]
      OH.PDoffT <- data.frame(OPH = na.omit(OH.PD[animparOH$ordercode, offsprinttested]), stringsAsFactors = F)
      OH.PDoffT$ID <- rownames(OH.PDoffT)
      OH.PDoffT <- OH.PDoffT[which(OH.PDoffT$OPH <= threshOMM), ]

      if (nrow(OH.PDoffT) == 0) {
        damchosen <- NA
        damchosenOH <- NA
        dampossib <- NA
        dampossibOH <- NA
      } else if (nrow(OH.PDoffT) == 1) {
        OH.PDoffT <- OH.PDoffT[order(OH.PDoffT$OPH, decreasing = F), ]
        damchosen <- OH.PDoffT$ID[1]
        damchosenOH <- OH.PDoffT$OPH[1]
        dampossib <- NA
        dampossibOH <- NA
      } else if (nrow(OH.PDoffT) > 1) {
        OH.PDoffT <- OH.PDoffT[order(OH.PDoffT$OPH, decreasing = F), ]
        damchosen <- OH.PDoffT$ID[1]
        damchosenOH <- OH.PDoffT$OPH[1]
        dampossib <- paste(c(OH.PDoffT$ID[-1]), collapse = "/")
        dampossibOH <- paste(c(OH.PDoffT$OPH[-1]), collapse = "/")
      }

      OHmdone <- cbind.data.frame(
        ID = offspring[i, 1], sire = sirechosen, OHsire = sirechosenOH,
        dam = damchosen, OHdam = damchosenOH, sirepossib = sirepossib,
        OHsirepossib = sirepossibOH, dampossib = dampossib, OHdampossib = dampossibOH
      )

      pedigreconst <- rbind.data.frame(pedigreconst, OHmdone, stringsAsFactors = F)

      if (i == 1) {
        write.table(OHmdone, paste(outfile, ".csv", sep = ""), quote = F, row.names = F, col.names = T, sep = ",")
      } else {
        write.table(OHmdone, paste(outfile, ".csv", sep = ""), quote = F, row.names = F, col.names = F, append = T, sep = ",")
      }
      if (i %% iterchecks.anim == 0) {
        cat("... offspring ", i, " ... out of ", nrow(offspring), " ... done\n")
      }
    }
    cat("\n")
    enddate <- paste("ends ...", date(), "...")
    cat(strdate, "\n", enddate)
    return(pedigreconst)
  } else if (matchchecks != F) {
    if (!file.exists(paste(matchchecks))) {
      cat("... The file does not exist in the folder !! ...")
      return()
    }
    cat("... checking known matches !! ...")
    matchanims <- read.table(matchchecks, header = T, stringsAsFactors = F, sep = ",")
    colnames(matchanims) <- c("ID", "matchanims")
    matchanimsall <- data.frame(ID = as.vector(unique(c(matchanims$ID, matchanims$matchanims))), stringsAsFactors = F)

    dat <- read.table(paste(outfile, ".raw", sep = ""), skip = 1)
    dat <- merge(dat, matchanimsall, by.x = 2, by.y = 1)
    cat("\n", nrow(dat), " ... animals and\n", ncol(dat), " markers remaining for OH analysis ...\n\n")
    ids <- data.frame(ID = as.vector(dat[, 1]), ordercode = seq_len(nrow(dat)), stringsAsFactors = F)
    dat <- data.matrix(dat[, -1:-6])
    dat[is.na(dat)] <- 9
    rownames(dat) <- ids$ID
    cat("... Genotypes imported and edited ...\n")

    matchanimsavail <- merge(matchanims, ids, by.x = 2, by.y = 1)
    OH.PD <- fastOH(genotype = dat)
    diag(OH.PD) <- NA
    cat("... OH pedigree is been generated ...")

    matchacheckconst <- data.frame(ID = character(), check = character(), OHwithcheck = numeric())
    cat("\n... assigning offspring to parents based on the OH counts ...\n")
    iterchecks.anim <- round(nrow(matchanimsavail) / 5, digits = 0)
    for (i in seq_len(nrow(matchanimsavail))) {
      OHid <- data.frame(ID = c(matchanimsavail$matchanims[i], matchanimsavail$ID[i]), stringsAsFactors = F)
      animparOH <- ids[ids$ID %in% OHid$ID, 2]
      offsprinttested <- ids[ids$ID %in% tail(OHid$ID, 1), 2]
      OH.PDoffT <- data.frame(OPH = na.omit(OH.PD[animparOH, offsprinttested]), stringsAsFactors = F)
      OH.PDoffT$ID <- rownames(OH.PDoffT)

      OHmmatchcheck <- cbind.data.frame(
        ID = matchanimsavail$ID[i], check = matchanimsavail$matchanims[i],
        OHwithcheck = OH.PDoffT$OPH
      )
      matchacheckconst <- rbind.data.frame(matchacheckconst, OHmmatchcheck, stringsAsFactors = F)
      if (i %% iterchecks.anim == 0) {
        cat("... ", i, " ... out of ", nrow(matchanimsavail), " match checks ... done\n")
      }
    }
    cat("\n")
    enddate <- paste("ends ...", date(), "...")
    cat(strdate, "\n", enddate)
    write.table(matchacheckconst, paste(outfile, "match.csv", sep = ""), quote = F, row.names = F, col.names = T, sep = ",")
    return(matchacheckconst)
  }
}