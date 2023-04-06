### Setup ###
rm(list = ls())
library(OrgMassSpecR) # The OrgMassSpecR version must be >= 0.5-1.

# Make directory structure including each step
cat("* Creating directory structure as needed...")
dir.create("Step 1 Input Spectra", showWarnings = FALSE)
dir.create("Step 2 Input Tables", showWarnings = FALSE)
dir.create("Step 3 Split MSP Files", showWarnings = FALSE)
dir.create("Step 4 MSP to CSV", showWarnings = FALSE)
dir.create("Step 5 Halogenation Filter", showWarnings = FALSE)
dir.create("Step 6 Additional Rules", showWarnings = FALSE)
dir.create("Step 7 Final Results", showWarnings = FALSE)
cat("OK\n")

# Check that theoretical distributions are present. 
# cat("* Checking for presence of Theoretical Distributions.csv...")
# if(file.exists("Theoretical Distributions.csv")) {
#   cat("OK\n") } else {
#     cat("ERROR\n")
#     stop("Working directory must contain the file Theoretical Distributions.csv")
# }

# Check if input spectra and input tables have the same name.
cat("* Checking if input spectra and input tables have the same name...")
step1Files <- sort(sub(".txt", "", list.files("Step 1 Input Spectra")))
step2Files <- sort(sub(".csv", "", list.files("Step 2 Input Tables")))
if(all(step1Files == step2Files)) {
  cat("OK\n") } else {
  cat("ERROR\n")
    stop("Input spectra (*.txt) in Step 1 folder and input 
       tables (*.csv) in Step 2 folder must have the same names.
       Recommend manually deleting folders for Steps 3-7 before re-running.") 
  }

### Function: Split MSP spectra into individual files ###
SplitMSP <- function(originFileName, outputFolderName) {
  
  a <- readLines(originFileName)
  
  # Make factor specifying each spectrum.
  spectrumFactor <- vector(mode = "integer", length = length(a)) 
  spectrumCount <- 1
  # Variable spectrumCount is the same for each line of the spectrum, then iterates
  # +1 if encounters a newline (defining the next spectrum).
  for(i in 1: length(a)) {
    if(a[i] != "") {
      spectrumFactor[i] <- spectrumCount
    } else {
      spectrumFactor[i] <- spectrumCount
      spectrumCount <- spectrumCount + 1
    }
  }
  
  # TODO: If end of file contains more than one newline, ignore.

  # Split vector `a` containing the spectra by the defined spectrumFactor.
  spectrumList <- split(a, f = spectrumFactor)
  
  # Make individual files, where the filename specifies the original order in 
  # the input file. Removes blank lines and writeLines adds a single blank line
  # to the bottom of the file.
  for(j in 1:length(spectrumList)) {
    individualSpectrum <- spectrumList[[j]]
    individualSpectrum <- individualSpectrum[individualSpectrum != ""]
    fullFilename <- paste(outputFolderName, "/Spectrum ", j, ".txt", sep = "")
    writeLines(individualSpectrum, fullFilename)
    # Reporting
    cat('Step 3: Splitting ', fullFilename, '...\n', sep = "")
  }
}

## ----Split spectra-------------------------------------------------------
cat("* Beginning Step 3...\n")
tmp <- list.files("Step 1 Input Spectra", full.names = TRUE)

for (i in 1:length(tmp)) {
  # Capture MSP file name to use as folder name for split files.
  f1 <- strsplit(tmp[i], split = c("/"), fixed = TRUE)[[1]][2]
  folder <- substr(f1, 1, nchar(f1) - 4) # remove .txt
  dir.create(paste("Step 3 Split MSP Files/", folder, " Spectra", sep = ""))
  SplitMSP(originFileName = tmp[i], 
           outputFolderName = paste("Step 3 Split MSP Files/", folder, " Spectra", sep = ""))
}

## ----Function to convert MSP to CSV--------------------------------------
rm(list = ls())

ConvertMSP <- function(file) {
  x <- ReadMspFile(file, skip = 2, comment.char = "", remove.placeholders = FALSE)
  # If reporting is needed:
  # cat('Making ', file, '\n', sep = "")
  f1 <- strsplit(file, split = c("/"), fixed = TRUE)[[1]][3] # filename
  f2 <- strsplit(file, split = c("/"), fixed = TRUE)[[1]][2] # sample directory
  fileName <- substr(f1, 1, nchar(f1) - 4) # remove .txt
  cat('Step 4: Converting ', fileName, '...\n', sep = "")
  write.csv(x, file = paste("Step 4 MSP to CSV/", f2, "/", fileName, ".csv", sep = ''), row.names = FALSE)
}

## ----Convert MSP to CSV--------------------------------------------------
cat("* Beginning Step 4...\n")
tmp <- list.dirs("Step 3 Split MSP Files", recursive = FALSE)

for (i in 1:length(tmp)) {
  inputSpectraDirectory <- tmp[i]
  spectra <- dir(inputSpectraDirectory)
  f <- strsplit(tmp[i], split = c("/"), fixed = TRUE)[[1]][2] # sample directory
  dir.create(paste("Step 4 MSP to CSV/", f, sep = ""))
  x <- do.call('rbind', lapply(paste(inputSpectraDirectory, '/', spectra, sep = ''), ConvertMSP))
}

## ----Filter functions----------------------------------------------------

# The isotope ratio filter.
ProcessSpectrum <- function(x) {
  
  # Normalize peak intensity to percentage of max peak and only examine m/z > 100 amu.
  x$percentIntensity <- with(x, intensity / max(intensity) * 100)
  x$Filter <- FALSE
  x <- x[x$mz > 100, ]
  
  for(i in (nrow(x) - 6):5) {
    
    if(x$percentIntensity[i] > 3) {
      
      # Calculate peak ratios.
      
      peakRatio1 <- x$percentIntensity[i + 1] / x$percentIntensity[i]
      peakRatio2 <- x$percentIntensity[i + 2] / x$percentIntensity[i]
      peakRatio3 <- x$percentIntensity[i + 3] / x$percentIntensity[i]
      peakRatio4 <- x$percentIntensity[i + 4] / x$percentIntensity[i]
      peakRatio5 <- x$percentIntensity[i + 5] / x$percentIntensity[i]
      peakRatio6 <- x$percentIntensity[i + 6] / x$percentIntensity[i]
      peakRatioa <- x$percentIntensity[i - 1] / x$percentIntensity[i]
      peakRatiob <- x$percentIntensity[i - 2] / x$percentIntensity[i]
      peakRatioc <- x$percentIntensity[i - 4] / x$percentIntensity[i]
      
      # Determine if peak ratios pass rules.
      
      if(peakRatio1 < 0.5 & peakRatio1 > 0.003 &   # Rule 1
         peakRatio2 > 0.225 & peakRatio2 <= 1 &    # Rule 2
         peakRatio1 < peakRatio2 &                 # Rule 3
         peakRatio3 < 0.5 & peakRatio3 > 0.001 &   # Rule 4
         peakRatio3 < peakRatio1 &                 # Rule 5
         peakRatio4 < peakRatio2 &                 # Rule 6
         peakRatio5 < peakRatio1 &                 # Rule 7
         peakRatio6 < peakRatio2 &                 # Rule 8
         peakRatioa < 0.5 &                        # Rule 9
         peakRatiob < 1 &                          # Rule 10
         peakRatioc < 1) {                         # Rule 11
        
        x$Filter[i] <- TRUE
        
      }
      
    } # encloses filter
    
  } # encloses for loop that iterates through each m/z
  
  return(any(x$Filter))  # Returns TRUE if compound is halogenated.
  
}

# Prepares filenames and calls ProcessSpectrum

FindHalogens <- function(spectrumFilename) {
  
  cat('Step 5: Processing spectrum', spectrumFilename, '...\n')
  
  x <- read.csv(spectrumFilename, stringsAsFactors = FALSE)
  
  # Shorten filename so it is the same for all spectra (used later).
  fileName1 <- strsplit(spectrumFilename, split = c("/"), fixed = TRUE)[[1]][3]
  fileName2 <- substr(fileName1, 1, nchar(fileName1) - 4) # remove extension
  
  y <- ProcessSpectrum(x)
  
  return(data.frame(SpectrumFilename = fileName2, Halogenated = y))
  
}

# Calls FindHalogens for all spectra in sample and merges table information (retention times, etc.)

ProcessHalogenationSearch <- function(spectraDirectory, inputTableName, sampleName, outputTableName) {
  
  spectra <- dir(spectraDirectory)
  x <- do.call('rbind', lapply(paste(spectraDirectory, '/', spectra, sep = ''), FindHalogens))
  chromaTOF_table <- read.csv(inputTableName, stringsAsFactors = FALSE)
  chromaTOF_table$SpectrumFilename <- paste("Spectrum", chromaTOF_table$PeakNumber)
  x <- merge(x, chromaTOF_table, by = 'SpectrumFilename', all = TRUE)
  x$Sample <- sampleName
  write.csv(x, file = outputTableName, row.names = FALSE)
  
}

## ----Apply filter--------------------------------------------------------
cat("* Beginning Step 5...\n")
tmp <- dir("Step 4 MSP to CSV") # "Sample 1 Spectra" "Sample 2 Spectra"
for (i in 1:length(tmp)) {
  sn <- sub(" Spectra", "", x = tmp[i], fixed = TRUE)
  ProcessHalogenationSearch(spectraDirectory = paste("Step 4 MSP to CSV/", tmp[i], sep = ""),
                            inputTableName = paste("Step 2 Input Tables/", sn, ".csv", sep = ""),
                            sampleName = sn,
                            outputTableName = paste("Step 5 Halogenation Filter/Step 5 ", sn, ".csv", sep = ""))
}

## ----Additional rule functions-------------------------------------------

ApplyFilter <- function(x) {

  x$Filter <- FALSE
  x <- x[x$mz > 100, ]

  for(i in (nrow(x) - 6):5) {

    if(x$percentIntensity[i] > 3) {

      # Calculate peak ratios.
      
      peakRatio1 <- x$percentIntensity[i + 1] / x$percentIntensity[i]
      peakRatio2 <- x$percentIntensity[i + 2] / x$percentIntensity[i]
      peakRatio3 <- x$percentIntensity[i + 3] / x$percentIntensity[i]
      peakRatio4 <- x$percentIntensity[i + 4] / x$percentIntensity[i]
      peakRatio5 <- x$percentIntensity[i + 5] / x$percentIntensity[i]
      peakRatio6 <- x$percentIntensity[i + 6] / x$percentIntensity[i]
      peakRatioa <- x$percentIntensity[i - 1] / x$percentIntensity[i]
      peakRatiob <- x$percentIntensity[i - 2] / x$percentIntensity[i]
      peakRatioc <- x$percentIntensity[i - 4] / x$percentIntensity[i]
      
      # Determine if peak ratios pass rules.
      
      if(peakRatio1 < 0.5 & peakRatio1 > 0.003 &   # Rule 1
         peakRatio2 > 0.225 & peakRatio2 <= 1 &    # Rule 2
         peakRatio1 < peakRatio2 &                 # Rule 3
         peakRatio3 < 0.5 & peakRatio3 > 0.001 &   # Rule 4
         peakRatio3 < peakRatio1 &                 # Rule 5
         peakRatio4 < peakRatio2 &                 # Rule 6
         peakRatio5 < peakRatio1 &                 # Rule 7
         peakRatio6 < peakRatio2 &                 # Rule 8
         peakRatioa < 0.5 &                        # Rule 9
         peakRatiob < 1 &                          # Rule 10
         peakRatioc < 1) {                         # Rule 11

        x$Filter[i] <- TRUE

      }

    }

  }

  return(x)

}

# Select (already identified) halogenated distributions. Remove miss-identified
# peaks within the same distribution. Select 5 most abundant distributions. 
# Select m/z window for the distributions.

SelectHalogenatedDistribution <- function(xFilter, x) { 
  
  # Remove miss-identified peaks within the same distribution. Sometimes another
  # peak within the same distribution is identified. Only the maximum peak 
  # within the distribution should be identified. The following code eliminates
  # any multiple identifications within the same distribution, although not
  # necessairly the maximum one. The maximum peak is identified again later when
  # plotting the distribution.
  
  xFilter$Diff <- c(diff(xFilter$mz), 100) # 100 is placeholder at the end
  xFilter <- xFilter[xFilter$Diff > 10, ] # exclude if within 10 mz units
  
  # Select the 5 most abundant distributions. Then return to the original order.
  if(nrow(xFilter) > 5) {
    
    xFilter <- head(xFilter[rev(order(xFilter$intensity)), ], n = 5) 
    xFilter <- xFilter[order(xFilter$mz), ]
    
  }
  
  # Select window of m/z around the experimental distributions and store in
  # distList.
  
  distList <- vector(mode = 'list')
  for (i in 1:nrow(xFilter)) {
    distList[[i]] <- x[x$mz %in% seq(from = xFilter$mz[i] - 10, to = xFilter$mz[i] + 10, by = 1), ]  
  }
  return(distList)
  
}

# Function to determine matching theoretical distributions. Inputs are
# distList[[i]] and tDistList

SimilarityScoreMatch <- function(xDist, tDistList, tDistAll) {
  
  resultStore <- data.frame(NULL)
  
  # Iterate through each theoretical distribution and calculate a similarity
  # score against the experimental spectrum.
  for (j in 1:length(tDistList)) {
    
    tDist <- tDistAll[tDistAll$DistributionLabel == tDistList[j], ] # current theoretical distribution
    
    # Align theoretical distribution to experimental distribution based on the
    # maximum peak (maxvalue). The theoretical distribution m/z window is
    # given the same size as the experimental m/z window.
    tDist$mz[tDist$percent == 100] <- xDist$mz[which.max(xDist$percentIntensity)]
    maxvalue <- xDist$mz[which.max(xDist$percentIntensity)]
    index <- which.max(tDist$mz)
    
    if(index > 1) {
      startValue <- maxvalue - (index - 1)
      endValue <- startValue + nrow(tDist) - 1
      tDist$mz <- seq(from = startValue, to = endValue, by = 1)
    }
    
    # Need special case when the max ion is the first one.
    if(index == 1) {
      startValue <- maxvalue
      endValue <- maxvalue + nrow(tDist) - 1
      tDist$mz <- seq(from = startValue, to = endValue, by = 1)
    }
    
    # Merge experimental and theoretical distributions by the m/z.
    cDist <- merge(xDist, tDist, by = 'mz', all.x = TRUE, all.y = FALSE)
    cDist$percent[is.na(cDist$percent)] <- 0
    cDist$DistributionLabel <- tDistList[j]
    
    # Calculate similarity score, add to dataframe, and store in resultStore.
    u <- (cDist$intensity.x/max(cDist$intensity.x)) * 100
    v <- cDist$percent
    cDist$SimilarityScore <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
    resultStore <- rbind(resultStore, cDist)
    
  }
  
  return(resultStore)
  
}

Report <- function(spectrumFilename) {
  cat('Step 6: Making table for', spectrumFilename, '...\n')
  x <- read.csv(spectrumFilename, stringsAsFactors = FALSE)
  x$percentIntensity <- with(x, intensity / max(intensity) * 100)
  
  # Identify halogenated distributions
  x <- ApplyFilter(x = x)
  
  # Isolate halogenated distributions. xFilter contains only ions identifed as
  # halogenated.
  xFilter <- x[x$Filter == TRUE, ]
  
  # Need to skip if no halogenation was detected.
  if(nrow(xFilter) > 0) {
    distList <- SelectHalogenatedDistribution(xFilter, x)
    
    # Iterate through experimental distributions.
    for (i in 1:length(distList)) {
      xDist <- distList[[i]]
    }
    
    # Identify Distributions -------------------------------------------
    
    tDistAll <- structure(list(mz = c(234, 235.01, 236, 237.01, 238.01, 311.91, 
312.92, 313.91, 314.92, 315.91, 316.91, 317.92, 389.83, 390.83,
391.82, 392.83, 393.82, 394.82, 395.82, 396.82, 397.83, 398.83,
467.74, 468.74, 469.73, 470.74, 471.73, 472.74, 473.73, 474.73, 
475.73, 476.73, 477.73, 545.65, 546.65, 547.64, 548.65, 549.64,
550.65, 551.64, 552.64, 553.64, 554.64, 555.64, 556.64, 557.64,
623.56, 624.56, 625.55, 626.56, 627.55, 628.56, 629.55, 630.55,
631.55, 632.55, 633.55, 634.55, 635.55, 636.55, 637.55, 701.47,
702.47, 703.47, 704.47, 705.46, 706.47, 707.46, 708.46, 709.46,
710.46, 711.46, 712.46, 713.46, 714.46, 715.45, 716.46, 779.38,
780.38, 781.38, 782.38, 783.37, 784.38, 785.37, 786.38, 787.37,
788.37, 789.37, 790.37, 791.37, 792.37, 793.36, 794.37, 795.36,
796.36, 857.29, 858.29, 859.29, 860.29, 861.28, 862.29, 863.28,
864.29, 865.28, 866.28, 867.28, 868.28, 869.28, 870.28, 871.27,
872.28, 873.27, 874.28, 875.27, 877.28, 935.2, 936.2, 937.2,
938.2, 939.19, 940.2, 941.19, 942.2, 943.19, 944.19, 945.19,
946.19, 947.19, 948.19, 949.18, 950.19, 951.18, 952.19, 953.18, 
954.18, 955.18, 956.18, 190.05, 191.06, 192.05, 193.06, 194.06,
224.02, 225.02, 226.01, 227.02, 228.01, 229.01, 230.02, 231.02,
257.98, 258.98, 259.97, 260.98, 261.97, 262.97, 263.97, 264.97,
265.97, 291.94, 292.94, 293.94, 294.94, 295.93, 296.94, 297.93,
298.93, 299.93, 300.93, 301.93, 325.9, 326.9, 327.9, 328.9, 329.89,
330.9, 331.89, 332.89, 333.89, 334.89, 335.88, 359.86, 360.86,
361.86, 362.86, 363.85, 364.86, 365.85, 366.85, 367.85, 368.85,
369.85, 370.85, 371.84, 393.82, 394.82, 395.82, 396.82, 397.82,
398.82, 399.81, 400.82, 401.81, 402.81, 403.81, 404.81, 405.8,
406.82, 427.78, 428.79, 429.78, 430.78, 431.78, 432.78, 433.77,
434.78, 435.77, 436.77, 437.77, 438.77, 439.77, 440.77, 441.76,
461.74, 462.75, 463.74, 464.74, 465.74, 466.74, 467.73, 468.74,
469.73, 470.73, 471.73, 472.73, 473.73, 474.73, 475.72, 495.7, 
496.71, 497.7, 498.7, 499.7, 500.7, 501.7, 502.7, 503.69, 504.7,
505.69, 506.69, 507.69, 508.69, 509.68, 510.69, 511.68, 513.68,
267.97, 268.97, 269.96, 270.97, 271.96, 272.96, 273.97, 274.97,
301.93, 302.93, 303.92, 304.93, 305.92, 306.92, 307.92, 308.92,
309.93, 335.89, 336.89, 337.89, 338.89, 339.88, 340.89, 341.88,
342.88, 343.88, 344.88, 369.85, 370.85, 371.85, 372.85, 373.84,
374.85, 375.84, 376.84, 377.84, 378.84, 379.84, 380.84, 403.81,
404.81, 405.81, 406.81, 407.8, 408.81, 409.8, 410.8, 411.8, 412.8,
413.8, 414.8, 415.79, 416.8, 345.88, 346.88, 347.87, 348.88,
349.87, 350.87, 351.87, 352.87, 353.88, 379.84, 380.84, 381.83,
382.84, 383.83, 384.84, 385.83, 386.83, 387.83, 388.83, 389.83,
413.8, 414.8, 415.8, 416.8, 417.79, 418.8, 419.79, 420.79, 421.79,
422.79, 423.79, 424.79, 425.79, 447.76, 448.76, 449.76, 450.76,
451.75, 452.76, 453.75, 454.75, 455.75, 456.75, 457.75, 458.75,
459.74, 461.75, 481.72, 482.72, 483.72, 484.72, 485.72, 486.72,
487.71, 488.72, 489.71, 490.71, 491.71, 492.71, 493.7, 494.71, 
495.7, 423.79, 424.79, 425.78, 426.79, 427.78, 428.79, 429.78,
430.78, 431.78, 432.78, 433.78, 457.75, 458.75, 459.75, 460.75,
461.74, 462.75, 463.74, 464.74, 465.74, 466.74, 467.74, 468.74,
469.74, 491.71, 492.71, 493.71, 494.71, 495.7, 496.71, 497.7,
498.7, 499.7, 500.7, 501.7, 502.7, 503.69, 504.7, 525.67, 526.67,
527.67, 528.67, 529.66, 530.67, 531.66, 532.67, 533.66, 534.66,
535.66, 536.66, 537.65, 538.66, 539.65, 559.63, 560.63, 561.63,
562.63, 563.63, 564.63, 565.62, 566.63, 567.62, 568.62, 569.62,
570.62, 571.62, 572.62, 573.61, 574.63, 575.61, 501.7, 502.7,
503.69, 504.7, 505.69, 506.7, 507.69, 508.69, 509.69, 510.69, 
511.69, 512.69, 513.69, 535.66, 536.66, 537.66, 538.66, 539.65,
540.66, 541.65, 542.65, 543.65, 544.65, 545.65, 546.65, 547.64,
548.65, 569.62, 570.62, 571.62, 572.62, 573.61, 574.62, 575.61,
576.62, 577.61, 578.61, 579.61, 580.61, 581.61, 582.61, 583.6,
584.61, 603.58, 604.58, 605.58, 606.58, 607.58, 608.58, 609.57,
610.58, 611.57, 612.57, 613.57, 614.57, 615.57, 616.57, 617.56,
618.57, 619.56, 620.56, 637.54, 638.54, 639.54, 640.54, 641.54,
642.54, 643.53, 644.54, 645.53, 646.53, 647.53, 648.53, 649.53,
650.53, 651.52, 652.53, 653.52, 654.52, 579.61, 580.61, 581.61,
582.61, 583.6, 584.61, 585.6, 586.6, 587.6, 588.6, 589.6, 590.6,
591.6, 592.6, 613.57, 614.57, 615.57, 616.57, 617.56, 618.57,
619.56, 620.57, 621.56, 622.56, 623.56, 624.56, 625.56, 626.56,
627.55, 628.56, 629.56, 647.53, 648.53, 649.53, 650.53, 651.53,
652.53, 653.52, 654.53, 655.52, 656.52, 657.52, 658.52, 659.52, 
660.52, 661.51, 662.52, 663.51, 681.49, 682.49, 683.49, 684.49,
685.49, 686.49, 687.48, 688.49, 689.48, 690.48, 691.48, 692.48,
693.48, 694.48, 695.47, 696.48, 697.47, 698.47, 715.45, 716.45,
717.45, 718.45, 719.45, 720.45, 721.44, 722.45, 723.44, 724.45, 
725.44, 726.44, 727.44, 728.44, 729.44, 730.44, 731.43, 732.43,
733.43, 734.43), intensity = c(4480L, 599L, 4342L, 553L, 26L,
2203L, 319L, 4410L, 584L, 2196L, 270L, 18L, 1155L, 164L, 3330L,
425L, 3291L, 425L, 1067L, 133L, 9L, 1L, 615L, 76L, 2221L, 292L,
3309L, 432L, 2219L, 306L, 465L, 62L, 3L, 288L, 32L, 1422L, 167L, 
2813L, 357L, 2756L, 327L, 1353L, 187L, 260L, 34L, 4L, 124L, 23L,
879L, 122L, 2150L, 267L, 2783L, 359L, 1990L, 267L, 774L, 109L,
139L, 13L, 1L, 66L, 7L, 491L, 67L, 1542L, 171L, 2421L, 285L, 
2465L, 261L, 1437L, 207L, 453L, 44L, 73L, 10L, 34L, 4L, 308L,
29L, 1023L, 136L, 2001L, 233L, 2386L, 299L, 1859L, 254L, 958L,
145L, 255L, 37L, 35L, 4L, 18L, 2L, 190L, 20L, 648L, 81L, 1524L,
197L, 2193L, 267L, 2152L, 264L, 1450L, 164L, 588L, 84L, 120L, 
14L, 23L, 1L, 10L, 2L, 89L, 15L, 436L, 57L, 1115L, 129L, 1893L,
234L, 2155L, 288L, 1725L, 233L, 961L, 138L, 373L, 40L, 87L, 13L,
6L, 1L, 6588L, 867L, 2276L, 258L, 11L, 5106L, 637L, 3249L, 437L,
489L, 75L, 6L, 1L, 3791L, 535L, 3738L, 460L, 1177L, 177L, 111L,
8L, 3L, 2895L, 407L, 3626L, 457L, 1867L, 242L, 424L, 48L, 31L,
2L, 1L, 2113L, 304L, 3536L, 485L, 2323L, 298L, 701L, 108L, 101L,
20L, 11L, 1678L, 224L, 3180L, 405L, 2612L, 320L, 1099L, 138L,
261L, 31L, 45L, 5L, 2L, 1236L, 183L, 2801L, 373L, 2755L, 333L,
1492L, 172L, 479L, 70L, 86L, 8L, 11L, 1L, 915L, 130L, 2494L,
324L, 2743L, 361L, 1781L, 240L, 707L, 92L, 164L, 21L, 23L, 4L,
1L, 702L, 98L, 2068L, 280L, 2615L, 323L, 2114L, 267L, 962L, 147L,
286L, 57L, 66L, 3L, 12L, 496L, 72L, 1772L, 231L, 2503L, 353L,
2208L, 288L, 1207L, 165L, 492L, 63L, 100L, 18L, 24L, 4L, 3L, 
1L, 3351L, 441L, 4393L, 571L, 1103L, 131L, 9L, 1L, 2579L, 342L,
4129L, 533L, 1872L, 245L, 265L, 31L, 4L, 1888L, 251L, 3783L,
506L, 2481L, 303L, 652L, 83L, 44L, 9L, 1512L, 211L, 3308L, 468L,
2737L, 308L, 1064L, 140L, 197L, 33L, 17L, 5L, 1117L, 143L, 2894L,
368L, 2950L, 363L, 1448L, 171L, 400L, 52L, 80L, 9L, 4L, 1L, 1692L,
234L, 3847L, 494L, 2779L, 349L, 531L, 70L, 4L, 1301L, 193L, 3388L,
443L, 2903L, 391L, 1084L, 144L, 132L, 20L, 1L, 988L, 127L, 2837L,
359L, 3056L, 423L, 1522L, 216L, 399L, 46L, 23L, 3L, 1L, 748L,
105L, 2416L, 318L, 2981L, 386L, 1915L, 248L, 680L, 75L, 108L,
12L, 7L, 1L, 574L, 71L, 1957L, 267L, 2876L, 365L, 2159L, 313L,
965L, 131L, 248L, 26L, 41L, 3L, 4L, 854L, 108L, 2838L, 374L,
3349L, 427L, 1561L, 210L, 254L, 24L, 1L, 657L, 93L, 2376L, 282L,
3187L, 426L, 1976L, 262L, 603L, 75L, 56L, 6L, 1L, 517L, 72L,
1898L, 213L, 3029L, 396L, 2345L, 277L, 876L, 145L, 197L, 21L, 
13L, 1L, 410L, 60L, 1492L, 196L, 2701L, 354L, 2475L, 315L, 1316L,
161L, 384L, 62L, 59L, 8L, 7L, 309L, 39L, 1289L, 192L, 2434L,
311L, 2511L, 316L, 1517L, 207L, 611L, 88L, 132L, 23L, 19L, 1L,
1L, 472L, 64L, 1835L, 242L, 3012L, 409L, 2411L, 328L, 971L, 129L,
104L, 22L, 1L, 329L, 50L, 1527L, 184L, 2792L, 338L, 2639L, 299L,
1291L, 157L, 315L, 45L, 32L, 2L, 251L, 33L, 1256L, 160L, 2483L,
326L, 2591L, 361L, 1587L, 202L, 547L, 79L, 99L, 14L, 10L, 1L,
175L, 33L, 961L, 148L, 2161L, 283L, 2587L, 341L, 1848L, 228L,
850L, 93L, 219L, 27L, 39L, 4L, 2L, 1L, 159L, 15L, 742L, 104L,
1896L, 232L, 2566L, 320L, 1969L, 268L, 1086L, 111L, 376L, 59L,
80L, 5L, 11L, 1L, 199L, 30L, 1191L, 144L, 2453L, 313L, 2784L,
359L, 1679L, 206L, 514L, 59L, 65L, 4L, 209L, 26L, 906L, 121L,
2158L, 283L, 2604L, 330L, 2032L, 229L, 783L, 99L, 173L, 28L,
17L, 1L, 1L, 121L, 18L, 726L, 87L, 1886L, 271L, 2561L, 334L,
2060L, 279L, 1103L, 136L, 299L, 40L, 64L, 11L, 4L, 98L, 14L, 
557L, 76L, 1613L, 219L, 2362L, 317L, 2239L, 291L, 1321L, 177L,
504L, 76L, 96L, 18L, 19L, 3L, 66L, 10L, 452L, 60L, 1395L, 160L,
2179L, 287L, 2238L, 313L, 1534L, 202L, 709L, 78L, 231L, 33L,
46L, 2L, 3L, 2L), percent = c(100, 13.37, 96.92, 12.34, 0.58,
49.95, 7.23, 100, 13.24, 49.8, 6.12, 0.41, 34.68, 4.92, 100,
12.76, 98.83, 12.76, 32.04, 3.99, 0.27, 0.03, 18.59, 2.3, 67.12,
8.82, 100, 13.06, 67.06, 9.25, 14.05, 1.87, 0.09, 10.24, 1.14,
50.55, 5.94, 100, 12.69, 97.97, 11.62, 48.1, 6.65, 9.24, 1.21,
0.14, 4.46, 0.83, 31.58, 4.38, 77.25, 9.59, 100, 12.9, 71.51,
9.59, 27.81, 3.92, 4.99, 0.47, 0.04, 2.68, 0.28, 19.92, 2.72,
62.56, 6.94, 98.22, 11.56, 100, 10.59, 58.3, 8.4, 18.38, 1.78,
2.96, 0.41, 1.42, 0.17, 12.91, 1.22, 42.88, 5.7, 83.86, 9.77,
100, 12.53, 77.91, 10.65, 40.15, 6.08, 10.69, 1.55, 1.47, 0.17,
0.82, 0.09, 8.66, 0.91, 29.55, 3.69, 69.49, 8.98, 100, 12.18,
98.13, 12.04, 66.12, 7.48, 26.81, 3.83, 5.47, 0.64, 1.05, 0.05, 
0.46, 0.09, 4.13, 0.7, 20.23, 2.65, 51.74, 5.99, 87.84, 10.86,
100, 13.36, 80.05, 10.81, 44.59, 6.4, 17.31, 1.86, 4.04, 0.6,
0.28, 0.05, 100, 13.16, 34.55, 3.92, 0.17, 100, 12.48, 63.63,
8.56, 9.58, 1.47, 0.12, 0.02, 100, 14.11, 98.6, 12.13, 31.05,
4.67, 2.93, 0.21, 0.08, 79.84, 11.22, 100, 12.6, 51.49, 6.67,
11.69, 1.32, 0.85, 0.06, 0.03, 59.76, 8.6, 100, 13.72, 65.7,
8.43, 19.82, 3.05, 2.86, 0.57, 0.31, 52.77, 7.04, 100, 12.74,
82.14, 10.06, 34.56, 4.34, 8.21, 0.97, 1.42, 0.16, 0.06, 44.13,
6.53, 100, 13.32, 98.36, 11.89, 53.27, 6.14, 17.1, 2.5, 3.07,
0.29, 0.39, 0.04, 33.36, 4.74, 90.92, 11.81, 100, 13.16, 64.93,
8.75, 25.77, 3.35, 5.98, 0.77, 0.84, 0.15, 0.04, 26.85, 3.75,
79.08, 10.71, 100, 12.35, 80.84, 10.21, 36.79, 5.62, 10.94, 2.18,
2.52, 0.11, 0.46, 19.82, 2.88, 70.8, 9.23, 100, 14.1, 88.21,
11.51, 48.22, 6.59, 19.66, 2.52, 4, 0.72, 0.96, 0.16, 0.12, 0.04,
76.28, 10.04, 100, 13, 25.11, 2.98, 0.2, 0.02, 62.46, 8.28, 100,
12.91, 45.34, 5.93, 6.42, 0.75, 0.1, 49.91, 6.63, 100, 13.38, 
65.58, 8.01, 17.23, 2.19, 1.16, 0.24, 45.71, 6.38, 100, 14.15,
82.74, 9.31, 32.16, 4.23, 5.96, 1, 0.51, 0.15, 37.86, 4.85, 98.1,
12.47, 100, 12.31, 49.08, 5.8, 13.56, 1.76, 2.71, 0.31, 0.14,
0.03, 43.98, 6.08, 100, 12.84, 72.24, 9.07, 13.8, 1.82, 0.1,
38.4, 5.7, 100, 13.08, 85.68, 11.54, 32, 4.25, 3.9, 0.59, 0.03,
32.33, 4.16, 92.83, 11.75, 100, 13.84, 49.8, 7.07, 13.06, 1.51,
0.75, 0.1, 0.03, 25.09, 3.52, 81.05, 10.67, 100, 12.95, 64.24,
8.32, 22.81, 2.52, 3.62, 0.4, 0.23, 0.03, 19.96, 2.47, 68.05,
9.28, 100, 12.69, 75.07, 10.88, 33.55, 4.55, 8.62, 0.9, 1.43,
0.1, 0.14, 25.5, 3.22, 84.74, 11.17, 100, 12.75, 46.61, 6.27,
7.58, 0.72, 0.03, 20.61, 2.92, 74.55, 8.85, 100, 13.37, 62, 8.22,
18.92, 2.35, 1.76, 0.19, 0.03, 17.07, 2.38, 62.66, 7.03, 100,
13.07, 77.42, 9.14, 28.92, 4.79, 6.5, 0.69, 0.43, 0.03, 15.18,
2.22, 55.24, 7.26, 100, 13.11, 91.63, 11.66, 48.72, 5.96, 14.22, 
2.3, 2.18, 0.3, 0.26, 12.31, 1.55, 51.33, 7.65, 96.93, 12.39,
100, 12.58, 60.41, 8.24, 24.33, 3.5, 5.26, 0.92, 0.76, 0.04,
0.04, 15.67, 2.12, 60.92, 8.03, 100, 13.58, 80.05, 10.89, 32.24,
4.28, 3.45, 0.73, 0.03, 11.78, 1.79, 54.69, 6.59, 100, 12.11,
94.52, 10.71, 46.24, 5.62, 11.28, 1.61, 1.15, 0.07, 9.69, 1.27,
48.48, 6.18, 95.83, 12.58, 100, 13.93, 61.25, 7.8, 21.11, 3.05,
3.82, 0.54, 0.39, 0.04, 6.76, 1.28, 37.15, 5.72, 83.53, 10.94,
100, 13.18, 71.43, 8.81, 32.86, 3.59, 8.47, 1.04, 1.51, 0.15,
0.08, 0.04, 6.2, 0.58, 28.92, 4.05, 73.89, 9.04, 100, 12.47,
76.73, 10.44, 42.32, 4.33, 14.65, 2.3, 3.12, 0.19, 0.43, 0.04,
7.15, 1.08, 42.78, 5.17, 88.11, 11.24, 100, 12.9, 60.31, 7.4,
18.46, 2.12, 2.33, 0.14, 8.03, 1, 34.79, 4.65, 82.87, 10.87,
100, 12.67, 78.03, 8.79, 30.07, 3.8, 6.64, 1.08, 0.65, 0.04,
0.04, 4.72, 0.7, 28.35, 3.4, 73.64, 10.58, 100, 13.04, 80.44,
10.89, 43.07, 5.31, 11.68, 1.56, 2.5, 0.43, 0.16, 4.15, 0.59, 
23.58, 3.22, 68.29, 9.27, 100, 13.42, 94.79, 12.32, 55.93, 7.49,
21.34, 3.22, 4.06, 0.76, 0.8, 0.13, 2.95, 0.45, 20.2, 2.68, 62.33,
7.15, 97.36, 12.82, 100, 13.99, 68.54, 9.03, 31.68, 3.49, 10.32,
1.47, 2.06, 0.09, 0.13, 0.09), DistributionLabel = c("Br1", "Br1",
"Br1", "Br1", "Br1", "Br2", "Br2", "Br2", "Br2", "Br2", "Br2",
"Br2", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3",
"Br3", "Br3", "Br4", "Br4", "Br4", "Br4", "Br4", "Br4", "Br4",
"Br4", "Br4", "Br4", "Br4", "Br5", "Br5", "Br5", "Br5", "Br5",
"Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br6",
"Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6",
"Br6", "Br6", "Br6", "Br6", "Br6", "Br7", "Br7", "Br7", "Br7",
"Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7",
"Br7", "Br7", "Br7", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8",
"Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8", 
"Br8", "Br8", "Br8", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9",
"Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9",
"Br9", "Br9", "Br9", "Br9", "Br9", "Br10", "Br10", "Br10", "Br10",
"Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Br10",
"Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Br10",
"Br10", "Br10", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl2", "Cl2",
"Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl3", "Cl3", "Cl3",
"Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl4", "Cl4", "Cl4",
"Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl5",
"Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5",
"Cl5", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6",
"Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl7", "Cl7", "Cl7", "Cl7",
"Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7",
"Cl7", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", 
"Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl9", "Cl9",
"Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9",
"Cl9", "Cl9", "Cl9", "Cl9", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10",
"Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10",
"Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Br1Cl1", "Br1Cl1", "Br1Cl1",
"Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl2", "Br1Cl2",
"Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2",
"Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3",
"Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4",
"Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4",
"Br1Cl4", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5",
"Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5",
"Br1Cl5", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", 
"Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2",
"Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2",
"Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3",
"Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl4",
"Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4",
"Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl5",
"Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5",
"Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5",
"Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1",
"Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl2", "Br3Cl2", "Br3Cl2",
"Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2",
"Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", 
"Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3",
"Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4",
"Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4",
"Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl5", "Br3Cl5", "Br3Cl5",
"Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5",
"Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5",
"Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1",
"Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl2",
"Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2",
"Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl3",
"Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3",
"Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3",
"Br4Cl3", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4",
"Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", 
"Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl5", "Br4Cl5",
"Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5",
"Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5",
"Br4Cl5", "Br4Cl5", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1",
"Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1",
"Br5Cl1", "Br5Cl1", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2",
"Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2",
"Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl3", "Br5Cl3",
"Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3",
"Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3",
"Br5Cl3", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4",
"Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", 
"Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl5", "Br5Cl5",
"Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5",
"Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5",
"Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5")), class = "data.frame", row.names = c(NA,
-606L))
    
    # Read in theoretical distributions
    tDistAll <- read.csv('Theoretical Distributions.csv', stringsAsFactors = FALSE)
    tDistAll$mz <- NA
    tDistList <- unique(tDistAll$DistributionLabel)
    bestMatchStore <- data.frame(NULL)
    
    # Iterate through experimental distributions to make an identificaiton of each
    # one. bestMatchStore holds the identification information for output.
    for (i in 1:length(distList)) {
      
      xDist <- distList[[i]] 
      resultStore <- SimilarityScoreMatch(xDist, tDistList, tDistAll)
      
      # Deterime best match by highest similarity score, and calculate the noise
      # of the experimental distribution.
      bestMatch <- resultStore[resultStore$SimilarityScore == max(resultStore$SimilarityScore), ]
      bestMatch$Noise <- with(bestMatch, round((abs(percent - percentIntensity) / percent) * 100))
      
      # Set up best match data frame.
      bestMatch$ID <- i
      bestMatch$Spectrum <- spectrumFilename
      bestMatchStore <- rbind(bestMatchStore, bestMatch)
      
    }
    
  }
  
  if(nrow(xFilter) > 0) {
    
    names(bestMatchStore) <- c('Expt_mz', 'Expt_Intenstiy', 'Expt_PercentIntensity', 'HalogenationFilter', 
                               'Theory_Intensity', 'Theory_PercentIntensity', 'DistributionLabel',
                               'SimilarityScore', 'Noise', 'ID', 'Spectrum')
    
    bestMatchStore <- bestMatchStore[bestMatchStore$HalogenationFilter == TRUE &
                                       bestMatchStore$Theory_PercentIntensity != 0, ]
  
    return(bestMatchStore)
    
  } else {
    
    return(NA)
    
  }
  
}


ProcessTableReports <- function(spectraDirectory, inputTableName, sampleName, outputTableName) {
  spectra <- read.csv(inputTableName, stringsAsFactors = FALSE)
  spectra <- spectra$SpectrumFilename[spectra$Halogenated == TRUE]
  x <- do.call('rbind', lapply(paste(spectraDirectory, '/', spectra, '.csv', sep = ''), Report))
  x$Sample <- sampleName
  tmp <- strsplit(x$Spectrum, split = "/", fixed = TRUE)
  tmp1 <- sapply(tmp, "[[", 3)
  x$SpectrumFilename <- sub(".csv", "", tmp1, fixed = TRUE)
  write.csv(x, file = outputTableName, row.names = FALSE)
}

## ----Apply functions to make additional rule tables----------------------
cat("* Beginning Step 6...\n")
tmp <- dir("Step 4 MSP to CSV") # "Sample 1 Spectra" "Sample 2 Spectra"
for (i in 1:length(tmp)) {
  sn <- sub(" Spectra", "", x = tmp[i], fixed = TRUE)
  ProcessTableReports(spectraDirectory = paste("Step 4 MSP to CSV/", tmp[i], sep = ""),
                      inputTableName = paste("Step 5 Halogenation Filter/Step 5 ", sn, ".csv", sep = ""),
                      sampleName = sn,
                      outputTableName = paste("Step 6 Additional Rules/Step 6 ", sn, ".csv", sep = ""))
}

## ----Function to apply additional rules----------------------------------
SubsetTables <- function(inputFile, chromaTOF_file, sampleName, outputFile) {
  
  # Read back in a apply additional filters.
  x <- read.csv(inputFile, stringsAsFactors = FALSE) 
  
  # Apply additonal filters.
  x <- x[x$Theory_PercentIntensity != 0, ]
  x <- x[x$DistributionLabel != 'Cl1', ]
  x <- x[x$SimilarityScore > 0.85, ]
  
  # Select only the spectra with > 1 halogenated fragment per spectrum.
  y <- setNames(aggregate(Expt_mz ~ Spectrum, data = x, FUN = length), c('Spectrum', 'Count'))
  countFilter <- y$Spectrum[y$Count > 1]
  x <- x[x$Spectrum %in% countFilter, ]
  
  # Merge back ChromaTOF table information.
  
  # Assemble ChromaTOF table.
  chromaTOF_table <- read.csv(chromaTOF_file, stringsAsFactors = FALSE)
  chromaTOF_table$SpectrumFilename <- paste("Spectrum", chromaTOF_table$PeakNumber)
  chromaTOF_table$Sample <- sampleName
  
  # Write table of only spectra passing the additional filters.
  tmp <- x[, c('Sample', 'SpectrumFilename')]
  tmp <- tmp[!duplicated(tmp), ]
  result2 <- merge(tmp, chromaTOF_table, by = c('Sample', 'SpectrumFilename'))
  
  cat('Step 7: Processing additional rules', inputFile, '...\n')
  
  write.csv(result2, file = outputFile, row.names = FALSE)
}

## ----Apply additional rules----------------------------------------------
cat("* Beginning Step 7...\n")
tmp <- dir("Step 6 Additional Rules")
for (i in 1:length(tmp)) {
  sn1 <- sub("Step 6 ", "", x = tmp[i], fixed = TRUE)
  sn2 <- sub(".csv", "", x = sn1, fixed = TRUE)
  SubsetTables(inputFile = paste("Step 6 Additional Rules/Step 6 ", sn2, ".csv", sep = ""),
                      chromaTOF_file = paste("Step 2 Input Tables/", sn2, ".csv", sep = ""),
                      sampleName = sn2,
                      outputFile = paste("Step 7 Final Results/Step 7 ", sn2, ".csv", sep = ""))
}
cat("* Processing complete.")
# To extract code:
# library(knitr)
# purl("Identify_Halogenation_Script.Rmd") 

