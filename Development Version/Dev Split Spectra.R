# Dev version for MSP spectrum

### Setup ###
rm(list = ls())
library(OrgMassSpecR) # The OrgMassSpecR version must be >= 0.5-1.

# Make directory structure including each step
cat("* Creating directory structure as needed...")
dir.create("Step 1 Input Spectra", showWarnings = FALSE)
dir.create("Step 2 Input Tables", showWarnings = FALSE)
dir.create("Step 3 Split MSP Files", showWarnings = FALSE)
dir.create("Step 4 MSP to CSV", showWarnings = FALSE)
cat("OK\n")

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
  browser()
  # TODO If end of file contains more than one newline, ignore.

  # Split vector `a` containing the spectra by the defined spectrumFactor.
  spectrumList <- split(a, f = spectrumFactor)
  
  # Make individual files, where the filename specifies the original order in 
  # the input file. Removes blank lines and writeLines adds a single blank line
  # to the bottom of the file.

  # TODO Add conversion to CSV here. Exclude first to lines.

  for(j in 1:length(spectrumList)) {
    individualSpectrum <- spectrumList[[j]]
    individualSpectrum <- individualSpectrum[individualSpectrum != ""]
    fullFilename <- paste(outputFolderName, "/Spectrum ", j, ".txt", sep = "")
    writeLines(individualSpectrum, fullFilename)
    # Reporting
    cat('Step 3: Splitting ', fullFilename, '...\n', sep = "")
  }
}

# Call function SplitMSP to spilt into individual MSP files. 
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