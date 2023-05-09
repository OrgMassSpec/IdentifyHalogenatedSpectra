# Dev version for MSP spectrum

rm(list = ls())
# library(OrgMassSpecR) # The OrgMassSpecR version must be >= 0.5-1.

# Make directory structure including each step
# cat("* Creating directory structure as needed...")
# dir.create("Step 1 Input Spectra", showWarnings = FALSE)
# dir.create("Step 2 Input Tables", showWarnings = FALSE)
# dir.create("Step 3 Split MSP Files", showWarnings = FALSE)
# dir.create("Step 4 MSP to CSV", showWarnings = FALSE)
# cat("OK\n")

# Split ChromaTOF-exported MSP into list of individual MSP spectra

# Read input MSP file into character vector of one line per element
input_sample_list <- list.files("Step 1 Input Spectra", full.names = TRUE) # input_sample_list is a vector of paths to the input MSP spectra files

# for (i in 1:length(input_sample_list)) { # Increment through path list
i <- 1 # For dev: work with single file
working_sample <- input_sample_list[i] # Select path from list
tmp_name <- strsplit(working_sample, split = '/', fixed = TRUE)[[1]][2] # Select file name from path
working_sample_name <- substr(tmp_name, 1, nchar(tmp_name) - 4) # Remove '.txt' from file name
working_sample_spectra <- readLines(working_sample) # Read in to chr vector of MSP file.

# Make vector defining the separate spectra
spectrum_factor <- vector(mode = "integer", length = length(working_sample_spectra)) # Empty vector to specify each spectrum
spectrum_count <- 1
for(i in 1: length(working_sample_spectra)) {
  if(working_sample_spectra[i] != "") { 
    spectrum_factor[i] <- spectrum_count # spectrum_factor same for lines within a spectrm
  } else {
    spectrum_factor[i] <- spectrum_count
    spectrum_count <- spectrum_count + 1 # spectrum_count iterates by +1 if it encounters a newline (defining the next spectrum)
  }
}

# TODO If end of file contains more than one newline, ignore

# Split spectra into list elements and convert to data frame
spectrum_list <- split(working_sample_spectra, f = spectrum_factor) # Split into list elements defined by spectrum_factor

  # > str(spectrum_list)
  # List of 3
  # $ 1: chr [1:10] "NAME: 1,1'-Bicyclohexyl" "Num Peaks: 64" "50.03646 3.46;51.04831 11.13;52.05490 6.61;...

results <- data.frame(subnominal_mz = NULL, intensity = NULL, sample = NULL, peak_number = NULL, nominal_mz = NULL) # Empty data frame
for(ii in 1:length(spectrum_list)) { # Get individual spectrum
  x <- strsplit(spectrum_list[[ii]], split = ";", fixed = TRUE) # Split out peaks (m/z and intensity pairs) into list subelements

  # > str(x)
  # List of 10
  #  $ : chr "NAME: 1,1'-Bicyclohexyl"
  #  $ : chr "Num Peaks: 64"
  #  $ : chr [1:10] "50.03646 3.46" "51.04831 11.13" "52.05490 6.61" "53.06468 49.07" ...
  #  $ : chr [1:10] "63.05270 1.24" "65.06774 14.55" "65.28902 1.35" "66.07368 11.33" ...

  x <- x[lengths(x) != 0] # Remove `character(0)` elements from list that results from newlines
  for(i in 3:length(x)) { # Iterate through rows containing peaks, skiping NAME and Num Peaks
    for(j in 1:length(x[[i]])) { # Iterate through peaks within a row
      y <- strsplit(x[[i]][j], split = "[[:space:]]")[[1]] # Vector of m/z intensity pair
      z <- y[y != ""] # Remove if empty
      nominal_mz <- round(as.numeric(z[1])) # Make nominal mass
      z <- c(z, working_sample_name, ii, nominal_mz) # Accurate m/z, intensity, sample name, spectrum number, nominal m/z

      # TODO NEXT Bin and add subnominal intensity. Look for repeating nominal m/z to identify? Add subnominal intensities.
      # TODO NEXT Add missing integer nominal_mz values to data frame.
      # TODO Try data.table:rbindlist instead of rbind
      # TODO Implement intensity threshold. Apply at end?

      results <- rbind(results, z)
      
      # TODO NEXT Rename headers

    }
  }
}

# Has doublet:
# 24     79.07735   39.75    PCB Mix v2    1    79
# 25     79.29673    5.60    PCB Mix v2    1    79


write.csv(results, file = 'Output.txt')
# for(i in 3:length(x)) {
#   for(j in 1:length(x[[i]])) {
#     y <- strsplit(x[[i]][j], split = "[[:space:]]")
#     z <- y[y != ""]
#     results <- rbind(results, z)
#   }
# }

# RemoveWhite <- function(x) {
#     y <- strsplit(x, split = "[[:space:]]")[[1]]
#     z <- y[y != ""]
#     return(as.numeric(z))
#   }
# results.list <- lapply(x, RemoveWhite)
 
# 
# y <- as.data.frame(do.call("rbind", x))
# values <- NULL
#   for (i in 1:ncol(y)) {
#     values <- c(values, y[, i])
#   }
# tall.format <- as.vector(values, mode = "character")
# 
# # lapply(y, FUN = )
# 
# RemoveWhite <- function(x) {
#     y <- strsplit(x, split = "[[:space:]]")[[1]]
#     z <- y[y != ""]
#     return(as.numeric(z))
#   }
#   results.list <- lapply(x, RemoveWhite)

# Make individual files, where the filename specifies the original order in 
# the input file. Removes blank lines and writeLines adds a single blank line
# to the bottom of the file.

# TODO Add conversion to CSV here. Exclude first two lines. Check last line for null.
# TODO Add nominal mass placeholders

# for(j in 1:length(spectrumList)) {
#   individualSpectrum <- spectrumList[[j]]
#   individualSpectrum <- individualSpectrum[individualSpectrum != ""]
#   fullFilename <- paste(outputFolderName, "/Spectrum ", j, ".txt", sep = "")
#   writeLines(individualSpectrum, fullFilename)
#   # Reporting
#   cat('Step 3: Splitting ', fullFilename, '...\n', sep = "")
# }


# Call function SplitMSP to spilt into individual MSP files. 
# cat("* Beginning Step 3...\n")
# tmp <- list.files("Step 1 Input Spectra", full.names = TRUE)
# 
# for (i in 1:length(tmp)) {
#   # Capture MSP file name to use as folder name for split files.
#   f1 <- strsplit(tmp[i], split = c("/"), fixed = TRUE)[[1]][2]
#   folder <- substr(f1, 1, nchar(f1) - 4) # remove .txt
#   dir.create(paste("Step 3 Split MSP Files/", folder, " Spectra", sep = ""))
#   SplitMSP(originFileName = tmp[i], 
#            outputFolderName = paste("Step 3 Split MSP Files/", folder, " Spectra", sep = ""))
# }

## ----Function to convert MSP to CSV--------------------------------------
# rm(list = ls())
# 
# ConvertMSP <- function(file) {
#   x <- ReadMspFile(file, skip = 2, comment.char = "", remove.placeholders = FALSE)
#   # If reporting is needed:
#   # cat('Making ', file, '\n', sep = "")
#   f1 <- strsplit(file, split = c("/"), fixed = TRUE)[[1]][3] # filename
#   f2 <- strsplit(file, split = c("/"), fixed = TRUE)[[1]][2] # sample directory
#   fileName <- substr(f1, 1, nchar(f1) - 4) # remove .txt
#   cat('Step 4: Converting ', fileName, '...\n', sep = "")
#   write.csv(x, file = paste("Step 4 MSP to CSV/", f2, "/", fileName, ".csv", sep = ''), row.names = FALSE)
# }
# 
# ## ----Convert MSP to CSV--------------------------------------------------
# cat("* Beginning Step 4...\n")
# tmp <- list.dirs("Step 3 Split MSP Files", recursive = FALSE)
# 
# for (i in 1:length(tmp)) {
#   inputSpectraDirectory <- tmp[i]
#   spectra <- dir(inputSpectraDirectory)
#   f <- strsplit(tmp[i], split = c("/"), fixed = TRUE)[[1]][2] # sample directory
#   dir.create(paste("Step 4 MSP to CSV/", f, sep = ""))
#   x <- do.call('rbind', lapply(paste(inputSpectraDirectory, '/', spectra, sep = ''), ConvertMSP))
# }