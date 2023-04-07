# Dev version for MSP spectrum

### Setup ###
rm(list = ls())
# library(OrgMassSpecR) # The OrgMassSpecR version must be >= 0.5-1.

# Make directory structure including each step
# cat("* Creating directory structure as needed...")
# dir.create("Step 1 Input Spectra", showWarnings = FALSE)
# dir.create("Step 2 Input Tables", showWarnings = FALSE)
# dir.create("Step 3 Split MSP Files", showWarnings = FALSE)
# dir.create("Step 4 MSP to CSV", showWarnings = FALSE)
# cat("OK\n")

# Split MSP into list of individual MSP

input_sample_list <- list.files("Step 1 Input Spectra", full.names = TRUE)
# for (i in 1:length(input_sample_list)) { # add closing bracket later
i <- 1
working_sample <- input_sample_list[i] 
working_sample_spectra <- readLines(working_sample)

# Make a factor specifying each spectrum
spectrum_factor <- vector(mode = "integer", length = length(working_sample_spectra)) 
spectrum_count <- 1
# Variable spectrum_count is the same for each line of a spectrum, then iterates by +1 if it encounters a newline (defining the next spectrum).
for(i in 1: length(working_sample_spectra)) {
  if(working_sample_spectra[i] != "") {
    spectrum_factor[i] <- spectrum_count
  } else {
    spectrum_factor[i] <- spectrum_count
    spectrum_count <- spectrum_count + 1
  }
}

# TODO If end of file contains more than one newline, ignore.

# Split vector `a` containing the spectra by the defined spectrumFactor.
spectrum_list <- split(working_sample_spectra, f = spectrum_factor)

# TODO NEXT spectrum_list[[1]] below will need to iterate through spectrum list
x <- strsplit(spectrum_list[[1]], split = ";", fixed = TRUE)
# Remove `character(0)` element from list that results from the newline
x <- x[lengths(x) != 0]

results <- data.frame(mz = NULL, intensity = NULL, sample = NULL)

for(i in 3:length(x)) {
  for(j in 1:length(x[[i]])) {
    y <- strsplit(x[[i]][j], split = "[[:space:]]")[[1]]
    z <- y[y != ""]
    z <- c(z, working_sample)

# TODO Try data.table:rbindlist instead of rbind

    results <- rbind(results, z)
  }
}

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