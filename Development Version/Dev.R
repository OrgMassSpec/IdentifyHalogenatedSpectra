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
working_sample_spectra <- readLines(working_sample) # Read file in to character vector 

# Make a vector that defines the separate spectra from the file
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

# Using spectrum_factor, split the spectra into list elements and convert to a data frame
spectrum_list <- split(working_sample_spectra, f = spectrum_factor) # Split into list elements defined by spectrum_factor

# REPORT number of spectra in file

N <- sum(lengths(spectrum_list) - 3) * 10 # Preallocate the data frame size with an estimate of the total number of m/z peaks in the file: find the number of rows in each spectrum, exclude 3 rows per spectrum for the headers and blank line, find the sum total, and multiply by 10 (the max number of allowed m/s peaks per row)
results <- data.frame(sample = rep("", N), peak_number = rep(NA, N), 
  subnominal_mz = rep(NA, N), nominal_mz = rep(NA, N), intensity = rep(NA, N)) # Empty data frame

row_num <- 1 # Track data frame row number
for(ii in 1:length(spectrum_list)) { # Get individual spectrum
# ii <- 1
  x <- strsplit(spectrum_list[[ii]], split = ";", fixed = TRUE) # Split out peaks (m/z and intensity pairs) into list subelements

  x <- x[lengths(x) != 0] # Remove `character(0)` elements from list that results from newlines
  
  for(i in 3:length(x)) { # Iterate through rows, skiping NAME and Num Peaks
    for(j in 1:length(x[[i]])) { # Iterate through m/z peaks within a row
      y <- strsplit(x[[i]][j], split = "[[:space:]]")[[1]] # Make a vector of the m/z intensity pair
      z <- y[y != ""] # Remove if empty
      nominal_mz <- round(as.numeric(z[1])) # Make nominal mass

      a <- list(working_sample_name, as.numeric(ii), as.numeric(z)[1], 
        as.numeric(nominal_mz), as.numeric(z)[2]) 

      # TODO Try data.table:rbindlist instead of rbind
      # TODO Implement intensity threshold. Apply at end?
      # TODO Report status to separate output file and terminal. 

      #results <- rbind(results, z, make.row.names = FALSE) # Slow?
      results[row_num, ] <- a
      row_num <- row_num + 1

      results <- na.omit(results) # Remove unused rows

    }
  }
}

results <- aggregate(results['intensity'], by = results[c('peak_number', 'nominal_mz')], sum) # Sum the intensities of duplicated nominal_mz values within a spectrum 
# duplicated(results[, c("peak_number", "nominal_mz")]) # To manually check for duplicates

# TODO resort results

# REPORT number of dulicated nominal m/z peaks per spectrum

write.csv(results, file = 'Output.txt')
