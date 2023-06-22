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

# Terms: 'sample' is a single MSP file generated from the analysis of a sample, The script will process all sample MSP files in the folder.  input_sample_list is a vector of paths to the input MSP spectra files.
input_sample_list <- list.files("Step 1 Input Spectra", full.names = TRUE) 

# for (i in 1:length(input_sample_list)) { # Increment through path list
# TODO complete all tasks for single file then use rm(list()) to clear memory except for functions before going to the next file?
i <- 1 # For dev: work with single file
working_sample <- input_sample_list[i] # Select path from list
tmp_name <- strsplit(working_sample, split = '/', fixed = TRUE)[[1]][2] # Select file name from path
working_sample_name <- substr(tmp_name, 1, nchar(tmp_name) - 4) # Remove '.txt' from file name
working_sample_spectra <- readLines(working_sample) # Read file in to character vector 

# Make a vector (spectrum_factor) that defines the separate spectra within the file. spectrum_factor is defined by spectrum_count, which iterates by + 1 if it encounters a newline that defines the next spectrum. 
spectrum_factor <- vector(mode = "integer", length = length(working_sample_spectra))
spectrum_count <- 1
for(i in 1: length(working_sample_spectra)) {
  if(working_sample_spectra[i] != "") { 
    spectrum_factor[i] <- spectrum_count 
  } else {
    spectrum_factor[i] <- spectrum_count
    spectrum_count <- spectrum_count + 1 
  }
}

# Using spectrum_factor, split the spectra into list elements as defined by spectrum_factor and convert to a data frame.
spectrum_list <- split(working_sample_spectra, f = spectrum_factor)

# TODO If end of file contains more than one newline, ignore
# REPORT number of spectra in file

# Preallocate the data frame size (N) with an estimate of the total number of m/z peaks in the file: find the number of rows in each spectrum, exclude 3 rows per spectrum for the headers and blank line, find the sum total, and multiply by 10 (the max number of allowed m/s peaks per row). Create the empty data frame (df).
N <- sum(lengths(spectrum_list) - 3) * 10 
df <- data.frame(sample = rep("", N), peak_number = rep(NA, N), 
  subnominal_mz = rep(NA, N), nominal_mz = rep(NA, N), intensity = rep(NA, N))

# Track data frame row number (row_num) and iterate through the individual spectra (spectrum_list_i) from spectrum_list. Split out peaks (m/z and intensity pairs) into list subelements (x). Remove `character(0)` elements from list that results from newlines.

n <- length(spectrum_list)

row_num <- 1 

for(spectrum_list_i in 1:length(spectrum_list)) { 

  cat('Processing spectrum', spectrum_list_i, 'of', n, '\n')

  x <- strsplit(spectrum_list[[spectrum_list_i]], split = ";", fixed = TRUE)
  x <- x[lengths(x) != 0] 

  # Within an individual spectrum, iterate through the rows (i), skiping NAME and Num Peaks. Then iterate through the m/z peaks (m/z and intensity pairs) within a row (j) and make a vector of the pair (mz_int). Remove if empty and use round() to deterine the nominal m/z. Construct a list (a) of the values and add this as a row to df.
  for(i in 3:length(x)) { 
    for(j in 1:length(x[[i]])) { 
      mz_int <- strsplit(x[[i]][j], split = "[[:space:]]")[[1]] 
      mz_int <- mz_int[mz_int != ""] 
      nominal_mz <- round(as.numeric(mz_int[1]))

      a <- list(working_sample_name, as.numeric(spectrum_list_i), as.numeric(mz_int)[1], 
        as.numeric(nominal_mz), as.numeric(mz_int)[2]) 

      df[row_num, ] <- a
      row_num <- row_num + 1

    }
  }
}

# Remove unused rows and rows with and intensity below the threshold. For the remaining data, sum the intensities of duplicated nominal_mz values within the spectrum. The subnominal_mz for the the m/z peak with the larger original intensity is kept, but with the original intensity replaced by the summed intensity. The smaller intensity m/z peak is deleted. Columns are then reordered.
df <- na.omit(df) 
df <- df[df$intensity > 50, ]    
df_aggregated <- aggregate(df['intensity'], by = df[c('peak_number', 'nominal_mz')], sum)
df_combined <- merge(df, df_aggregated, by = c('peak_number', 'nominal_mz'))
df_ordered <- df_combined[order(df_combined$peak_number, df_combined$nominal_mz, -df_combined$intensity.x), ] 
df_complete <- df_ordered[!duplicated(df_ordered[c('peak_number', 'nominal_mz', 'intensity.y')]), ]
rownames(df_complete) <- NULL
df_complete <- df_complete[, c('sample', 'peak_number', 'subnominal_mz', 'nominal_mz', 'intensity.y')]
names(df_complete)[names(df_complete) == 'intensity.y'] <- 'intensity'

# TODO Make Intensity threshold settable
# REPORT status to separate output file and terminal. 
# REPORT number of dulicated nominal m/z peaks per spectrum
# REPORT number of m/z peaks under intensity threshold

# TODO Double check with more complex examples

write.csv(df_complete, file = 'Output.txt')
