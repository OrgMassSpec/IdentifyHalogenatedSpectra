# Development version of Identify Halogenated Spectra

rm(list = ls())
# library(OrgMassSpecR) # The OrgMassSpecR version must be >= 0.5-1.

# Make directory structure including each step
# cat("* Creating directory structure as needed...")
# dir.create("Step 1 Input Spectra", showWarnings = FALSE)
# dir.create("Step 2 Input Tables", showWarnings = FALSE)
# dir.create("Step 3 Split MSP Files", showWarnings = FALSE)
# dir.create("Step 4 MSP to CSV", showWarnings = FALSE)
# cat("OK\n")

### Split ChromaTOF-exported MSP file into list of individual spectra

# Terms: 'sample' is a single MSP file generated from the analysis of a sample. Process all sample MSP files in the folder, where input_sample_list is a vector of paths to the input MSP spectra files.
input_sample_list <- list.files("Step 1 Input Spectra", full.names = TRUE) 

# for (i in 1:length(input_sample_list)) { # Increment through path list
# TODO complete all tasks for single file then use rm(list()) to clear memory except for functions before going to the next file?
i <- 1 # For dev: work with single file
working_sample <- input_sample_list[i] # Select path from list
tmp_name <- strsplit(working_sample, split = '/', fixed = TRUE)[[1]][2] # Select file name from path
working_sample_name <- substr(tmp_name, 1, nchar(tmp_name) - 4) # Remove '.txt' from file name # TODO Use or delete?
working_sample_spectra <- readLines(working_sample) # Read file in to character vector 

cat('Processing ', tmp_name, '... \n', sep = '')

# Make a vector (spectrum_factor) that defines the separate spectra within the MSP file, where spectrum_factor is defined by spectrum_count, which iterates by + 1 if it encounters a newline that defines the next spectrum. 
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

# Preallocate the data frame size (N) with an estimate of the total number of m/z peaks in the file: find the number of rows in each spectrum, exclude 3 rows per spectrum for the headers and blank line, find the sum total, and multiply by 10 (the max number of allowed m/z peaks per row). Create the empty data frame (df).
N <- sum(lengths(spectrum_list) - 3) * 10 
df <- data.frame(peak_number = rep(NA, N), subnominal_mz = rep(NA, N), nominal_mz = rep(NA, N), intensity = rep(NA, N))

# Track data frame row number (row_num) and iterate through the individual spectra (spectrum_list_i) from spectrum_list. Split out peaks (m/z and intensity pairs) into list subelements (x). Remove `character(0)` elements from list that results from newlines.
n <- length(spectrum_list)
row_num <- 1 
cat('Total number of spectra:', n, '\n')
cat('Reading in spectrum:\n')

for(spectrum_list_i in 1:length(spectrum_list)) { 

  cat(sprintf('\r'), spectrum_list_i, sep = '')

  x <- strsplit(spectrum_list[[spectrum_list_i]], split = ";", fixed = TRUE)
  x <- x[lengths(x) != 0] 

  # Within an individual spectrum, iterate through the rows (i), skiping NAME and Num Peaks. Then iterate through the m/z peaks (m/z and intensity pairs) within a row (j) and make a vector of the pair (mz_int). Remove if empty and use round() to deterine the nominal m/z. Construct a list (a) of the values and add this as a row to df.
  for(i in 3:length(x)) { 
    for(j in 1:length(x[[i]])) { 
      mz_int <- strsplit(x[[i]][j], split = "[[:space:]]")[[1]] 
      mz_int <- mz_int[mz_int != ""] 
      nominal_mz <- round(as.numeric(mz_int[1]))
      a <- list(as.numeric(spectrum_list_i), as.numeric(mz_int)[1], as.numeric(nominal_mz), as.numeric(mz_int)[2]) 
      df[row_num, ] <- a
      row_num <- row_num + 1

    }
  }
}

### Check for duplicated nominal m/z values and combine

# Remove unused rows (NA). For the remaining data, sum the intensities of duplicated nominal_mz values within the spectrum. The subnominal_mz for the the m/z peak with the larger original intensity is kept, but with the original intensity replaced by the summed intensity. The smaller intensity m/z peak is deleted. Columns are then reordered.
df <- na.omit(df) 
df_aggregated <- aggregate(df['intensity'], by = df[c('peak_number', 'nominal_mz')], sum)
df_combined <- merge(df, df_aggregated, by = c('peak_number', 'nominal_mz'))
df_ordered <- df_combined[order(df_combined$peak_number, df_combined$nominal_mz, -df_combined$intensity.x), ] 
df_read <- df_ordered[!duplicated(df_ordered[c('peak_number', 'nominal_mz', 'intensity.y')]), ]
rownames(df_read) <- NULL
df_read <- df_read[, c('peak_number', 'subnominal_mz', 'nominal_mz', 'intensity.y')]
names(df_read)[names(df_read) == 'intensity.y'] <- 'intensity'

# TODO Make Intensity threshold settable
# REPORT status to separate output file and terminal. 
# REPORT number of dulicated nominal m/z peaks per spectrum
# REPORT number of m/z peaks under intensity threshold
# TODO Double check with more complex examples

### Remove noise spectra

rm(list = setdiff(ls(), "df_read"))
cat(' -> Completed...\n')
cat('Memory used to read data frame:', format(object.size(df_read), units = 'Mb'), '\n')

# Remove spectra with less than 10 ions above m/z 100. These are 'noise spectra'. 
n_removed <- 0
for(i in unique(df_read$peak_number)) {
  tmp <- df_read[df_read$peak_number == i, ]
  if(nrow(tmp[tmp$nominal_mz > 100, ]) < 10) {
    df_read <- df_read[df_read$peak_number != i, ]
    n_removed <- n_removed + 1
  }    
}
cat('Number of noise spectra removed:', n_removed, '\n')
cat('New total number of spectra:', length(unique(df_read$peak_number)), '\n')

### Make data frame for input into filters

# Add in placeholders for all nominal m/z values
# Construct df with all nominal mz values and merge. Find min and max values for each peak_number and build df using rep().  
cat('Adding in missing nominal m/z values... \n')
z <- data.frame(NULL)
for(i in unique(df_read$peak_number)) {
  tmp <- df_read[df_read$peak_number == i, ]
  x_min <- min(tmp$nominal_mz)
  x_max <- max(tmp$nominal_mz)
  nominal_values <- x_min:x_max
  y <- data.frame(peak_number = i, nominal_mz = nominal_values)
  z <- rbind(z, y)
}

df <- merge(df_read, z, by = c('peak_number', 'nominal_mz'), all = TRUE)
df$intensity[is.na(df$intensity)] <- 0
rownames(df_read) <- NULL
cat('Memory used by full nominal m/z data frame:', format(object.size(df), units = 'Mb'), '\n')

#write.csv(df_read, file = 'Output.txt')

### Halogenation filter 1
# TODO Cite paper
rm(list = setdiff(ls(), "df"))

# Normalize peak intensity to percentage of max peak and only examine m/z > 100 amu.
cat('Applying filter 1 to spectrum (original ID#):\n')
df$filter_1 <- FALSE
y <- data.frame(NULL)

for(i in unique(df$peak_number)) {
  cat(sprintf('\r'), i, sep = '')
  x <- df[df$peak_number == i, ]
  x$intensity <- with(x, intensity / max(intensity) * 100)
  x <- x[x$nominal_mz > 100, ]

  for(i in (nrow(x) - 6):5) {
    
    if(x$intensity[i] > 3) {
      
      # Calculate peak ratios.
      
      peakRatio1 <- x$intensity[i + 1] / x$intensity[i]
      peakRatio2 <- x$intensity[i + 2] / x$intensity[i]
      peakRatio3 <- x$intensity[i + 3] / x$intensity[i]
      peakRatio4 <- x$intensity[i + 4] / x$intensity[i]
      peakRatio5 <- x$intensity[i + 5] / x$intensity[i]
      peakRatio6 <- x$intensity[i + 6] / x$intensity[i]
      peakRatioa <- x$intensity[i - 1] / x$intensity[i]
      peakRatiob <- x$intensity[i - 2] / x$intensity[i]
      peakRatioc <- x$intensity[i - 4] / x$intensity[i]
      
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
        
        x$filter_1[i] <- TRUE 
      }
    } 
  } 
  x <- x[, c('peak_number', 'nominal_mz', 'filter_1')]
  y <- rbind(y, x)

}

# Dataframe df holds the results. Remove and rename columns.
df <- merge(df, y, by = c('peak_number', 'nominal_mz'))
df <- subset(df, select = -filter_1.x)
names(df)[names(df) == 'filter_1.y'] <- 'filter_1'
cat(' -> Completed...\n') 
rm(list = setdiff(ls(), "df"))
cat('Number of halogenated spectra per filter 1:', length(unique(df$peak_number[df$filter_1 == TRUE])), '\n')
cat('Memory used by data frame with filter 1 results:', format(object.size(df), units = 'Mb'), '\n')

### Halogenation filter 2
# TODO filter 2

# Is this needed? Wouldn't similarity score be poor against a second peak in the distribution anyways? Computationally it is easier to just calculate the similarity score against all filter_1 == TRUE. 

  #SelectHalogenatedDistribution <- function(xFilter, x) { 
  #
  ## xFilter is the full data frame and x is the results of filter 1
  #
  #
  #  # Remove miss-identified peaks within the same distribution. Sometimes another
  #  # peak within the same distribution is identified. Only the maximum peak 
  #  # within the distribution should be identified. The following code eliminates
  #  # any multiple identifications within the same distribution, although not
  #  # necessairly the maximum one. The maximum peak is identified again later when
  #  # plotting the distribution.
  #  
  #  xFilter$Diff <- c(diff(xFilter$mz), 100) # 100 is placeholder at the end
  #  xFilter <- xFilter[xFilter$Diff > 10, ] # exclude if within 10 mz units
  #  
  #  # Select the 5 most abundant distributions. Then return to the original order.
  #  if(nrow(xFilter) > 5) {
  #    
  #    xFilter <- head(xFilter[rev(order(xFilter$intensity)), ], n = 5) 
  #    xFilter <- xFilter[order(xFilter$mz), ]
  #    
  #  }
  #  
  #  # Select window of m/z around the experimental distributions and store in
  #  # distList.
  #  
  #  distList <- vector(mode = 'list')
  #
  #
## NEED THIS
#  # Slices -10 to + 10 around the filter_1 = true result? This then goes to the match against the theoretical distributions.
#  for (i in 1:nrow(xFilter)) {
#    distList[[i]] <- x[x$mz %in% seq(from = xFilter$mz[i] - 10, to = xFilter$mz[i] + 10, by = 1), ]  
#  }
#  return(distList)
#  
#}
#
## Isolate halogenated distributions. xFilter contains only ions identifed as
#  # halogenated.
#  xFilter <- x[x$Filter == TRUE, ]
#  
#  # Need to skip if no halogenation was detected.
#  if(nrow(xFilter) > 0) {
#    distList <- SelectHalogenatedDistribution(xFilter, x)
#    
#    # Iterate through experimental distributions.
#    for (i in 1:length(distList)) {
#      xDist <- distList[[i]]
#    }
#    
#
#
## Implement this.
#SimilarityScoreMatch <- function(xDist, tDistList, tDistAll) {
#  
#  resultStore <- data.frame(NULL)
#  
#  # Iterate through each theoretical distribution and calculate a similarity
#  # score against the experimental spectrum.
#  for (j in 1:length(tDistList)) {
#    
#    tDist <- tDistAll[tDistAll$DistributionLabel == tDistList[j], ] # current theoretical distribution
#    
#    # Align theoretical distribution to experimental distribution based on the
#    # maximum peak (maxvalue). The theoretical distribution m/z window is
#    # given the same size as the experimental m/z window.
#
#
#    # Where is it getting xDist. 
#    tDist$mz[tDist$percent == 100] <- xDist$mz[which.max(xDist$percentIntensity)]
#    maxvalue <- xDist$mz[which.max(xDist$percentIntensity)]
#    index <- which.max(tDist$mz)
#    
#    if(index > 1) {
#      startValue <- maxvalue - (index - 1)
#      endValue <- startValue + nrow(tDist) - 1
#      tDist$mz <- seq(from = startValue, to = endValue, by = 1)
#    }
#    
#    # Need special case when the max ion is the first one.
#    if(index == 1) {
#      startValue <- maxvalue
#      endValue <- maxvalue + nrow(tDist) - 1
#      tDist$mz <- seq(from = startValue, to = endValue, by = 1)
#    }
#    
#    # Merge experimental and theoretical distributions by the m/z.
#    cDist <- merge(xDist, tDist, by = 'mz', all.x = TRUE, all.y = FALSE)
#    cDist$percent[is.na(cDist$percent)] <- 0
#    cDist$DistributionLabel <- tDistList[j]
#    
#    # Calculate similarity score, add to dataframe, and store in resultStore.
#    u <- (cDist$intensity.x/max(cDist$intensity.x)) * 100
#    v <- cDist$percent
#    cDist$SimilarityScore <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
#    resultStore <- rbind(resultStore, cDist)
#    
#  }
#  
#  return(resultStore)
#  
#}
#
## Move this code to start of filter 2?
#
#  # Read in theoretical distributions
#    tDistAll <- read.csv('Theoretical Distributions.csv', stringsAsFactors = FALSE)
#    tDistAll$mz <- NA
#    tDistList <- unique(tDistAll$DistributionLabel)
#    bestMatchStore <- data.frame(NULL)
#    
#    # Iterate through experimental distributions to make an identificaiton of each
#    # one. bestMatchStore holds the identification information for output.
#    for (i in 1:length(distList)) {
#      
#      xDist <- distList[[i]] 
#      resultStore <- SimilarityScoreMatch(xDist, tDistList, tDistAll)
#      
#      # Deterime best match by highest similarity score, and calculate the noise
#      # of the experimental distribution.
#      bestMatch <- resultStore[resultStore$SimilarityScore == max(resultStore$SimilarityScore), ]
#      bestMatch$Noise <- with(bestMatch, round((abs(percent - percentIntensity) / percent) * 100))
#      
#      # Set up best match data frame.
#      bestMatch$ID <- i
#      bestMatch$Spectrum <- spectrumFilename
#      bestMatchStore <- rbind(bestMatchStore, bestMatch)
#      
#    }
#    
#  }

# Print similarity score in results data frame along with TRUE FALSE result.