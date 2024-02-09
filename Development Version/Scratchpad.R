# df, t_distributions
# write.csv(df, file = 'df.csv', row.names = FALSE)
# write.csv(t_distributions, file = 'theory.csv', row.names = FALSE)
# source('Scratchpad.R', echo = FALSE)
options(width = 120)
rm(list = ls())

df <- read.csv('df.csv', stringsAsFactors = FALSE)
t_distributions <- read.csv('theory.csv', stringsAsFactors = FALSE) # theoretical distributions

# Starting after Filter 1

# Iterate through the passing filter 1 results. For each, determine the m/z range for selection and extract the experimental distribution. Convert intensities to the percent intensity within the distribution. This conversion also changes intensity values that are NA to percent_intensity values that are 0.

filter_1 <- df[df$filter_1 == TRUE,] # passing rows of filter 1
distribution_names <- unique(t_distributions$distribution_name) # theoretical distribution names
filter_2 <- data.frame(NULL) # to collect filter 2 results
monoisotopic_peaks <- data.frame(NULL) # to collect monoisotopic peaks

cat('Applying filter 2 to spectrum (count, original ID#):\n')
m <- nrow(filter_1)

for(i in 1:nrow(filter_1)) {

    n <- filter_1$spectrum_number[i]
    cat(sprintf('\r'), i, " of ", m, ", ", n, sep = '')

    mz_range <- data.frame(spectrum_number = n, nominal_mz = seq(from = filter_1$nominal_mz[i] - 10, to = filter_1$nominal_mz[i] + 13, by = 1)) # set up mz range
    e_distribution <- merge(mz_range, df, by = c('spectrum_number', 'nominal_mz'), all.x = TRUE, all.y = FALSE) # fill in missing rows in the experimental distribution
    e_distribution$intensity[is.na(e_distribution$intensity)] <- 0

    # Calculate a similarity score for the experimental distribution against every theoretical distribution. Both have 24 nominal m/z values with 10 below the maximum intensity peak and 13 above.

    match_results <- data.frame(NULL)

    for(j in 1:length(distribution_names)) {

      t_distribution <- t_distributions[t_distributions$distribution_name == distribution_names[j], ] # theoretical distribution
      match_df <- cbind(e_distribution, t_distribution) # bound distributions
      match_df <- match_df[!is.na(match_df$theoretical_mz), ] # exclude m/z range outside of theoretical distribution
      match_df$percent_intensity <- with(match_df, intensity / max(intensity) * 100) # scale the experimental distribution intensities
      u <- match_df$percent_intensity # experimental intensity vector
      v <- match_df$theoretical_percent_intensity # theoretical intensity vector
      similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
      match_df$similarity_score <- similarity_score # match results for a single theoretical distribution
      match_results <- rbind(match_results, match_df) # collect all similarity score results

    }

    best_match <- match_results[match_results$similarity_score == max(match_results$similarity_score) & match_results$percent_intensity == 100, ] # select the best match
    filter_2 <- rbind(filter_2, best_match) # collect the best matches in the sample
    monoisotopic_peak <- t_distribution[1, ] # first m/z is theoretically the monoisotopic peak
    monoisotopic_peaks <- rbind(monoisotopic_peaks, monoisotopic_peak)
}

filter_2 <- filter_2[c('spectrum_number', 'nominal_mz', 'percent_intensity', 'distribution_name', 'similarity_score')]
monoisotopic_peaks <- monoisotopic_peaks[c('spectrum_number', 'nominal_mz')]
monoisotopic_peaks$monoisotopic_peak <- TRUE 
cat(' -> Completed...\n') 

df <- merge(df, filter_2, by = c('spectrum_number', 'nominal_mz'), all.x = TRUE, all.y = FALSE)
df <- merge(df, monoisotopic_peaks, by = c('spectrum_number', 'nominal_mz'), all.x = TRUE, all.y = FALSE)
# Apply similarity score threshold. Need to assess what this should be. 
cat('Total number of halogenated spectra passing filter 2:', nrow(filter_2), '\n') # TODO apply threshold, currently the same as for filter 1
cat('Memory used by data frame with filter 1 and 2 results:', format(object.size(df), units = 'Mb'), '\n')

