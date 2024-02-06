# df, t_distributions
# write.csv(df, file = 'df.csv', row.names = FALSE)
# write.csv(t_distributions, file = 'theory.csv', row.names = FALSE)
options(width = 200)
rm(list = ls())

df <- read.csv('df.csv', stringsAsFactors = FALSE)
t_distributions <- read.csv('theory.csv', stringsAsFactors = FALSE)

# Starting after Filter 1

# Iterate through the passing filter 1 results. For each, determine the m/z range for selection and extract the experimental distribution. Convert intensities to the percent intensity within the distribution. This conversion also changes intensity values that are NA to percent_intensity values that are 0.

filter_1 <- df[df$filter_1 == TRUE,] # passing rows of filter 1
distribution_names <- unique(t_distributions$distribution_name) # theoretical distribution names
filter_2 <- data.frame(NULL) # filter 2 results

for(i in 1:nrow(filter_1)) {
    mz_range <- data.frame(nominal_mz = seq(from = filter_1$nominal_mz[i] - 10, to = filter_1$nominal_mz[i] + 13, by = 1)) # set up mz range
    e_distribution <- df[df$spectrum_number == filter_1$spectrum_number[i] & df$nominal_mz %in% mz_range$nominal_mz, ] # extract the experimental distribution
    e_distribution$percent_intensity <- with(e_distribution, intensity / max(intensity) * 100) # scale the distribution

    # Calculate a similarity score for the experimental distribution against every theoretical distribution. The experimental distribution tmp_experimental and the theoretical distribution tmp_theoretical both have 24 nominal m/z values with 10 below the maximun peak and 13 above the maximum peak. Nominal m/z values without experimental intensities are remove before the similarity score is calculated.

    match_results <- data.frame(NULL)

    for(j in 1:length(distribution_names)) {
      t_distribution <- t_distributions[t_distributions$distribution_name == distribution_names[j], ] # theoretical distribution
      match_df <- cbind(e_distribution, t_distribution) # bound distributions
      match_df <- match_df[!is.na(match_df$theoretical_mz), ] # exclude m/z range outside of theoretical distribution
      u <- match_df$percent_intensity # experimental intensity vector
      v <- match_df$theoretical_percent_intensity # theoretical intensity vector
      similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
      match_df$similarity_score <- similarity_score # match results for a single theoretical distribution
      match_results <- rbind(match_results, match_df) # combined results
    }

    best_match <- match_results[match_results$similarity_score == max(match_results$similarity_score) & match_results$percent_intensity == 100, ] 
  }

# TODO merge results data frame (a) with df and apply threshold for passing filter 2


