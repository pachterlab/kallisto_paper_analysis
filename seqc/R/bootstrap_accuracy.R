library("cowplot")
library("sleuth")
library("parallel")

# base_dir <- "../results/paired"
# sample_id <- dir(base_dir)
#
# kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, "kallisto"))
#
# kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, "kallisto_elen"))
#
# kal <- read_kallisto(kal_dirs[1])
kallisto_bootstrap_path  <- file.path('../results/A/BGI', 1, 'kallisto_bootstrap', "abundance.h5")
kal <- read_kallisto_h5(kallisto_bootstrap_path)

options(mc.cores = 40)

bootstrap_summaries <- function(kal) {
  num_bs <- length(kal$bootstrap)

  all_summary <- sleuth:::summarize_bootstrap(kal, 'est_counts')
  mclapply(2:num_bs,
    function(i) {
      message('iteration: ', i)
      cur_kal <- kal
      cur_kal$bootstrap <- cur_kal$bootstrap[1:i]
      cur_summary <- sleuth:::summarize_bootstrap(cur_kal, 'est_counts')

      cur_summary
    })
}

bootstrap_accuracy <- function(summaries) {
  final <- summaries[[length(summaries)]]
  all_accuracy <- lapply(seq_along(summaries),
    function(i) {
      x <- summaries[[i]]
      data.frame(i = i + 1,
        rmse = rmse(x, final),
        mrd = median(
          abs(mamabear:::relative_difference(
            x$bs_sd_est_counts,
            final$bs_sd_est_counts,
            na_zeroes = TRUE)),
          na.rm = TRUE),
        re = median(
          abs(relative_error(x$bs_sd_est_counts, final$bs_sd_est_counts)),
          na.rm = TRUE)
          )
    })
  bind_rows(all_accuracy)
}

rmse <- function(estimate, truth) {
  sqrt(mean((estimate$bs_sd_est_counts - truth$bs_sd_est_counts)^2))
}

relative_error <- function(estimate, truth) {
  which_filter <- truth > 0

  (estimate[which_filter] - truth[which_filter]) / truth[which_filter]
}

all_summaries <- bootstrap_summaries(kal)



#####
all_errors <- bootstrap_accuracy(all_summaries)
saveRDS(all_summaries, file = "all_summaries.rds")
saveRDS(all_errors, file = "all_errors.rds")

ggplot(all_errors, aes(i, error)) +
  geom_point()

ggplot(all_errors, aes(i, mrd)) +
  geom_point()
