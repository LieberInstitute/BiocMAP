#   Read in table of SGE job info for all tested pipelines and perform stats/
#   visualizations

#   Compute max number of concurrent CPUs used at any time point
bins = 10000
bin_size_time = (max(job_df$end_time) - min(job_df$start_time)) / bins
bin_size = as.numeric(bin_size_time)

print(paste('Using bin size of', bin_size_time))

total_cpus = rep.int(0, bins)
for (i in 1:bins) {
    bin_start = min(job_df$start_time) + (i - 1) * bin_size
    bin_end = bin_start + bin_size
    
    total_cpus[i] = job_df |>
        #   A given job is counted only if it it running throughout the full bin
        filter(start_time <= bin_start, end_time > bin_end) |>
        summarize(total_cpus = sum(slots)) |>
        pull(total_cpus)
}

ggplot(job_df, aes(x = start_time)) + geom_histogram(aes(weight = slots))
