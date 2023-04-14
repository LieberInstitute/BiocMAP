#   Read in table of SGE job info for all tested pipelines and perform stats/
#   visualizations


job_df_path = here(
    'BiocMAP_benchmark', 'benchmark_stats', 'combined_job_info.csv'
)

job_df = read.csv(job_df_path) |>
    as_tibble() |>
    mutate(
        start_time = strptime(start_time, format = "%Y-%m-%d %H:%M:%S"),
        end_time = strptime(end_time, format = "%Y-%m-%d %H:%M:%S"),
        qsub_time = strptime(qsub_time, format = "%Y-%m-%d %H:%M:%S")
    )

################################################################################
#   Stats and visualizations
################################################################################

#   Based on https://stackoverflow.com/a/28938694. First, merge overlapping
#   jobs to form the minimal set of disjoint time periods, allowing for 3-minute
#   gaps (Nextflow/Snakemake has a small overhead of processing between
#   submitted jobs). Here, 'qsub_time' is used, since in extreme cases, a job
#   will take a couple hours before starting once submitted.
#
#   Next, sum up the wallclock time for these disjoint intervals to determine
#   the total wallclock time for each pipeline software.
print('Total wallclock times:')
wallclock_df = job_df |>
    mutate(end_time = end_time + 180) |>
    arrange(software, qsub_time, end_time) |>
    group_by(software) |>
    mutate(
        index = c(
            0,
            cumsum(
                as.numeric(lead(qsub_time)) > cummax(as.numeric(end_time))
            )[-n()]
        )
    ) |>
    group_by(software, index) |>
    summarize(
        start_time = first(start_time), end_time = last(end_time) - 180
    ) |>
    #   Now add up disjoint intervals to compute total wallclock time for each
    #   pipeline software
    group_by(software) |>
    summarize(
        active_wallclock_days = as.numeric(
            sum(end_time - start_time), units = "days"
        )
    ) |>
    ungroup()

#   Compute other statistics
stats_df = job_df |>
    group_by(software) |>
    summarize(
        total_mem_TB_hours = sum(mem) / 1e12 / 3600,
        vmem_frac_used = sum(maxvmem) / sum(requested_h_vmem),
        cpu_hours = sum(as.numeric(cpu |> str_remove('s'))) / 3600,
        practical_wallclock_days = max(end_time) - min(start_time)
    ) |>
    left_join(wallclock_df)

print('Stats:')
print(stats_df)

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
