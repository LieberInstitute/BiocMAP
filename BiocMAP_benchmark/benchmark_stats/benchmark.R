#   Read in table of SGE job info for all tested pipelines and perform stats/
#   visualizations
library('tidyverse')
library('here')

job_df_path = here(
    'BiocMAP_benchmark', 'benchmark_stats', 'combined_job_info.csv'
)

stats_df_path = here(
    'BiocMAP_benchmark', 'benchmark_stats', 'benchmark_stats.csv'
)

bin_size_s = 30

metrics = c(
    'Wallclock Duration (Days)' = 'active_wallclock_days',
    'CPU Hours' = 'cpu_hours',
    'Peak Total CPUs' = 'max_concurrent_cpus',
    'TB Hours' = 'total_mem_TB_hours',
    'Peak Total Memory (GB)' = 'max_concurrent_vmem'
)

################################################################################
#   Functions
################################################################################

#   Given job_df, time bins of [[bins_size_s]] seconds, a column name in job_df
#   to sum ('resource_name'), and one of the values in job_df$software
#   ('software_name'), compute the maximal concurrent usage of the given
#   resource and return.
get_max_concurrent_resource = function(
        software_name, job_df, bin_size_s, resource_name
) {
    #   Infer an integer number of time bins given the bin size
    bins = job_df |>
        filter(software == software_name) |>
        summarize(
            bins = as.integer(
                as.numeric(
                    max(end_time) - min(start_time), units = "hours"
                ) * 3600 / bin_size_s
            )
        ) |>
        pull(bins)
    
    #   The time of the start of the first job for this software
    pipeline_start = job_df |>
        filter(software == software_name) |>
        summarize(a = min(start_time)) |>
        pull(a)
    
    current_max = 0
    
    for (i in 1:bins) {
        bin_start = pipeline_start + (i - 1) * bin_size_s
        bin_end = bin_start + bin_size_s
        
        this_total = job_df |>
            #   A given job is counted only if it it running throughout the full
            #   bin
            filter(start_time <= bin_start, end_time > bin_end) |>
            summarize(total = sum(get({{ resource_name }}))) |>
            pull(total)
        
        if (this_total > current_max) {
            current_max = this_total
        }
    }
    
    return(current_max)
}

################################################################################
#   Compute statistics
################################################################################

job_df = read.csv(job_df_path) |>
    as_tibble() |>
    mutate(
        start_time = strptime(start_time, format = "%Y-%m-%d %H:%M:%S"),
        end_time = strptime(end_time, format = "%Y-%m-%d %H:%M:%S"),
        qsub_time = strptime(qsub_time, format = "%Y-%m-%d %H:%M:%S")
    ) |>
    filter(software %in% c('BiocMAP', 'MethylSeq'))

#   Based on https://stackoverflow.com/a/28938694. First, merge overlapping
#   jobs to form the minimal set of disjoint time periods, allowing for 3-minute
#   gaps (Nextflow/Snakemake has a small overhead of processing between
#   submitted jobs). Here, 'qsub_time' is used, since in extreme cases, a job
#   will take a couple hours before starting once submitted.
#
#   Next, sum up the wallclock time for these disjoint intervals to determine
#   the total wallclock time for each pipeline software.
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

#   Compute the total GPU hours used by BiocMAP
biocmap_gpu_hours = job_df |>
    filter(software == "BiocMAP", grepl('AlignReads', jobname)) |>
    #   BiocMAP is configured to use 3 GPUs per process
    summarize(
        gpu_hours = as.numeric(
            sum(end_time - start_time) * 3, units = "hours"
        )
    ) |>
    pull(gpu_hours)

#   Add GPU hours and maximum concurrent CPU and vmem usage for each pipeline
stats_df = tibble(
    software = c('BiocMAP', 'MethylSeq'),
    #   non-BiocMAP pipelines don't use GPUs
    gpu_hours = c(biocmap_gpu_hours, 0),
    max_concurrent_cpus = sapply(
        c('BiocMAP', 'MethylSeq'),
        get_max_concurrent_resource, job_df, bin_size_s, 'slots'
    ),
    max_concurrent_vmem = sapply(
        c('BiocMAP', 'MethylSeq'),
        get_max_concurrent_resource, job_df, bin_size_s, 'maxvmem'
    )
)

#   Compute other statistics
stats_df = job_df |>
    group_by(software) |>
    summarize(
        total_mem_TB_hours = sum(mem) / 1e12 / 3600,
        vmem_frac_used = sum(maxvmem) / sum(requested_h_vmem),
        cpu_hours = sum(as.numeric(cpu |> str_remove('s'))) / 3600,
        practical_wallclock_days = as.numeric(
            max(end_time) - min(start_time), units = "days"
        )
    ) |>
    left_join(wallclock_df) |>
    left_join(stats_df)

print('Stats:')
print(stats_df)

write.csv(stats_df, stats_df_path, quote = FALSE, row.names = FALSE)

################################################################################
#   Visualizations
################################################################################

vis_df = stats_df |>
    select(
        - all_of(c('vmem_frac_used', 'practical_wallclock_days', 'gpu_hours'))
    ) |>
    #   Convert max concurrent vmem to GB
    mutate(max_concurrent_vmem = max_concurrent_vmem / 1e9) |>
    pivot_longer(cols = -software) |>
    #   Re-order metrics and use human-readable names
    mutate(
        name = factor(
            names(metrics)[match(name, metrics)],
            levels = names(metrics),
            ordered = TRUE
        )
    )

#   Bar chart of various key metrics
ggplot(vis_df) +
    geom_col(aes(x = software, y = value, fill = software)) +
    facet_wrap(~name, scales = "free_y", nrow = 2) +
    labs(x = NULL, y = NULL, fill = "Software") +
    theme_bw()
