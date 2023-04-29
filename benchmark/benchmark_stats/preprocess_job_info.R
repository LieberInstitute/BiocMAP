#   For each pipeline tested, extract SGE job details into a tibble. Ultimately
#   form a single tibble of job info and write to CSV. This CSV will be used
#   downstream to compute statistics and perform visualizations.

library('here')
library('sgejobs')

source(here('BiocMAP_benchmark', 'benchmark_stats', 'sge_job_info.R'))

log_path_wg = here('BiocMAP_benchmark', 'wg_blimp', 'logs', 'run.log')
log_path_bio_first = here('BiocMAP_benchmark', 'BiocMAP_dir', '.nextflow.log.1')
log_path_bio_second = here('BiocMAP_benchmark', 'BiocMAP_dir', '.nextflow.log')
log_path_methyl = here('BiocMAP_benchmark', 'methylseq', '.nextflow.log')

out_path = here('BiocMAP_benchmark', 'benchmark_stats', 'combined_job_info.csv')

################################################################################
#   Read in job info for wg-blimp
################################################################################

log_text = readLines(log_path_wg)

job_info = log_text[grep('^Submitted job .* with external jobid', log_text)] |>
    str_extract(" 'Your job [0-9]+ \\(.*\\) ")

job_df_wg = tibble(
    job_id = job_info |>
        str_extract('job [0-9]+ \\(') |>
        str_remove_all('[ job\\(]'),
    software = 'wg-blimp'
)

################################################################################
#   Read in job info for BiocMAP
################################################################################

log_text = c(readLines(log_path_bio_first), readLines(log_path_bio_second))

job_df_bio = tibble(
    job_id = log_text[grep("submitted process", log_text)] |>
        ss('> jobId: ', 2) |>
        ss(';', ),
    software = 'BiocMAP'
)

################################################################################
#   Read in job info for methylseq
################################################################################

log_text = readLines(log_path_methyl)

job_df_methyl = tibble(
    job_id = log_text[grep("submitted process", log_text)] |>
        ss('> jobId: ', 2) |>
        ss(';', ),
    software = 'MethylSeq'
)

################################################################################
#   Extract SGE job info from the combined table
################################################################################

job_df = rbind(job_df_wg, job_df_bio, job_df_methyl)

job_files = sapply(
    job_df$job_id, get_job_file, month_limit = 3, fail_on_error = FALSE
)

print('The following jobs could not be found:')
print(job_df[is.na(job_files),])

job_df = job_df[!is.na(job_files),]

job_df = accounting_read(
        job_df$job_id, accounting_files = job_files[!is.na(job_files)]
    ) |>
    accounting_parse() |>
    cbind(job_df) |>
    as_tibble()

print('The following jobs had nonzero exit stati:')
job_df |>
    filter(exit_status != "0") |> 
    select(all_of(c('job_id', 'software'))) |>
    print()

job_df = job_df |>
    #   Take only jobs that finished successfully
    filter(exit_status == "0") |>
    mutate(
        software = as.factor(software),
        #   Get requested mem_free and h_vmem. Note this code assumes these two
        #   things are always requested, and that the amount is always in GB!
        requested_h_vmem = category |>
            str_extract('h_vmem=[0-9]+G') |>
            str_remove_all('h_vmem=|G') |>
            as.numeric() * 1e9,
        requested_mem_free = category |>
            str_extract('mem_free=[0-9]+G') |>
            str_remove_all('mem_free=|G') |>
            as.numeric() * 1e9
    ) |>
    #   Remove duplicate columns
    select(!all_of(c('input_id', 'jobnumber', 'category')))



write.csv(job_df, out_path, row.names = FALSE, quote = FALSE)
