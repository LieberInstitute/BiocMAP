#   For each pipeline tested, extract SGE job details into a tibble. Ultimately
#   form a single tibble of job info and write to CSV. This CSV will be used
#   downstream to compute statistics and perform visualizations.

library('here')

source(here('BiocMAP_benchmark', 'benchmark_stats', 'sge_job_info.R'))

var_names = c(
    'maxvmem', 'ru_maxrss', 'exit_status', 'cpu', 'start_time', 'end_time',
    'slots'
)

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
    job_name = job_info |>
        str_extract('"snakejob.*\\.sh"') |>
        str_remove_all('"snakejob\\.|\\.sh"'),
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
    job_name = log_text[grep("submitted process", log_text)] |>
        ss('process ', 2) |>
        ss(' \\(') |>
        ss(' >'),
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
    job_name = log_text[grep("submitted process", log_text)] |>
        ss('process ', 2) |>
        ss(' \\(') |>
        ss(' >') |>
        ss(':METHYLSEQ:', 2),
    software = 'MethylSeq'
)

################################################################################
#   Extract SGE job info from the combined table
################################################################################

job_df = rbind(job_df_wg, job_df_bio, job_df_methyl)

qacct_list = lapply(
    job_df$job_id,
    function (job_id) get_job_info(
        job_id, month_limit = 3, fail_on_error = FALSE
    )
)

#   Extract variables of interest from the 'qacct' output
for (var_name in var_names) {
    job_df[[var_name]] = sapply(
        qacct_list,
        function(x) {
            if (any(is.na(x))) {
                return(NA)
            } else {
                return(
                    x[grep(paste0('^', var_name), x)] |>
                        str_remove_all(paste0('(', var_name, ' *| *$)'))
                )
            }
        }
    )
}

#   Print atypical jobs (those with NA SGE info or nonzero exit stati)
print('The following jobs had NA SGE info or nonzero exit stati:')
print(job_df |> filter((is.na(exit_status)) | (exit_status != "0")))

#   Take only jobs that finished successfully; extract start and end times in
#   'chron' format
job_df = job_df |>
    filter((!is.na(exit_status)) & (exit_status == "0")) |>
    mutate(
        start_time = to_chron(start_time),
        end_time = to_chron(end_time),
        slots = as.numeric(slots),
        software = as.factor(software)
    )

#   Clean up memory variables: measure in GB
to_numeric_memory = function(mem_string) {
    suffices = mem_string |> str_extract('[[:alpha:]]{2}$')
    coeffs = as.numeric(mem_string |> ss('[[:alpha:]]{2}$'))
    
    mem_num = rep(0, length(mem_string))
    
    for (i in 1:length(mem_string)) {
        if (suffices[i] == 'KB') {
            mem_num[i] = 1e3 * coeffs[i]
        } else if (suffices[i] == 'MB') {
            mem_num[i] = 1e6 * coeffs[i]
        } else if (suffices[i] == 'GB') {
            mem_num[i] = 1e9 * coeffs[i]
        } else if (suffices[i] == 'TB') {
            mem_num[i] = 1e12 * coeffs[i]
        } else {
            stop('Unknown memory unit')
        }
    }
    
    return(mem_num)
} 

job_df$maxvmem |> to_numeric_memory()

write.csv(job_df, out_path, row.names = FALSE, quote = FALSE)