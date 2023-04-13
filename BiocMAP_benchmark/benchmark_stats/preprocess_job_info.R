#   For each pipeline tested, extract SGE job details into a tibble. Ultimately
#   form a single tibble of job info and write to CSV. This CSV will be used
#   downstream to compute statistics and perform visualizations.

library('here')

source(here('BiocMAP_benchmark', 'benchmark_stats', 'sge_job_info.R'))

var_names = c(
    'maxvmem', 'ru_maxrss', 'exit_status', 'cpu', 'start_time', 'end_time',
    'slots'
)

log_path = here('BiocMAP_benchmark', 'wg_blimp', 'logs', 'run.log')
log_text = readLines(log_path)

job_info = log_text[grep('^Submitted job .* with external jobid', log_text)] |>
    str_extract(" 'Your job [0-9]+ \\(.*\\) ")

job_df = tibble(
    job_name = job_info |>
        str_extract('"snakejob.*\\.sh"') |>
        str_remove_all('"snakejob\\.|\\.sh"'),
    job_id = job_info |>
        str_extract('job [0-9]+ \\(') |>
        str_remove_all('[ job\\(]')
)

a = lapply(
    job_df$job_id,
    function (job_id) get_job_info(
        job_id, month_limit = 3, fail_on_error = FALSE
    )
)

#   Extract variables of interest from the 'qacct' output
for (var_name in var_names) {
    job_df[[var_name]] = sapply(
        a,
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

#   Take only jobs that finished successfully; extract start and end times in
#   'chron' format
job_df = job_df |>
    filter((!is.na(exit_status)) & (exit_status == "0")) |>
    mutate(
        start_time = to_chron(start_time),
        end_time = to_chron(end_time),
        slots = as.numeric(slots)
    )
