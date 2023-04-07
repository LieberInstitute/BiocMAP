library('here')
library('tidyverse')

source(here('BiocMAP_benchmark', 'benchmark_stats', 'sge_job_info.R'))

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
