#   Given an SGE job ID 'job_id', search accounting files dated at most
#   [month_limit] months back and return a character vector effectively
#   containing the lines of output from running 'qacct -j [job_id]'

get_job_info = function(job_id, month_limit = 3) {
    acct_dir = '/cm/shared/apps/sge/sge-8.1.9/default/common'
    
    user = system('whoami', intern=TRUE)
    
    #--------------------------------------------------------------------------#
    #  Determine which accounting files should potentially be searched through
    #--------------------------------------------------------------------------#
    
    all_acct_files = list.files(acct_dir, pattern = 'accounting.*\\.txt')
    
    date_now = system('date +%m/%Y', intern=TRUE)
    date_file = ss(all_acct_files, '_', 2)
    
    time_now = as.numeric(ss(date_now, '/', 2)) * 12 + # years from 0 as months
        as.numeric(ss(date_now, '/', 1))        # leftover months
    time_file = as.numeric(substr(date_file, 1, 4)) * 12 + # years from 0 as months
        as.numeric(substr(date_file, 5, 6))        # leftover months
    
    age_file = time_now - time_file # (rounded) age in months of file
    
    valid_acct_files = all_acct_files[age_file <= month_limit]
    valid_acct_files = file.path(acct_dir, valid_acct_files) # make full path
    
    #  The newest accounting file is named differently
    valid_acct_files = c(valid_acct_files, file.path(acct_dir, 'accounting'))
    
    #--------------------------------------------------------------------------#
    #  Form an unparsed vector of job information combined across relevant files
    #--------------------------------------------------------------------------#
    
    job_info = c()
    i = 1
    while((length(job_info) == 0) && i <= length(valid_acct_files)) {
        acct_file = valid_acct_files[i]
        
        #  Since an error is returned if the job is not found, we check that
        #  manually before running 'qacct'
        cmd = paste(
            'grep', user, acct_file, '| cut -d ":" -f 6 | grep', job_id,
            '| wc -l'
        )
        num_lines = system(cmd, intern=TRUE)
        
        if (num_lines > 0) {
            cmd = paste0('qacct -j ', job_id, ' -f ', acct_file, ' -u ', user)
            job_info = system(cmd, intern=TRUE)
        }
        
        i = i + 1
    }
    
    if (length(job_info) == 0) {
        stop("No job found in this time period for this ID.")
    }
    
    return(job_info)
}
