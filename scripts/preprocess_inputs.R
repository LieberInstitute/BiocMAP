## Required libraries
library('devtools')
library('jaffelab')

###################################################
#  Functions
###################################################

run_command = function(command) {
    print(paste0("Running command: '", command, "'..."))
    system(command)
    print("Done.")
}

get_value = function(rules, key) {
    this_line = rules[grep(key, ss(rules, '='))]
    if (length(this_line) != 1) {
        stop(paste0("Key '", key, "' had ", length(this_line), " values in rules.txt."))
    }
    
    return(gsub(' ', '', ss(this_line, '=', 2)))
}

form_links = function(rules, key, ids, end_name) {
    value = get_value(rules, key)
    if (value != "NA") {
        pieces = ss(value, '[id]', fixed=TRUE)
        
        #  Get the paths, substituting in actual ids for each '[id]'
        paths = ''
        for (i in 1:(length(pieces)-1)) {
            paths = paste0(pieces[i], paths)
        }
        paths = paste0(paths, pieces[length(pieces)])
        stopifnot(length(paths) == length(ids))
        
        #  Symbolically link each file into the current working directory
        for (i in 1:length(ids)) {
            command = paste0('ln -s ', paths[i], ' ', ids[i], end_name)
            run_command(command)
        }
    }
}

###################################################
#  Main
###################################################

rules = readLines('rules.txt')

man_path = file.path(get_value(rules, 'base_dir'), get_value(rules, 'manifest'))
manifest = read.table(man_path, header = FALSE, stringsAsFactors = FALSE)
ids = manifest[ncol(manifest)]

#  Arioc SAMs and logs
form_links(rules, 'sam', ids, '.sam')
form_links(rules, 'arioc_log', ids, '_arioc.log')

#  XMC logs
form_links(rules, 'xmc_log', ids, '_xmc.log')

#  BME logs
form_links(rules, 'bme_log', ids, '_bme.log')

#  Trim Galore reports
form_links(rules, 'trim_report_r1', ids, '_trim_report_r1.txt')
form_links(rules, 'trim_report_r2', ids, '_trim_report_r2.txt')
