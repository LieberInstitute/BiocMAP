library('jaffelab')

filenames = system('ls *.log', intern=TRUE)
ids = basename(ss(filenames, '.', fixed=TRUE))
percs = lapply(filenames, function(f) {
  raw_chars = system(paste0('cat ', f, ' | grep "C methylated in C.. context:" | cut -d ":" -f 2'), intern=TRUE)
  return(as.numeric(sub("[:space:]|%", "", raw_chars)))
})
mat = matrix(unlist(percs), dimnames = list(ids, c("perc_M_CpG", "perc_M_CHG", "perc_M_CHH")), byrow=TRUE, ncol=3)

save(mat, file='BME_metrics.rda')