library('jaffelab')
library('getopt')

#  Query "GPU usage" and "GPU memory used" for the installed GPUs to ultimately
#  select a GPU or GPUs to use for alignment (by modifying "gpuMask" in the
#  Arioc config)

spec = matrix(
    c(
        "usage_thres", "u", 1, "integer", "GPU usage cutoff to determine availability",
        "max_gpus", "m", 1, "integer", "Max number of GPUs to use for at once"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

raw_info = system('nvidia-smi -q', intern=TRUE)

#  Parse output of 'nvidia-smi -q' into a matrix, where rows are GPU usage and
#  GPU memory usage (2 rows) and columns are the installed GPUs
usage = sapply(
    grep('Utilization', raw_info),
    function(x) {
        temp = gsub(' ', '', raw_info[(x+1):(x+2)])
        stopifnot(all(ss(temp, ':', 1) == c('Gpu', 'Memory')))
        
        temp = ss(temp, ':', 2)
        temp = as.numeric(ss(temp, '%', 1))
        
        return(temp)
    }
)

#  The integer indices of "available" GPUs
installed_gpus = 0:(ncol(usage) - 1)
good_gpus = installed_gpus[sapply((installed_gpus + 1), function(i) all(usage[,i] < opt$usage_thres))]
good_gpus = good_gpus[1:opt$max_gpus]

if (length(good_gpus) == 0) {
    stop("All GPUs are occupied.")
}

#  Write the available GPU indices to a file, which will be read and parsed
#  into a command that sets CUDA_VISIBLE_DEVICES
writeLines(as.character(good_gpus), con='open_gpus.txt')
