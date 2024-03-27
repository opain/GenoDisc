# Convert opt list into data.frame for log file
opt_to_df<-function(opt){
    opt<-stack(opt)
    names(opt)<-c('Value','Parameter')
    opt<-opt[,c('Parameter','Value')]
    return(opt)
}

# Create log file with standard header including script name and opt
log_header<-function(log_file, opt, script){
    sink(file = log_file, append = F)
        cat0(
            '#################################################################\n',
            '# ',script,'\n',
            '# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)\n',
            '#################################################################\n'
            )
        cat0('---------------\n')
        print.data.frame(opt_to_df(opt), row.names = F, quote = F, right = F)
        cat0('---------------\n')
        cat0('Analysis started at ',as.character(start.time),'\n')
    sink()
}

# Add to an existing log file
log_add<-function(log_file = NULL, message, sep = '\n'){
    if(is.null(log_file)){
        cat(message, sep = sep)
    } else {
        sink(file = log_file, append = T)
            cat(message, sep = sep)
        sink()
    }
}

# cat function with sep = ''
cat0 <- function(..., sep = '', file = "", append = FALSE) {
  cat(..., sep = sep, file = file, append = append)
}
