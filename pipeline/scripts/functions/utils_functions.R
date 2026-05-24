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

# Create function to read in parameters in the config file
read_param <- function(config, param, return_obj = T, quiet = F){
  library(yaml)
  
  # Read in the config file
  config_file <- read_yaml(config)
  
  if(all(names(config_file) != param)){
    # Check default config file. Resolve via the global option `pipeline_dir`
    # (set by each script after parse_args) so this works when CWD is not the
    # pipeline folder; fall back to CWD-relative for backward compatibility.
    .default_config_path <- if (!is.null(getOption('pipeline_dir'))) {
      file.path(getOption('pipeline_dir'), 'config.yaml')
    } else {
      'config.yaml'
    }
    config_file <- read_yaml(.default_config_path)
    
    if(all(names(config_file) != param)){
      if(quiet == F){
        cat(param, 'parameter is not present in user specified config file or default config file.\n')
      }
      return(NULL)
    } else {
      if(quiet == F){
        cat(param, 'parameter is not present in user specified config file, so will use value in default config file.\n')
      }
    }
  }
  
  # Identify value for param
  file <- config_file[[param]]
  file[file == 'NA']<-NA
  
  # If resdir, and NA, set to 'resources'
  if(param == 'resdir'){
    if(is.na(file)){
      file <- 'resources'
    }
  }
  
  # If refdir, and NA, set to '<resdir>/data/ref'
  if(param == 'refdir'){
    if(is.na(file)){
      resdir <- read_param(config = config, param = 'resdir', return_obj = F, quiet = quiet)
      file <- paste0(resdir, '/data/ref')
    }
  }
  
  if(return_obj){
    if(!is.na(file)){
      obj <- fread(file)
    } else {
      obj <- NULL
    }
    return(obj)
  } else {
    file <- file[order(file)]
    return(file)
  }
}
# Make a function to source all files in a directory
source_all <- function(directory) {
  # List all .R files in the specified directory
  r_files <- list.files(directory, pattern = "\\.R$", full.names = TRUE)
  
  # Source each file
  for (file in r_files) {
    source(file)
  }
}
