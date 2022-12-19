log_msg <- function(msg, verbose = TRUE, log_time = TRUE, sample_id = ""){
    if(verbose == TRUE){
        if(sample_id != ""){
           cat(sample_id)
           cat(": ")
        }
        cat(msg)
        if(log_time == TRUE){
            cat(" - ")
            cat(format(Sys.time(), format="%Y-%m-%d %H:%M:%S", usetz = FALSE))
        }
        cat("\n")
    }
}
