handle_message <- function(cnd){cat("\n\n######## MESSAGE ###########\n\n");tree <- rlang::trace_back(); print(cnd); print(tree);cat("\n\n###################\n\n")}
handle_warning <- function(cnd){cat("\n\n######## WARNING ###########\n\n");tree <- rlang::trace_back(); print(cnd); print(tree);cat("\n\n###################\n\n")}

suppress_certain_warning <- function(.expr, .f, ...) {
  eval.parent(substitute(
  withCallingHandlers( .expr, warning = function(w) {
    cm <- conditionMessage(w)
    if(is.character(.f)){
      cond_vec <- logical(length(.f))
      for (i in seq_len(length(.f))){
        cond_vec[i] <- grepl(.f[i], cm)
      }
      cond <- any(cond_vec)
    }else{
      cond <- rlang::as_function(.f)(cm,...)
    }
    if (cond) {
      invokeRestart("muffleWarning")
    }
  })
  ))
}

suppress_certain_message <- function(.expr, .f, ...) {
  eval.parent(substitute(
  withCallingHandlers( .expr, message = function(w) {
    cm <- conditionMessage(w)
    cond <- 
      if(is.character(.f)){
      cond_vec <- logical(length(.f))
      for (i in seq_len(length(.f))){
        cond_vec[i] <- grepl(.f[i], cm)
      }
      cond <- any(cond_vec)
    }else{
      cond <- rlang::as_function(.f)(cm,...)
    }
    if (cond) {
      invokeRestart("muffleMessage")
    }
  })
  ))
}