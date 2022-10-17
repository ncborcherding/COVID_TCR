plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}

get.medians<- function(input, group2) {
  input <- cbind.data.frame(group2, input)
  num <- ncol(input)-1
  med <- input %>%
    group_by(group2) %>%
    dplyr::summarise(across(seq_len(all_of(num)), median))
  med <- as.data.frame(med[,seq_len(num) + 1])
  rownames(med) <- paste0("median.", unique(group2))
  med <- as.data.frame(t(med))
  return(med)
}

normalize <- function(x){
  (x- min(x)) /(max(x)-min(x))
}

"%!in%" <- Negate("%in%")