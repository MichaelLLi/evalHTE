# function to call Julia and run optimization
#' @import JuliaCall
#' @importFrom JuliaCall julia_setup
#' @importFrom JuliaCall julia_source
#' @importFrom JuliaCall julia_eval
run_optimization <- function() {
     tryCatch({
       if (!requireNamespace("JuliaCall", quietly = TRUE)) {
         stop("JuliaCall package is required but not installed.")
       }
       
       julia_path <- "/Users/alana/.julia/juliaup/julia-1.10.5+0.aarch64.apple.darwin14/bin/julia"
       if (!file.exists(julia_path)) {
         stop("Julia executable not found at the specified path.")
       }
       
       JuliaCall::julia_setup(JULIA_HOME = dirname(dirname(julia_path)))
       
       # Source the Julia file
       JuliaCall::julia_source("optimization.jl")
       
       # Run the optimization
       result <- JuliaCall::julia_eval("run_optimization()")
       
       # Extract and return results as a list
       list(
         alphas = result[[1]],
         new_law = result[[2]],
         normal_law = result[[3]]
       )
     }, error = function(e) {
       warning("Failed to run Julia optimization: ", e$message)
       return(NULL)
     })
   }


# helper function to plot the GATE estimates
gate_ggplot <- function(data) {
    ggplot(data, aes(
      x = .data$group, y = .data$estimate,
      ymin = .data$lower , ymax = .data$upper, color = .data$algorithm)) +
    ggdist::geom_pointinterval(
      width = 0.5,
      position = position_dodge(0.5),
      interval_size_range = c(0.8, 1.5),
      fatten_point = 2.5) +
    theme_bw() +
    # add color theme
    scale_color_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank()) +
    labs(x = "Group", y = "GATE estimate") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#4e4e4e")
}


# rename the columns of the data frame with the interaction terms
rename_interaction_terms <- function(interaction_df){
  colnames(interaction_df) <- gsub(":", "_", colnames(interaction_df))
  colnames(interaction_df) <- gsub("\\*", "_", colnames(interaction_df))
  colnames(interaction_df) <- gsub("\\(", "", colnames(interaction_df))
  colnames(interaction_df) <- gsub("\\)", "", colnames(interaction_df))
  colnames(interaction_df) <- gsub("\\+", "_", colnames(interaction_df))
  return(interaction_df)
}



# function to convert formula and create new variables
convert_formula <- function(user_formula, data, treatment){

  # get the outcome variable name from the formula
  outcome <- all.vars(user_formula)[1]

  # get the covariates from the formula
  interaction_df <- model.matrix(user_formula, data)
  interaction_df <- rename_interaction_terms(interaction_df)

  # remove variable Intercept from covariates list by name
  covariates_vec <- colnames(interaction_df)
  covariates_vec <- covariates_vec[!covariates_vec %in% c("Intercept", paste0(treatment))]
  # combine the interaction_df with the original data frame
  new_data = data %>% dplyr::select(all_of(outcome))
  combined_data <- cbind(new_data, interaction_df)

  return(list(data = combined_data, covariates = covariates_vec, outcome = outcome))
}

