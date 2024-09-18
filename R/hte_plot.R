#' Plot the GATE estimate
#' @import ggplot2
#' @import ggthemes
#' @importFrom stats sd
#' @importFrom rlang .data
#' @param x An table object. This is typically an output of \code{evaluate_hte()} function.
#' @param ... Further arguments passed to the function.
#' @importFrom ggplot2 .data
#' @importFrom ggdist geom_pointinterval
#' @return A plot of ggplot2 object.
#' @export
plot.hte <- function(x, ...){

# parameters
estimate = x

# format output
bind_rows(estimate$out_algs$qoi$GATE) %>%
  mutate(
    estimate = gate,
    lower = gate - qnorm(0.975)*sd,
    upper = gate + qnorm(0.975)*sd,
    algorithm = gsub("_", " ", alg)) -> data

# plot GATE estimates
out = gate_ggplot(data)

return(out)

}

#' @export
plot_CI <- function(x, ...) {
  UseMethod("plot_CI")
}

#' Plot the uniform confidence interval 
#' @import ggplot2
#' @import ggthemes
#' @importFrom stats sd
#' @importFrom rlang .data
#' @importFrom scales percent
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map
#' @param x An object of \code{evaluate_hte()} class. This is typically an output of \code{evaluate_hte()} function.
#' @param alpha Significance level. Default is 0.05.
#' @param ... Further arguments passed to the function.
#' @return A plot of ggplot2 object.
#' @export
plot_CI.hte <- function(
  x, 
  alpha = 0.05,
  ...){

# parameters
estimate = x
estimate_algs = estimate$out_algs
estimate_user = estimate$out_user
data_algs = tibble()
data_user = tibble()

# -----------------------------------------
# estimate HTE from ML algorithms
# -----------------------------------------

# run optimization in Julia
# results <- run_optimization()
results <- test

if (is.null(results)) {
  warning("Julia optimization failed. Using default values.")
  min_uniform_beta_0 <- 1.2
  min_uniform_beta_1 <- 0.68
  min_pointwise_score <- 1.96
} else {

  # get the optimal parameters
  min_uniform_beta_0 <- results$beta_0[results$alphas == alpha]
  min_uniform_beta_1 <- results$beta_1[results$alphas == alpha]
  min_pointwise_score <- results$normal_law[results$alphas == alpha]
}

# get the estimate from ML algorithms
if(length(estimate_algs) != 0){

  # parameters
  fit = estimate_algs$qoi
  cv = estimate_algs$cv

  # format output under cross validation -----------------------------------------
  if(cv == TRUE){
    print("Not supported under cross-validation")
  }

  # format output under sample splitting -----------------------------------------
  if(cv == FALSE){

    # parameters
    data = estimate_algs$df$data
    algorithms = estimate_algs$df$algorithms


    Tcv = estimate_algs$estimates[['Tcv']] %>% as.numeric()
    Ycv = estimate_algs$estimates[['Ycv']] %>% as.numeric()

    purrr::map(fit$URATE, ~.x) %>%
      bind_rows() %>% 
      mutate(
        RATEmin = rate - min_uniform_beta_1*sd - min_uniform_beta_0*sd[length(sd)]*length(sd)/seq(1, length(sd)),
        # get the z-score 
        z_alpha = qnorm(1-alpha),
        RATEpoint = rate - z_alpha*sd, 
        # RATEpoint = rate - min_pointwise_score*sd, 
        fraction = rep(seq(1,length(Ycv))/length(Ycv), length(algorithms)),
        type = lapply(algorithms, function(x)rep(x,length(Ycv))) %>% unlist
  ) %>%
    rename(
    `GATE estimate` = rate, 
    `Uniform lower band` = RATEmin,
    `Pointwise lower band` =RATEpoint) %>%
    tidyr::pivot_longer(
      ., cols = c(
        "GATE estimate", "Uniform lower band", "Pointwise lower band"),
      names_to = "Type",
      values_to = "value"
    ) -> data_algs
  }
}

# -----------------------------------------
# get HTE from the user-defined function
# -----------------------------------------
if(length(estimate_user) != 0){

   # parameters
  fit = estimate_user$qoi
  cv = estimate_user$cv

  Tcv = estimate_user$estimates[['Tcv']] %>% as.numeric()
  Ycv = estimate_user$estimates[['Ycv']] %>% as.numeric()

  fit$AUPEC %>%
    bind_rows() %>%
    mutate(
      RATEmin = rate - min_uniform_beta_1*sd - min_uniform_beta_0*sd[length(sd)]*length(sd)/seq(1, length(sd)),
      # get the z-score of alpha 
      z_alpha = qnorm(1-alpha),
      RATEpoint = rate - z_alpha*sd,
      # RATEpoint = rate - min_pointwise_score*sd, 
      fraction = rep(seq(1,length(Ycv))/length(Ycv), 1),
      type = lapply("user-defined", function(x)rep(x,length(Ycv))) %>% unlist) %>%
    rename(
    `GATE estimate` = rate, 
    `Uniform lower band` = RATEmin,
    `Pointwise lower band` =RATEpoint) %>%
    tidyr::pivot_longer(
      ., cols = c(
        "GATE estimate", "Uniform lower band", "Pointwise lower band"),
      names_to = "Type",
      values_to = "value"
    ) -> data_user
}

# dataframe for plotting
data <- bind_rows(data_algs, data_user)

# plot   
ggplot(data, aes(x=fraction, y=value)) +
  geom_line(alpha=0.8, aes(color = Type)) +  
  scale_colour_few("Dark") +
  scale_linewidth_continuous(range = c(0.5, 1.5)) + 
  xlab("Maximum Proportion Treated") +
  ylab("GATE Estimates") +
  facet_wrap(~type) +
  scale_x_continuous(labels=scales::percent) +
  theme_few() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dotted") + 
  theme(
    legend.position = "right",
    text = element_text(size=13.5),
    axis.text = element_text(size=10),
    strip.text = element_text(size = 13.5)
  ) -> out
  return(out)
}
