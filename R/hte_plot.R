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
    lower = gate - 1.96*sd,
    upper = gate + 1.96*sd,
    algorithm = gsub("_", " ", alg)) -> data

# plot GATE estimates
out = gate_ggplot(data)

return(out)

}


#' Plot the uniform confidence interval
#' @import ggplot2
#' @import ggthemes
#' @importFrom stats sd
#' @importFrom rlang .data
#' @importFrom scales percent
#' @importFrom tidyr pivot_longer
#' @param x An object of \code{evaluate_hte()} class. This is typically an output of \code{evaluate_hte()} function.
#' @param ... Further arguments passed to the function.
#' @return A plot of ggplot2 object.
#' @export
plot_CI.hte <- function(x, ...){

# parameters
estimate = x
estimate_algs = estimate$out_algs
estimate_user = estimate$out_user
data_algs = tibble()
data_user = tibble()

# -----------------------------------------
# estimate HTE from ML algorithms
# -----------------------------------------

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

    # graphLabels <- data.frame(
    #   type = algorithms,
    #   Pval = map(
    #     fit$AUPEC, ~.x[c('aupec', 'sd')]) %>%
    #     bind_rows() %>%
    #     mutate(Pval = paste0("AUPEC = ", round(aupec, 2), " (s.e. = ", round(sd, 2), ")")) %>% pull(Pval))

    Tcv = estimate_algs$estimates[['Tcv']] %>% as.numeric()
    Ycv = estimate_algs$estimates[['Ycv']] %>% as.numeric()

    map(fit$URATE, ~.x) %>%
      bind_rows() %>%
      mutate(
        RATEmin = rate - 1.2*sd - 0.68*sd[length(sd)]*length(sd)/seq(1, length(sd)),
        RATEpoint = rate - 1.67*sd, 
        fraction = rep(seq(1,length(Ycv))/length(Ycv), length(algorithms)),
        type = lapply(algorithms, function(x)rep(x,length(Ycv))) %>% unlist) %>%
    rename(
    `GATE estimate` = rate, 
    `Minimum-area lower band` = RATEmin,
    `Pointwise lower band` =RATEpoint) %>%
    tidyr::pivot_longer(
      ., cols = c(
        "GATE estimate", "Minimum-area lower band", "Pointwise lower band"),
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

  # graphLabels <- data.frame(
  #   type = "user-defined HTE",
  #   Pval = map(
  #     fit$AUPEC, ~.x[c('aupec', 'sd')]) %>%
  #     bind_rows() %>%
  #     mutate(Pval = paste0("AUPEC = ", round(aupec, 2), " (s.e. = ", round(sd, 2), ")")) %>% pull(Pval))

  fit$AUPEC %>%
    bind_rows() %>%
    mutate(
      RATEmin = rate - 1.2*sd - 0.68*sd[length(sd)]*length(sd)/seq(1, length(sd)),
      RATEpoint = rate - 1.67*sd, 
      fraction = rep(seq(1,length(Ycv))/length(Ycv),1),
      type = lapply("user-defined", function(x)rep(x,length(Ycv))) %>% unlist) %>%
    rename(
    `GATE estimate` = rate, 
    `minimum-area lower band` = RATEmin,
    `pointwise lower band` =RATEpoint) %>%
    tidyr::pivot_longer(
      ., cols = c(
        "GATE estimate", "Minimum-area lower band", "Pointwise lower band"),
      names_to = "Type",
      values_to = "value"
    ) -> data_user
}

# dataframe for plotting
data <- bind_rows(data_algs, data_user) 

# plot   
ggplot(data, aes(x=fraction,y=value)) +
  geom_line(alpha=0.8, aes(color = Type)) +
  scale_colour_few("Dark")+
  xlab("Maximum Proportion Treated")+
  ylab("GATE Estimates")+
  facet_wrap(~type)+
  scale_x_continuous(labels=scales::percent)+
  scale_y_continuous(
    limits = c(quantile(data["Minimum-area lower band",], 0.01, na.rm = TRUE), quantile(data["Minimum-area lower band",],  0.99, na.rm = TRUE)))+
  theme_few()+
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "dotted") +  
  theme(
    legend.position = "right",
    text = element_text(size=13.5),
        axis.text = element_text(size=10),
        strip.text = element_text(size = 13.5)) -> out

  return(out)
}