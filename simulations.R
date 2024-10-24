source("functions.R")
  

########################################################
### Monte Carlo estimation of recovery probabilities ###
########################################################


## Dimension of the process
d <- 10

## Stability index (vector of size d)
alpha <- 3/2

## Marginal masses (2xd matrix) c_1^+,...,c_d^+ and c_1^-,...,c_d^-
c <- matrix(1, nrow = 2, ncol = d)


##########################
### Hüsler--Reiss setup ###
##########################

lower_thr <- 1 
upper_thr <- 2


#############################################
### Parameters for approximate simulation ###
#############################################

## epsilon threshold to be used for approximate simulation
epsilon <- 10^(-3)

## Average number of jumps larger than epsilon
avg_no_jumps <- 10000 #10000


## Function for estimating probability of recovering the tree.
## K is the number of simulations to base each probability estimate on,
## qvec is the values of q to use for the structure learning, and n is
## the number of times the process is sampled (i.e. sampling frequency is 1/n).

## Setting the values to be used
K <- 500
qvec <- seq(.5, .999, .01)
n_range <- c(500, 750, 1000, 2000)

## Estimation of probabilities
results <- matrix(nrow = length(qvec), ncol = length(n_range))
for (j in 1:length(n_range)) {
  print(j)
  matches <- matrix(NA, nrow = K, ncol = length(qvec))
  for (k in 1:K) {
      T0 <- generate_random_tree(d)
      Gamma <- complete_Gamma(Gamma = runif(d - 1, lower_thr, upper_thr), graph = T0)
      data_tmp <- sim_X(n = n_range[j], rate = avg_no_jumps, tree = T0, epsilon = epsilon, Gamma = Gamma, alpha = alpha, c = c)
    for (i in 1:length(qvec)) {
      matches[k, i] <- length(E(graph.intersection(T0, levy_mst(incs = diff(data_tmp), p = qvec[i])))) == d - 1
    }
  }
  results[, j] <- colMeans(matches)
}


custom_labels <- c("S1" = "500", "S2" = "750", "S3" = "1000", "S4" = "2000")

## Load saved data of simulation study:
# load("data/sim_study_d30_K500.Rdata")

dat_sim <- tibble(X = qvec) %>%
  bind_cols(tibble(S1 = results[, 1])) %>%
  bind_cols(tibble(S2 = results[, 2])) %>%
  bind_cols(tibble(S3 = results[, 3])) %>%
  bind_cols(tibble(S4 = results[, 4])) %>%
  pivot_longer(cols = -X, names_to = "sample_size", values_to = "Value")

# Plot
gg_sim <- ggplot(dat_sim, aes(x = X, y = Value, color = sample_size)) +
  geom_line() +
  scale_color_manual(values = c(
    "S1" = my_palette$light_blue,
    "S2" = my_palette$yellow,
    "S3" = my_palette$green,
    "S4" = my_palette$pink
  ), 
  labels =custom_labels) +
  theme_minimal() +
  labs(
    x = TeX("Probability threshold $q$"),
    y = "Recovery proportion",
    color = "Sample size"
  ) +
  theme(
    plot.margin = unit(c(1, 4.5, 1, 1), "lines"),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    panel.border = element_rect(size = 1 / 8, fill = NA),
    axis.title = element_text(size = 13),   # Axis titles
    axis.text = element_text(size = 13),    # Axis ticks
    legend.text = element_text(size = 13),  # Legend text
    legend.title = element_text(size = 13)  # Legend title
  )

gg_sim

#  save(results, file = "sim_study_d10_K500.Rdata")

save_myplot(
  plt = gg_sim,
  plt_nm = paste0(figure_dest_folder, "sim_study_d30", ".pdf"),
  width = 5,
  height = 5,
  cairo = FALSE
)







