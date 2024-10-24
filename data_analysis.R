######################
### Financial Data ###
######################

source("functions.R")
## Loading data 
log_prices <- read.csv("data/data_all.csv")

#Threshold
p <- 0.95

## Look at all stocks
full_tree <- levy_mst(diff(data.matrix(log_prices[,-1])), p) 
plot(full_tree, asp = 0, vertex.size = 4, vertex.label = names(log_prices[,-1]))
# zm()


#############################
### Selecting some stocks ###
#############################

selected_tickers <- c("JPM", "WFC", "BAC", "USB", "KO", "PEP", "MDLZ", "PG", "CVX", "OXY", "XOM", "APA", "AAPL", "MSFT", "CSCO", "TXN")

selected_colors <- c(rep(my_palette$light_blue, times=4), rep(my_palette$yellow, times=4), rep(my_palette$green, times=4), rep(my_palette$pink, times=4))

selected_data <- log_prices[selected_tickers]
selected_names <- names(selected_data)
selected_incs <- diff(data.matrix(selected_data))
names(selected_incs) <- selected_names
n <- nrow(selected_incs)
d <- ncol(selected_incs)



#####################
### Tree learning ###
#####################

## Estimating and plotting the tree
tree <- levy_mst(selected_incs, p)

layout_fixed <- layout_with_kk(tree)

pdf(file = "figures/tree16.pdf", width = 10, height = 10)
plot(tree, 
     layout = layout_fixed, 
     asp = 0, 
     vertex.size = 20,            # Increase node size
     vertex.label = selected_names, 
     vertex.label.cex = 1.5,      # Increase label size
    edge.width = 3         
)
dev.off()

pdf(file = "figures/tree16_sectors.pdf", width = 10, height = 10)
plot(tree, 
     layout = layout_fixed, 
     vertex.color= selected_colors,
     asp = 0, 
     vertex.size = 20,            # Increase node size
     vertex.label = selected_names, 
     vertex.label.cex = 1.5,      # Increase label size
     edge.width = 3          
)
dev.off()


chi_hat <- emp_chi_levy(data = selected_incs, p=p)
Gamma_hat <- chi2Gamma(chi_hat)



Gamma_tree <- complete_Gamma(Gamma = Gamma_hat, graph = tree)
chi_tree <- Gamma2chi(Gamma_tree)
plot(chi_hat, chi_tree, xlim=c(0,1), ylim=c(0,1))
abline(0,1)

## Estimation of the probabilities on the the quadrant (asymmetry) of each edge
proba_quad <- emp_quad_prop(data = selected_incs, tree = tree, p=.8)


edge_list <- as_edgelist(tree, names = FALSE)  # Returns edges as a 2-column matrix


data_hist <- data.frame(value = proba_quad[edge_list])

gg_hist <- ggplot(data_hist, aes(x = value)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.02, color = "black", fill = "blue", alpha = 0.4) +
  labs(
    title = "",
    x = "",
    y = "Density"
  ) +
  theme_minimal() 

gg_hist

#######################
## Plotting 3 stocks
######################
stock1 <- "USB"
stock2 <- "WFC"
stock3 <- "CSCO"


custom_labels <- c("S1" = "U.S. Bancorp", "S2" = "Wells Fargo", "S3" = "Cisco Systems")


dat_levy <- tibble(X = log_prices[, "date"]) %>%
  bind_cols(tibble(S1 = exp(selected_data[, stock1]))) %>%
  bind_cols(tibble(S2 = exp(selected_data[, stock2]))) %>%
  bind_cols(tibble(S3 = exp(selected_data[, stock3]))) %>%
  pivot_longer(cols = -X, names_to = "Stock", values_to = "Value")

# Plot
gg1 <- ggplot(dat_levy, aes(x = as.Date(X), y = Value, color = Stock)) +
  geom_line() +
  scale_color_manual(values = c(
    "S1" = my_palette$light_blue,
    "S2" = my_palette$yellow,
    "S3" = my_palette$green
  ), 
  labels =custom_labels) +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Price",
    color = "Stock"
  ) +
  theme(
    plot.margin = unit(c(1, 4.5, 1, 1), "lines"),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    # panel.border = element_blank(),
    panel.border = element_rect(size = 1 / 8, fill = NA),
    axis.title = element_text(size = 13),   # Axis titles
    axis.text = element_text(size = 13),    # Axis ticks
    legend.text = element_text(size = 13),  # Legend text
    legend.title = element_text(size = 13)  # Legend title
  )

gg1

save_myplot(
  plt = gg1,
  plt_nm = paste0(figure_dest_folder, "stocks", ".pdf"),
  width = 5,
  height = 5,
  cairo = FALSE
)

############
## Simulation of three stocks with empirical marginals
############



set.seed(45)
Xsim <- matrix(data.matrix(selected_data)[1,], nrow=nrow(selected_data), ncol =d, byrow=TRUE) + sim_X(n = nrow(selected_data), rate = 10*nrow(selected_data), epsilon = 1e-3, Gamma = Gamma_hat, alpha = 1, c = 1, tree = tree, proba_quad = proba_quad, emp_incs = selected_incs)


i <- which(selected_names == stock1)
j <- which(selected_names == stock2)
k <- which(selected_names == stock3)

dat_sim <- tibble(X = 1:nrow(log_prices)) %>%
  bind_cols(tibble(X1 = exp(Xsim[,i]))) %>%
  bind_cols(tibble(X2 = exp(Xsim[,j]))) %>%
  bind_cols(tibble(X3 = exp(Xsim[,k]))) %>%
  pivot_longer(cols = -X, names_to = "Stock", values_to = "Value")

# Plot
gg2 <- ggplot(dat_sim, aes(x = X, y = Value, color = Stock)) +
  geom_line() +
  scale_color_manual(values = c(
    "X1" = my_palette$light_blue,
    "X2" = my_palette$yellow,
    "X3" = my_palette$green
  )) +
  theme_minimal() +
  labs(
    x = "Index",
    y = "Value",
    color = "Coordinate"
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

gg2


save_myplot(
  plt = gg2,
  plt_nm = paste0(figure_dest_folder, "levy_simu", ".pdf"),
  width = 5,
  height = 5,
  cairo = FALSE
)

sim_incs <- diff(Xsim)
colnames(sim_incs) <- selected_names


g1 <- plot_stocks(data=selected_incs, stock1=stock1, stock2, save_plot = FALSE, suffix = "orig")
g2 <- plot_stocks(data=sim_incs, stock1=stock1, stock2, save_plot = FALSE, suffix = "simu")
g3 <- plot_stocks(data=selected_incs, stock1=stock1, stock3, save_plot = FALSE, suffix = "orig")
g4 <- plot_stocks(data=sim_incs, stock1=stock1, stock3, save_plot = FALSE, suffix = "simu")
g5 <- plot_stocks(data=selected_incs, stock1=stock2, stock3, save_plot = FALSE, suffix = "orig")
g6 <- plot_stocks(data=sim_incs, stock1=stock2, stock3, save_plot = FALSE, suffix = "simu")

ggsave(
  filename = paste0(figure_dest_folder, "incs_scatter", ".pdf"),
  plot = grid.arrange(g1, g3, g5, g2, g4, g6, ncol = 3),
  width = 15,
  height = 10
)





###################
### Subsampling ###
###################

set.seed(1245)
## Some setup
sub_size <- round(n / 2) #Size of random subsamples
N <- 300 #Number of subsampling iterations
edge_counts <- matrix(0, nrow = d, ncol = d)
q_range <- c(0.925, 0.95) #q is sampled uniformly in this range

## Function for updating edge_counts

index_mat <- replicate(N, sample(1:n, size = sub_size, replace = FALSE))
q_vec <- runif(N, min = q_range[1], max = q_range[2])

for(k in 1:N){
  tree_tmp <- levy_mst(incs = selected_incs[index_mat[,k],], p = q_vec[k])
  edges <- apply(as_edgelist(tree_tmp),1,toString)
  for (i in 1:(d-1)) {
    for (j in (i+1):d) {
      edge_counts[i,j] <- edge_counts[i,j] + paste(i, j, sep = ", ") %in% edges
    }
  }
}


## Plotting fully connected graph with edge widths proportional to edge_counts
weights_fully <- matrix(nrow=0,ncol=3)
for (i in 1:(d-1)) {
  for (j in (i+1):d) {
    weights_fully <- rbind(weights_fully, c(i, j, edge_counts[i,j]))
  }
}


graph_fully <- graph.data.frame(weights_fully[,1:2],directed=FALSE)
E(graph_fully)$weight <- weights_fully[,3]

non_positive_weights <- weights_fully[,3] <= 0
if (any(non_positive_weights)) {
  print("Non-positive weights found:")
  print(weights_fully[non_positive_weights, 3])
  positive_weights <- weights_fully[,3] > 0
  graph_fully <- subgraph.edges(graph_fully, E(graph_fully)[positive_weights])
  weights_fully <- weights_fully[positive_weights,]
  E(graph_fully)$weight <- weights_fully[,3]
  
}



pdf(file = "figures/tree16_resample.pdf", width = 10, height = 10)
plot(graph_fully, 
     layout = layout_fixed, 
     vertex.color= selected_colors,
     asp = 0, 
     vertex.size = 20,            # Increase node size
     vertex.label = selected_names, 
     vertex.label.cex = 1.5,      # Increase label size
     edge.width = 2*weights_fully[,3] / (N / 5)
)
dev.off()