require(igraph)
require(graphicalExtremes)
library(dplyr)
library(ggplot2)
library(tidyr)
library(latex2exp) 
library(gridExtra)


figure_dest_folder = "figures/"


#################
### Functions ###
#################


## Function for simulating n jumps from (thresholded) Levy measure 
## n: number of samples
## epsilon: threshold above which is simulated from the exponent measure
## Gamma: variogram matrix for the HR dependence 
## alpha: vector of marginal tail indices 
## c: 2xd matrix giving the marginal coefficients for each variable and positive and negative side (first row for positive, second row for negative half lines)
sim_jump <- function(n, epsilon, Gamma, alpha = 1, c = 1, proba_quad = NULL, tree = NULL) {
  d <- ncol(Gamma)  
  if(length(alpha)==1) alpha <- rep(alpha, times = d)
  if(!is.matrix(c)) c <- matrix(c, nrow = 2, ncol = d)

  jumps_std <- epsilon * rmpareto(n, model = "HR", d=d, par = Gamma)

  if(!is.null(proba_quad)){
    signs <- t(sapply(1:n, FUN = function(i) assign_quadrant(tree=tree, probabilities = proba_quad)))
  }
  else{
    signs <- matrix(sample(c(-1,1),size=n*d, replace=TRUE), nrow=n, ncol=d)
  }
  jumps_std <- signs * jumps_std 
 
  jumps <- sapply(1:d, function(i)  signs[,i] * (c[1 + as.numeric(signs[,i] == -1),i] * abs(jumps_std[,i]))^(1/alpha[i]))
  return(jumps)
}



##################
### Simulation ###
##################

## Simulation of the compound Poisson process.
## n is number of subintervals, rate is that of Poisson process (i.e. avg. no. of jumps)
sim_CPP <- function(n, rate, epsilon, Gamma, alpha = 1, c = 1, proba_quad = NULL, tree = NULL, emp_incs = NULL) { 
  d <- ncol(Gamma)
  ## Number of jumps in subintervals
  N <- rpois(n,rate/n) 
  ## Total number of jumps
  M <- sum(N)
  J <- sim_jump(n = M, epsilon = epsilon, Gamma = Gamma, alpha = alpha, c = c, proba_quad = proba_quad, tree = tree)
  cJ <- apply(J,2,cumsum)
  cN <- cumsum(N)
  Xtmp <- cJ[cN,]

  # Filling up with zeroes the process before any jump has occurred
  X <- rbind(matrix(rep(0,ncol(Xtmp)*(n-nrow(Xtmp))),ncol=ncol(Xtmp)), Xtmp)

  # Empirical transformation of marginals
  if(!is.null(emp_incs)){
    incs_X <- diff(rbind(rep(0, times=d),X))
    incs_trans <- matrix(NA, nrow=nrow(incs_X), ncol=d)
    for(j in 1:d){
      trans_tmp <- prepare_empirical_functions(data = emp_incs[,j])
      incs_trans[,j] <- unlist(trans_tmp$inverse_cdf(rank(incs_X[,j]) / (nrow(incs_X)+1)))
    }
    X_emp <- apply(incs_trans,2,cumsum)
    return(X_emp)
  }
  return(X)
}

## Approximate simulation of the Lévy process
## n is number of subintervals
sim_X <- function(n, rate, epsilon, Gamma, alpha = 1, c = 1, proba_quad = NULL, tree = NULL, emp_incs = NULL) {
  ## Linear drift
  # d <- ncol(Gamma)
  # drift <- matrix(0,n,d)
  ## Brownian motion
  # U <- matrix(0,n,d) #epsilon^(1-alpha/2)*t(chol(Sigma))%*%matrix(rnorm(n*d,sd=1/sqrt(n)),nrow=d)    CHANGE BACK!!!
  # BM <- t(apply(U,1,cumsum))
  ## Adding the drift, BM and CPP
  X <- sim_CPP(n = n, rate = rate, epsilon = epsilon, Gamma = Gamma, alpha = alpha, c = c, tree = tree, proba_quad = proba_quad, emp_incs = emp_incs)

  return(X)
}


## Function to plots pairs of stocks
plot_stocks <- function(data, stock1, stock2, save_plot = FALSE, suffix="") {
  df <- data.frame(
    x = data[, stock1], 
    y = data[, stock2] 
  )

  # Create the scatter plot
  gg <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = "black", alpha=.25) + 
    labs(x = stock1, y = stock2) + 
    theme_minimal() + 
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, size = 1), 
      aspect.ratio = 1 
    )
  if (save_plot) {
    save_myplot(
      plt = gg,
      plt_nm = paste0(figure_dest_folder, "scatter_", stock1, "_", stock2, "_", suffix,".pdf"),
      width = 5,
      height = 5,
      cairo = FALSE
    )
  } else {
    gg
  }
}





prepare_empirical_functions <- function(data) {
  n <- length(data)
  cdf_func <- function(x) {   
      return(ecdf(data)(x))
  }
  quantile_func <- function(q) {
      central_quantile_index <- min(max(1,floor(q * (n + 1))),n)
      return(sort(data)[central_quantile_index])
    }  
  return(list(cdf = Vectorize(cdf_func), inverse_cdf = Vectorize(quantile_func)))
}







################################
### Functions for estimation ###
################################

emp_chi_levy <- function(data, p = NULL){
    if (!is.matrix(data)) {
        stop("The data should be a matrix")
    }
    if (ncol(data) <= 1) {
        stop("The data should be a matrix with at least two columns.")
    }

    n <- nrow(data)
    d <- ncol(data)

    chi <- matrix(0, d, d)

  for (s1 in c(-1,1)) {
    for (s2 in c(-1,1)) {

    data_unif1 <- matrix(apply(s1 * data, 2, graphicalExtremes:::unif), nrow(data), 
        ncol(data))
    data_unif2 <- matrix(apply(s2 * data, 2, graphicalExtremes:::unif), nrow(data), 
        ncol(data))

    ind1 <- data_unif1 > p
    ind2 <- data_unif2 > p
    
    ind_mat <- matrix(colSums(ind1), byrow = TRUE, ncol = d, nrow = d)
    chi_tmp <- crossprod(ind1, ind2)/(1/2 * (ind_mat + t(ind_mat)))

    chi <- chi + chi_tmp
    }
  }
  return(chi/2)
}

emp_quad_prop <- function(data, tree, p){
    n <- nrow(data)
    d <- ncol(data)

    proba_pos <- proba_neg <- matrix(0, d, d)
  
    
    for (s1 in c(-1,1)) {
    for (s2 in c(-1,1)) {

    data_unif1 <- matrix(apply(s1 * data, 2, graphicalExtremes:::unif), nrow(data), 
        ncol(data))
    data_unif2 <- matrix(apply(s2 * data, 2, graphicalExtremes:::unif), nrow(data), 
        ncol(data))

    ind1 <- data_unif1 > p
    ind2 <- data_unif2 > p
    
    ind_mat <- matrix(colSums(ind1), byrow = TRUE, ncol = d, nrow = d)
    chi_tmp <- crossprod(ind1, ind2)/(1/2 * (ind_mat + t(ind_mat)))

    chi_tmp[1,2]

    if(s1*s2>0) proba_pos <- proba_pos + chi_tmp
    if(s1*s2<0) proba_neg <- proba_neg + chi_tmp
    }
  }
  return(proba_pos/(proba_neg + proba_pos))
}


assign_quadrant <- function(tree, probabilities) {
  # Get the number of nodes
  d <- vcount(tree)
  
  # Initialize a vector to store the values assigned to each node
  node_values <- numeric(d)
  
  # Assign the first node a value of +1 or -1 with probability 1/2
  node_values[1] <- sample(c(-1, 1), 1)
  
  # Use a queue to process nodes iteratively
  queue <- 1
  
  while (length(queue) > 0) {
    node <- queue[1]
    queue <- queue[-1]
    
    # Get the neighbors of the current node
    neighbors <- neighbors(tree, node, mode = "out")
    
    for (neighbor in neighbors) {
      # Only process the neighbor if it hasn't been assigned a value yet
      if (node_values[neighbor] == 0) {
        parent_value <- node_values[node]
        p_ij <- probabilities[node, neighbor]
        
        # Assign value to the neighbor based on the parent's value
        if (parent_value == -1) {
          node_values[neighbor] <- sample(c(-1, 1), 1, prob = c(p_ij, 1 - p_ij))
        } else if (parent_value == 1) {
          node_values[neighbor] <- sample(c(1, -1), 1, prob = c(p_ij, 1 - p_ij))
        }
        
        # Add the neighbor to the queue for further processing
        queue <- c(queue, neighbor)
      }
    }
  }
  
  return(node_values)
}




## Function for learning the tree
levy_mst <- function(incs, p) {
  d <- ncol(incs)
  chi <- emp_chi_levy(incs, p=p) 
  weight_matrix <- -chi

  graph.full <- igraph::make_full_graph(d)
  mst <- igraph::mst(graph = graph.full, weights = weight_matrix[igraph::ends(graph.full, 
        igraph::E(graph.full))], algorithm = "prim")
  mst
}




my_palette <- list(
  "red" = "#D55E00",
  "blue" = "#0072B2", 
  "green" = "#009E73",
  "yellow" = "#E69F00",
  "pink" = "#CC79A7",
  "light_blue" = "#56B4E9",
  "grey" = "#999999",
  "background" = "#332288"
)


save_myplot <- function(plt, plt_nm,
                        width, height, 
                        width_pdf = 50, height_pdf = 50,
                        crop = TRUE, cairo = TRUE) {
  
  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }
  
  if (cairo) {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"), 
           device = cairo_pdf, family = "Arial")
  } else {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"))
  }
  
  if (crop){
    knitr::plot_crop(plt_nm)
  } 
}


