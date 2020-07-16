## Mark Wilbur analyzing IPMs in R

library(MASS)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)
library(fields)
library(popbio)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

format_for_heat = function(mat, row_nm, col_nm, var_names){
  
  # Function converts a 2D matrix for plotting in ggplot
  
  tmat = mat
  colnames(tmat) = col_nm
  rownames(tmat) = row_nm
  long_mat = melt(tmat)
  colnames(long_mat) = var_names
  
  return(long_mat)
  
}

set_discretized_values = function(min_size, max_size, bins){
  
  # Calculates the necessary parameters to use the midpoint rule to evaluate
  # the IPM model
  
  # Parameters
  # ----------
  # min_size : The lower bound of the integral
  # max_size : The upper bound of the integral
  # bins : The number of bins in the discretized matrix
  
  # Returns
  # -------
  # list
  # min_size, max_size, bins, bnd (edges of discretized kernel), y (midpoints),
  # h (width of cells)
  
  
  # Set the edges of the discretized kernel
  bnd = min_size+c(0:bins)*(max_size-min_size) / bins
  
  # Set the midpoints of the discretizing kernel. Using midpoint rule for evaluation
  y = 0.5 * (bnd[1:bins] + bnd[2:(bins + 1)])
  
  # Width of cells
  h = y[2] - y[1]
  
  return(list(min_size=min_size, max_size=max_size, bins=bins, bnd=bnd, y=y,
              h=h))
  
}

get_full_P = function(P, row1, col1, min_size, y, plot_it=T){
  
  # Function takes in the fungal kernel with a zero class and sticks
  # the zero row and column onto it
  
  # Parameters
  # ----------
  # P : the transition matrix without a zero class
  # row1 : vector specifying the probability of any stage transition to 0
  # col1 : vector specifying the of zero transitioning to any stage
  # min_size : Lower bound of the kernels
  # y : midpoints of each stage, used for plotting
  
  # Returns
  # -------
  # : matrix
  #   The transition matrix with a 0 class
  
  full_P = matrix(NA, nrow=dim(P)[1] + 1, ncol=dim(P)[2] + 1)
  full_P[1, ] = row1
  full_P[ , 1] = col1
  full_P[2:dim(full_P)[1], 2:dim(full_P)[2]] = P
  
  # if(plot_it){
  
  #   par(mfrow=c(1, 1))
  
  #   image.plot(c(min_size - 1, y), c(min_size - 1, y), t(full_P),
  #     main="Full transition matrix", xlab="Load at t", ylab="Load at t + 1")
  
  #   abline(0,  1, lwd=3)
  
  # }
  
  return(full_P)
  
}

absorption_times = function(P){
  
  # Calculate the mean time to absoprtion (death) given a transition matrix P
  # Calculate the variance as well
  
  # Parameters
  # ----------
  # P : transition matrix
  #
  # Returns
  # -------
  # list: mean=mean absorption time, var=var absorption time
  
  I = diag(dim(P)[1])
  
  # Calculate the fundamental matrix
  N = solve((I - P))
  
  # Exepected time to absorption starting in transient state x
  mean_absorb_time = colSums(N)
  var_absorb_time = colSums((2*(N %*% N) - N)) - mean_absorb_time^2
  
  return(list(mean=mean_absorb_time, var=var_absorb_time))
  
}

get_elasticity_and_sensitivity = function(P, h, y, plot_it=T, save="temp"){
  
  # Calculate the elasticities and sensitivities of the population growth rate
  # to transition probabilities in the matrix P
  
  # Parameters
  # ----------
  # P : transition matrix
  # h : cell width
  # y : midpoints of cells, more plotting
  
  # Returns
  # -------
  # list: elas= elasticity matrix, sens= sensitivity matrix
  
  
  lam=Re(eigen(P)$values[1])
  w.eigen = Re(eigen(P)$vectors[,1])
  stable.dist = w.eigen / sum(w.eigen)
  v.eigen = Re(eigen(t(P))$vectors[,1])
  repro.val = v.eigen / v.eigen[1]
  
  v.dot.w = sum(stable.dist*repro.val) * h
  
  sens = sensitivity(P)# outer(repro.val,stable.dist) / v.dot.w
  elas = elasticity(P)# matrix(as.vector(sens)*as.vector(P) / lam, nrow=length(y))
  
  if(plot_it){
    
    myPalette = colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
    
    long_elas = format_for_heat(t(elas), y, y, c('size', 'sizeNext', 'value'))
    long_elas$id = "Elasticity"
    
    # Just plotting elasticity
    zplot = ggplot(long_elas, aes(x=size, y=sizeNext, fill=value)) + geom_tile()
    zplot = zplot + scale_fill_gradientn(colours = myPalette(100), name="Elasticity")
    zplot = zplot + theme_bw() + xlab("Log zoospore load at t") + ylab("Log zoospore load at t + 1")
    zplot = zplot + facet_wrap(~id)
    
    ggsave(paste(save, ".pdf", sep=""), height=5, width=7)
    
  }
  
  return(list(elas=elas, sens=sens))
  
}

get_stable_dist = function(P){
  
  # Get the stable distribution of a P matrix based on Caswell 2001
  
  w.eigen = Re(eigen(P)$vectors[,1])
  stable_dist = w.eigen / sum(w.eigen)
  return(stable_dist)
  
}

get_vm_ratio = function(P, y){
  
  # Get the variance to mean ratio of the stable distribution given the full
  # transition matrix P. y is the step size
  
  stable_dist = get_stable_dist(P)
  
  cond_stable_dist = stable_dist[2:length(stable_dist)] / (1 - stable_dist[1])
  cond_mean = sum(cond_stable_dist * y)
  cond_var = sum(cond_stable_dist * y^2) - (cond_mean)^2
  
  return(cond_var / cond_mean)
  
}

