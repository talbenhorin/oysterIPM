## Description
## ------------
## This script builds the host-parasite IPM model to estimate mean time to extinction,
## population growth rate, and stable distributions.

## oyster/IPM

#############################################################################

# Source the relevant functions (CHECK YOUR WORKING DIRECTORY!)
source("eviction_fiunctions.R")
source("IPM_functions_for_R.R")

# Read in the fitted parameters that were obtained by running
# ipm_vital_rate_parameter_analysis.R

# Parameters for linear growth function
params_no_26_l = readRDS("../results/IPM_parameters_no_26_linear.rds")

# Parameters for nonlinear growth function
params_no_26_nl = readRDS("../results/IPM_parameters_no_26_nonlinear.rds")

# Decision parameter. If FALSE, considers a non-linear growth curve
linear = TRUE

# Specify parameters based on temperature relationship
if(linear){
  tparams = params_no_26_l
} else{
  tparams = params_no_26_nl
}

all_temps = 4:20 # Vector of temperatures
#all_temps = c(4, 12, 20)
days = 3 # length of time between swabs

# Lower and upper bounds (chosen to minimize eviction)
min_size = -5
max_size = 18

# Number of cells in the discretized matrix. Adding more just gives more resolution
bins = 100

# Arrays to save results
stable_dist_results = array(NA, dim=c(length(all_temps), bins + 1))
sd_means = array(NA, dim=length(all_temps))
sd_vars = array(NA, dim=length(all_temps))
absorption_results = array(NA, dim=c(length(all_temps), bins + 1))
absorption_var_results = array(NA, dim=c(length(all_temps), bins + 1))
lambdas = array(NA, dim=length(all_temps))
evict_dlambdas = array(NA, dim=length(all_temps))
evict_max = array(NA, dim=length(all_temps))
elasticity_results = array(NA, dim=c(length(all_temps), bins, bins))

# For loops to perform calculations for multiple temperatures
for(i in 1:length(all_temps)){
  
  desired_temp = all_temps[i]
  
  # Specific parameters for each temperature
  params = set_temp_params(desired_temp, tparams, linear)
  
  ##############################################################################
  
  # Exploring the IPM model
  
  # In this section we analyze the IPM model by discretizing it and treating it
  # like a matrix population model. Using functions defined in
  # "IPM_functions_for_R"
  
  # Get the relevant matrix parameters: midpoints, cell edges, and cell width
  matrix_params = set_discretized_values(min_size, max_size, bins)
  bnd = matrix_params$bnd
  y = matrix_params$y
  h = matrix_params$h
  
  # Calculate and plot eviction.  Eviction is most important when e is large
  # and. Using the script provided by Williams et al. for checked for eviction
  
  k_xpx = function(xp, x, params){
    
    # Define the full kernel (for eviction calculations)
    
    return(g_xpx(xp, x, params) * s_x(x, params) * (1 - r_x(x, params)))
    
  }
  
  evict_results = evictionMeasuresFC.Iter(g_xpx, k_xpx, s_x, min_size, max_size, params)
  evict_dlambdas[i] = evict_results$dlambda
  evict_max[i] = max(evict_results$evict)
  #plot(evict_results$y, evict_results$evict)
  
  # Get the kernel of the IPM without the zero class
  kernels = get_the_kernel(g_xpx, s_x, r_x, bins, y, params, h, plot_it=T,
                           save=paste("../results/kernel_plots_temp", all_temps[i], sep=""))
  P = kernels$P
  
  ##############################################################################
  
  # Add on the uninfected class
  
  # Compute the uninfected to infected column.  The zero to zero transition is
  # given by surviving with a load of zero ($s(0)$) and not becoming infected (1
  # - $p_{inf}$). The other transitions are given by surviving the zero class
  # ($s(0)$) and becoming infected ($p_{inf}$) and gaining a given load $x'$
  # ($\phi(x')$), which needs to be calculated with the midpoint rule.
  
  # The zero -> zero transition
  zero_trans = params$class_zero_surv * (1 - params$prob_inf)
  
  # The zero to nonzero transitions: Transitioning to zero to infected
  prob_clump = h * dnorm(y, mean=params$clump_mean, sd=params$clump_sd)
  ztnz = params$class_zero_surv * params$prob_inf * prob_clump
  col1 = c(zero_trans, ztnz)
  
  
  # Compute the infected to uninfected row. This is just the probability of
  # surviving with load $x$ ($s(x)$) and then lossing an infection $r(x)$.  We
  # just need to multiply $s(x) r(x)$
  
  # The loss probabilities: Transitioning from non-zero to zero
  nztz = kernels$S * kernels$R
  row1 = c(zero_trans, nztz)
  
  # Now we just need to stick these extra columns onto the **P** matrix. The
  # most likely thing that is going to happen is you are going to stay in the
  # uninfected class.
  
  full_P = get_full_P(P, row1, col1, min_size, y, plot_it=F)
  
  ##############################################################################
  
  # Now we need to get some of the predictions from the model
  
  # Get the mean absorption times and the variance
  absorp_times = absorption_times(full_P)
  absorption_results[i, ] = absorp_times$mean
  absorption_var_results[i, ] = absorp_times$var
  
  # Get the sensitivities and elasticities of lambda on each tranistion
  # probability.
  
  # The dominant eigenvalue
  lam = Re(eigen(full_P)$values[1])
  lambdas[i] = lam
  
  # Get the elasticity and sensitivity
  es_res = get_elasticity_and_sensitivity(P, h, y,
                                          plot_it=T,
                                          save=paste("../results/elas_temp",
                                                     all_temps[i], sep=""))
  
  elasticity_results[i, , ] = es_res$elas
  
  # Stable distribution of individuals
  stable_dist = get_stable_dist(full_P)
  stable_dist_results[i, ] = stable_dist
  
  # Get mean and variance conditional on infection
  cond_stable_dist = stable_dist[2:length(stable_dist)] / (1 - stable_dist[1])
  cond_mean = sum(cond_stable_dist * y)
  cond_var = sum(cond_stable_dist * y^2) - (cond_mean)^2
  sd_means[i] = cond_mean
  sd_vars[i] = cond_var
  
} # End for loop for all temperatures