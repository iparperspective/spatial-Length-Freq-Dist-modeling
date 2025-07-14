#' Spatial LFD Simulation and Modeling
#' 
#' This script simulates spatially correlated abundance-at-length distributions and 
#' fits hierarchical Bayesian models using INLA-SPDE approach.
#' 
#' Author: Iosu Paradinas
#' Date: 09/07/2025
#' 
#' Required packages and setup ----

# Load required libraries
library(INLA)
library(inlabru)
library(tidyverse)
library(raster)
library(scales)
library(sf)
library(viridis)

# functions used from the book "Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA"
source("spde-book-functions.R")

# Helper functions ----
logit <- function(x) log(x / (1 - x))
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Configuration parameters ----
config <- list(
  # Spatial correlation parameters
  rho_size = 0.75,        # Size correlation coefficient
  range = 50,             # Spatial range
  params = c(variance = 0.5, kappa = 1),
  
  # Grid parameters
  dims = 101,             # Grid dimensions (101x101)
  n_size = 16,            # Number of size classes
  
  # Mesh parameters
  max_edge_factor = 0.05, # Max edge as fraction of dims
  bound_outer = 0.3,      # Outer boundary
  cutoff = 0.5,           # Mesh cutoff
  
  # Covariate parameters
  min_covar = 0,          # Minimum covariate value
  max_covar = 30,         # Maximum covariate value
  
  # Population parameters
  pop_age_means = c(1000, 750, 500, 100), # Mean population by age
  cv = 0.15,              # Coefficient of variation
  
  # Survey parameters
  n_survey = 400,         # Number of survey points
  n_predict = 2000,       # Number of prediction points
  n_samples = 2000        # Number of posterior samples
)

# Main simulation function ----
simulate_spatial_population <- function(config) {
  
  # Create coordinate grid
  coords <- expand.grid(1:config$dims, 1:config$dims)
  names(coords) <- c("x", "y")
  coords <- as.matrix(coords)
  n_total <- config$dims^2
  
  # Create mesh
  mesh <- inla.mesh.2d(
    loc.domain = cbind(c(0, config$dims, config$dims, 0, 0),
                      c(0, 0, config$dims, config$dims, 0)),
    max.edge = c(1, 5) * config$dims * config$max_edge_factor,
    cutoff = config$cutoff,
    offset = c(config$dims * config$max_edge_factor, -config$bound_outer)
  )
  
  # Generate Length Frequency Distribution (LFD)
  pop_age_sim <- rnorm(4, config$pop_age_means, config$pop_age_means * config$cv)
  
  LFD <- c(
    rnorm(pop_age_sim[1], 10, 1),
    rnorm(pop_age_sim[2], 15, 2),
    rnorm(pop_age_sim[3], 25, 4),
    rnorm(pop_age_sim[4], 30, 8)
  )
  
  # Discretize LFD into size classes
  lfd_seq <- seq(min(LFD), max(LFD), length.out = config$n_size)
  lfd_density <- density(LFD)
  dens <- lfd_density$y[sapply(lfd_seq, function(a, b) {
    which.min(abs(a - b))
  }, lfd_density$x)]
  
  # Generate covariate means and standard deviations
  means_cov <- 10 + cumsum(seq(0.05, 1.5, length.out = config$n_size))
  sd_cov <- round(seq(4, 8, length.out = config$n_size), 1)
  
  # Simulate correlated spatial fields
  x_k <- book.rspde(coords, 
                    range = config$range / config$params[2],
                    sigma = sqrt(config$params[1]), 
                    n = config$n_size, 
                    mesh = mesh,
                    return.attributes = FALSE)
  
  x <- x_k
  for (j in 2:config$n_size) {
    x[, j] <- config$rho_size * x[, j - 1] + 
              sqrt(1 - config$rho_size^2) * x_k[, j]
  }
  
  scaling_size <- rescale(dens, to = c(0, 1))
  
  # Generate effects and responses for each size class
  results <- list()
  
  for (i in 1:config$n_size) {
    # Generate covariate profile
    covar_profile <- round(rescale(
      5 * coordinates(coords)[, 1] + coordinates(coords)[, 1]^3,
      to = c(config$min_covar, config$max_covar)
    )) + 1
    
    # Calculate effects
    effect_covar <- rescale(
      dnorm(covar_profile, means_cov[i], sd_cov[i]),
      to = c(-1, scaling_size[i])
    )
    
    space_effect <- rescale(x[, i], to = c(-1, scaling_size[i]))
    
    # Occurrence model
    intercept_occur <- dens[i] * 5
    nu_occur <- intercept_occur + effect_covar + space_effect
    mu_occur <- inv_logit(nu_occur)
    y_occur <- rbinom(n_total, size = 1, prob = mu_occur)
    
    # Count model
    intercept_count <- dens[i] * 5
    nu_count <- intercept_count + effect_covar + space_effect
    mu_count <- exp(nu_count)
    y_count <- rpois(n_total, lambda = mu_count)
    
    # Store results
    results[[i]] <- list(
      covar_profile = covar_profile,
      effect_covar = effect_covar,
      space_effect = space_effect,
      mu_occur = mu_occur,
      y_occur = y_occur,
      mu_count = mu_count,
      y_count = y_count
    )
  }
  
  return(list(
    coords = coords,
    mesh = mesh,
    results = results,
    dens = dens,
    LFD = LFD,
    config = config
  ))
}

# Data processing function ----
process_simulation_data <- function(sim_results) {
  coords <- sim_results$coords
  results <- sim_results$results
  n_size <- sim_results$config$n_size
  n_total <- nrow(coords)
  
  # Create data frames
  sim_data <- data.frame(
    lon = rep(coordinates(coords)[, 1], n_size),
    lat = rep(coordinates(coords)[, 2], n_size),
    idx = rep(1:n_total, n_size),
    size_num = rep(1:n_size, each = n_total),
    Size = rep(paste0("Size_", 1:n_size), each = n_total)
  )
  
  # Add effects and responses
  sim_data$covar <- unlist(lapply(results, function(x) x$covar_profile))
  sim_data$CovarEffect <- unlist(lapply(results, function(x) x$effect_covar))
  sim_data$SpatialEffect <- unlist(lapply(results, function(x) x$space_effect))
  sim_data$Occurrence <- unlist(lapply(results, function(x) x$y_occur))
  sim_data$Count <- unlist(lapply(results, function(x) x$y_count))
  sim_data$`Mean Occur` <- unlist(lapply(results, function(x) x$mu_occur))
  sim_data$`Mean Abundance` <- unlist(lapply(results, function(x) x$mu_count))
  sim_data$Count_occur <- ifelse(sim_data$Count > 0, 1, 0)
  
  return(sim_data)
}

# Visualization functions ----
plot_covariate_effects <- function(sim_data, filename = NULL) {
  p <- sim_data %>%
    group_by(size_num) %>%
    slice(1) %>%
    ggplot(aes(x = covar, y = CovarEffect, color = factor(size_num))) +
    geom_line(linewidth = 1) +
    theme_bw() +
    xlab("Covariate") +
    ylab("Effect") +
    guides(color = guide_legend(title = "Size"))
  
  if (!is.null(filename)) ggsave(filename, plot = p)
  return(p)
}

plot_covariate_space <- function(sim_data, filename = NULL) {
  p <- sim_data %>%
    filter(size_num == 1) %>%
    ggplot(aes(x = lon, y = lat, fill = covar)) +
    geom_tile() +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    scale_fill_viridis() +
    guides(fill = guide_legend(title = "Covariate"))
  
  if (!is.null(filename)) ggsave(filename, plot = p)
  return(p)
}

plot_spatial_effects <- function(sim_data, filename = NULL) {
  p <- sim_data %>%
    ggplot(aes(x = lon, y = lat, fill = SpatialEffect)) +
    geom_tile() +
    facet_wrap(~size_num) +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    scale_fill_viridis() +
    guides(fill = guide_legend(title = "Spatial Effect"))
  
  if (!is.null(filename)) ggsave(filename, plot = p)
  return(p)
}

plot_count_data <- function(sim_data, survey_data, filename = NULL) {
  p <- sim_data %>%
    ggplot(aes(x = lon, y = lat, fill = Count)) +
    geom_tile() +
    facet_wrap(~size_num) +
    geom_point(data = survey_data, aes(x = lon, y = lat), 
               color = "red", inherit.aes = FALSE) +
    scale_fill_viridis() +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    guides(fill = guide_legend(title = "Count"))
  
  if (!is.null(filename)) ggsave(filename, plot = p)
  return(p)
}

# Model fitting function ----
fit_length_SDM <- function(survey_sf, mesh, covar_seq, n_size) {
  # Create covariate mapper
  mapper_covar <- bru_mapper(
    INLA::inla.mesh.1d(covar_seq, boundary = "free"),
    indexed = FALSE
  )
  
  # Define Matern prior
  matern <- inla.spde2.pcmatern(
    mesh,
    prior.range = c(20, 0.1),
    prior.sigma = c(0.5, 0.5)
  )
  
  # Define model components
  components <- Count ~ Intercept(1) +
    covar(covar, model = "rw2", mapper = mapper_covar, scale.model = TRUE,
          hyper = list(prec = list(prior = "pc.prec", param = c(4, 0.01))),
          group = size_num,
          ngroup = n_size,
          control.group = list(model = "ar1")) +
    size(Size, model = "factor_contrast") +
    S_size(geometry, model = matern,
           group = size_num,
           ngroup = n_size,
           control.group = list(model = "ar1"))
  
  # Fit model
  fit <- bru(
    components = components,
    data = survey_sf,
    family = "nbinomial",
    options = list(
      verbose = FALSE,
      control.inla = list(int.strategy = "eb")
    )
  )
  
  return(fit)
}

# Prediction function ----
make_predictions <- function(fit, sim_data_sf, n_samples = 2000, n_predict = 2000) {
  # Sample prediction locations
  prediction_data <- sim_data_sf %>%
    filter(idx %in% sample(unique(sim_data_sf$idx), n_predict))
  
  # Make predictions
  predictions <- predict(
    fit,
    newdata = prediction_data,
    formula = ~exp(Intercept + covar + S_size + size),
    n.samples = n_samples,
    include = c("Intercept", "covar", "S_size", "size")
  )
  
  # Add coordinates
  predictions$lon <- st_coordinates(predictions)[, 1]
  predictions$lat <- st_coordinates(predictions)[, 2]
  
  return(predictions)
}


  # 1. Simulate data
  sim_results <- simulate_spatial_population(config)
  
  # 2. Process data
  sim_data <- process_simulation_data(sim_results)
  
  # 3. Create survey data
  survey_indices <- sample(unique(sim_data$idx), config$n_survey)
  survey_data <- sim_data %>% filter(idx %in% survey_indices)
  survey_sf <- st_as_sf(survey_data, coords = c("lon", "lat"))
  
  # 4. Fit model
  covar_seq <- seq(min(sim_data$covar), max(sim_data$covar), length.out = config$n_size)
  fit <- fit_length_SDM(survey_sf, sim_results$mesh, covar_seq, config$n_size)
  
  # 5. Make predictions
  sim_data_sf <- st_as_sf(sim_data, coords = c("lon", "lat"))
  predictions <- make_predictions(fit, sim_data_sf, config$n_samples, config$n_predict)
  
  # 6. Generate plots
  plot_covariate_effects(sim_data, paste0(output_dir, "covariate_effects.png"))
  plot_covariate_space(sim_data, paste0(output_dir, "covariate_space.png"))
  plot_spatial_effects(sim_data, paste0(output_dir, "spatial_effects.png"))
  plot_count_data(sim_data, survey_data, paste0(output_dir, "count_data.png"))
  



  # Access results
  sim_data <- results$simulation_data
  model_fit <- results$model_fit
  predictions <- results$predictions
  
  # Print model summary
  summary(model_fit)
  
  # Quick visualization of predictions vs true values
  predictions %>%
    ggplot(aes(x = median, y = `Mean Abundance`)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~size_num) +
    theme_bw() +
    labs(title = "Predicted vs Simulated Abundances",
         x = "Predicted", y = "Simulated")
