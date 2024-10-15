#library(MASS)
#library(ggplot2)
#library(patchwork)
# Define the function
bias.test.custom <- function(result, 
                                  k_vectors = NULL, 
                                  n_perms = 1e5, 
                                  save_plots = FALSE, 
                                  output_folder = ".") {
  
  # Ensure that the required components exist in the result object
  required_components <- c("X", "Z", "vinv", "eta.hat")
  missing_components <- setdiff(required_components, names(result))
  if(length(missing_components) > 0){
    stop(paste("The result object is missing the following components:", 
               paste(missing_components, collapse = ", ")))
  }
  
  # Extract fixed effect names and count
  fixed_names <- colnames(result$X)
  n_fixed <- length(fixed_names)
  
  # Validate k_vectors if provided
  if(!is.null(k_vectors)){
    if(!is.list(k_vectors)){
      k_vectors <- list(k_vectors)
    }
    # Check each k vector
    for(i in seq_along(k_vectors)){
      k <- k_vectors[[i]]
      if(!is.numeric(k)){
        stop(paste("k_vectors[[", i, "]] is not numeric.", sep = ""))
      }
      if(length(k) != n_fixed){
        stop(paste("Length of k_vectors[[", i, "]] (", length(k), 
                   ") does not match the number of fixed effects (", n_fixed, ").", sep = ""))
      }
    }
  }
  
  # Generate k_list and corresponding names
  if(is.null(k_vectors)){
    # Generate standard basis vectors
    k_list <- lapply(1:n_fixed, function(j) {
      k <- rep(0, n_fixed)
      k[j] <- 1
      return(k)
    })
    # Assign names as fixed effect names
    names(k_list) <- fixed_names
  } else {
    # Use custom k_vectors
    k_list <- k_vectors
    # Generate descriptive names based on non-zero coefficients
    k_names <- sapply(k_list, function(k) {
      terms <- which(k != 0)
      if(length(terms) == 0){
        return("k_all_zero")
      }
      # Create strings like "1 * fixed_name1 + -1 * fixed_name3"
      term_strings <- sapply(terms, function(idx) {
        coef <- k[idx]
        # Format coefficient: avoid "+ -1", instead use "-1"
        if(coef < 0){
          return(paste0(coef, " * ", fixed_names[idx]))
        } else {
          return(paste0(coef, " * ", fixed_names[idx]))
        }
      })
      # Combine terms with " + "
      full_name <- paste(term_strings, collapse = " + ")
      return(full_name)
    })
    names(k_list) <- k_names
  }
  
  # Initialize storage vectors
  n_tests <- length(k_list)
  permutation_percentile_rec <- numeric(n_tests)
  names(permutation_percentile_rec) <- names(k_list)
  
  nu_prime_eta_rec <- numeric(n_tests)
  names(nu_prime_eta_rec) <- names(k_list)
  
  mean_nu_eta_ttest_rec <- numeric(n_tests)
  names(mean_nu_eta_ttest_rec) <- names(k_list)
  
  mean_eta_hat_rec <- numeric(n_tests)
  names(mean_eta_hat_rec) <- names(k_list)
  
  p_value_nu_eta_test_rec <- numeric(n_tests)
  names(p_value_nu_eta_test_rec) <- names(k_list)
  
  p_value_eta_hat_test_rec <- numeric(n_tests)
  names(p_value_eta_hat_test_rec) <- names(k_list)
  
  # Initialize a list to store plots
  plot_list <- list()
  
  # Create output folder if saving plots
  if(save_plots){
    if(!dir.exists(output_folder)){
      dir.create(output_folder, recursive = TRUE)
    }
  }
  
  # Loop over each k vector
  for(j in 1:n_tests){
    
    fe_name <- names(k_list)[j]
    k <- k_list[[j]]
    
    # Compute nu
    nu <- k %*% ginv(as.matrix(t(result$X) %*% result$vinv %*% result$X)) %*% 
      t(result$X) %*% result$vinv %*% result$Z
    
    # Compute observed value
    observed.value <- as.numeric(nu %*% result$eta.hat)
    
    # Initialize permutation record
    permutation.record <- numeric(n_perms)
    
    # Perform the permutation test
    for(p in 1:n_perms){
      shuffled_eta_hat <- sample(result$eta.hat)
      permutation.record[p] <- as.numeric(nu %*% shuffled_eta_hat)
    }
    
    # Compute permutation percentile
    permutation_percentile_rec[j] <- ecdf(c(permutation.record, observed.value))(observed.value)
    
    # Record nu'eta
    nu_prime_eta_rec[j] <- observed.value
    
    # Compute means
    mean_nu_eta_ttest_rec[j] <- mean(permutation.record)
    mean_eta_hat_rec[j] <- mean(result$eta.hat)
    
    # Perform t-tests
    p_value_nu_eta_test_rec[j] <- t.test(permutation.record)$p.value
    p_value_eta_hat_test_rec[j] <- t.test(as.vector(result$eta.hat))$p.value
    
    # **Plotting Section**
    permutation_record<-0
    # Create a data frame for plotting
    plot_data <- data.frame(
      permutation_record = permutation.record
    )
    
    # Define dynamic x-axis limits with a 5% margin
    combined_values <- c(permutation.record, observed.value)
    x_min <- min(combined_values) - 0.05 * abs(min(combined_values))
    x_max <- max(combined_values) + 0.05 * abs(max(combined_values))
    
    # Generate the histogram with ggplot2
    p <- ggplot(plot_data, aes(x = permutation_record)) +
      geom_histogram(binwidth = (x_max - x_min) / 100, fill = "lightblue", color = "black") +
      geom_vline(xintercept = observed.value, color = "red", linewidth = 1.2) +  # Updated 'linewidth'
      labs(
        title = paste("Test for:", fe_name),
        x = expression(nu %*% eta),
        y = "Frequency"
      ) +
      coord_cartesian(xlim = c(x_min, x_max)) +  # Use coord_cartesian instead of xlim
      theme_minimal()
    
    # Display the plot
    print(p)
    
    # Save the plot if required
    if(save_plots){
      # Create a valid filename by replacing non-alphanumeric characters with underscores
      filename <- paste0("permutation_test_", gsub("[^A-Za-z0-9]", "_", fe_name), ".png")
      ggsave(filename, plot = p, path = output_folder, width = 8, height = 6)
      cat("Saved plot for fixed effect:", fe_name, "as", file.path(output_folder, filename), "\n")
    }
    
    # Store the plot in the list
    plot_list[[j]] <- p
    
    # Progress Indicator
    cat("Completed permutation test and plot for fixed effect:", fe_name, "\n")
    
  }
  
  # **Combine all plots into a single multi-panel figure using patchwork**
  if(length(plot_list) > 0){
    # Determine grid layout based on the number of plots
    n_cols_combined <- ceiling(sqrt(n_tests))
    n_rows_combined <- ceiling(n_tests / n_cols_combined)

# Combine plots
combined_plot <- wrap_plots(plot_list, ncol = n_cols_combined, nrow = n_rows_combined) + 
  plot_annotation(title = "Permutation Test Results for All Fixed Effects")

# Display the combined plot
print(combined_plot)

# Save the combined plot if required
if(save_plots){
  ggsave("combined_permutation_tests.png", combined_plot, 
         path = output_folder, width = 4 * n_cols_combined, 
         height = 4 * n_rows_combined)
  cat("Saved combined permutation tests plot as", 
      file.path(output_folder, "combined_permutation_tests.png"), "\n")
}
  }
  
  # **Compile the results into a data frame**
  permutation_results <- data.frame(
    Fixed_Effect = names(k_list),
    Permutation_Percentile = permutation_percentile_rec,
    Nu_Prime_Eta = nu_prime_eta_rec,
    Mean_Nu_Eta = mean_nu_eta_ttest_rec,
    Mean_Eta_Hat = mean_eta_hat_rec,
    P_Value_Nu_Eta_Test = p_value_nu_eta_test_rec,
    P_Value_Eta_Hat_Test = p_value_eta_hat_test_rec,
    stringsAsFactors = FALSE
  )
  
  # Print the results data frame
  print(permutation_results)
  
  # Save the results data frame if required
  if(save_plots){
    write.csv(permutation_results, file = file.path(output_folder, "permutation_test_results.csv"), 
              row.names = FALSE)
    cat("Saved permutation test results as", 
        file.path(output_folder, "permutation_test_results.csv"), "\n")
  }
  
  # Return the results and plots
  return(list(permutation_results = permutation_results, plot_list = plot_list))
}


