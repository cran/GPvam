bias.test.custom <- function(result, 
                             k_vectors = NULL, 
                             n_perms = 1e5
                          ) {
  

  required_components <- c("X", "Z", "vinv", "eta.hat", "G", "num.teach", "persistence")
  missing_components <- setdiff(required_components, names(result))
  if(length(missing_components) > 0){
    stop(paste("The result object is missing the following components:", 
               paste(missing_components, collapse = ", ")))
  }
  

  if(!(result$persistence %in% c("CP", "VP", "ZP"))){
    warning("The function is currently not enabled for persistence types other than 'CP', 'VP', or 'ZP'.")
  }
  

  fixed_names <- colnames(result$X)
  n_fixed <- length(fixed_names)
  

  if(!is.null(k_vectors)){
    if(!is.list(k_vectors)){
      k_vectors <- list(k_vectors)
    }
  
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
  
  if(is.null(k_vectors)){
    k_list <- lapply(1:n_fixed, function(j) {
      k <- rep(0, n_fixed)
      k[j] <- 1
      return(k)
    })
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
        if(coef == 1){
          return(fixed_names[idx])
        } else if(coef == -1){
          return(paste0("- ", fixed_names[idx]))
        } else {
          return(paste0(coef, " * ", fixed_names[idx]))
        }
      })
      # Combine terms with " + "
      full_name <- paste(term_strings, collapse = " + ")
      # Clean up signs
      full_name <- gsub("\\+ -", "- ", full_name)
      return(full_name)
    })
    names(k_list) <- k_names
  }
  
  # Initialize storage vectors
  n_tests <- length(k_list)
  permutation_p_value_rec <- numeric(n_tests)
  names(permutation_p_value_rec) <- names(k_list)
  
  nu_prime_eta_rec <- numeric(n_tests)
  names(nu_prime_eta_rec) <- names(k_list)
  
  # Initialize a list to store plots
  plot_list <- list()
  
  
  # Identify indices for each variance component in G
  num_teach <- result$num.teach
  cum_teach <- c(0, cumsum(num_teach))
  n_eta <- length(result$eta.hat)
  if(cum_teach[length(cum_teach)] != n_eta){
    stop("The sum of num.teach does not match the length of eta.hat.")
  }
  
  # Permutation indices for each group
  group_indices <- lapply(1:length(num_teach), function(i){
    (cum_teach[i] + 1):cum_teach[i + 1]
  })
  sym_mat <- symmpart(t(result$X) %*% result$vinv %*% result$X)
  inv_sym_mat <- chol2inv(chol(sym_mat))
  common_part <- inv_sym_mat %*% t(result$X) %*% result$vinv %*% result$Z
  # Loop over each k vector
  for(j in 1:n_tests){
    
    fe_name <- names(k_list)[j]
    k <- k_list[[j]]
    
    # Compute nu
    nu <- k %*% common_part
    
    # Compute observed value
    observed.value <- as.numeric(nu %*% result$eta.hat)
    
    # Initialize permutation record
    permutation.record <- numeric(n_perms)
    
    # Perform the permutation test with group-wise permutations
    for(p in 1:n_perms){
      shuffled_eta_hat <- result$eta.hat
      # Permute within each group
      for(g in seq_along(group_indices)){
        idx <- group_indices[[g]]
        shuffled_eta_hat[idx] <- sample(result$eta.hat[idx])
      }
      permutation.record[p] <- as.numeric(nu %*% shuffled_eta_hat)
    }
    
    # Compute permutation p-value
    permutation_p_value_rec[j] <- mean(abs(permutation.record) >= abs(observed.value))
    
    # Record nu'eta
    nu_prime_eta_rec[j] <- observed.value
    
    # Plotting Section
    plot_data <- data.frame(
      permutation_record = permutation.record
    )
    
    # Define dynamic x-axis limits with a 5% margin
    combined_values <- c(permutation.record, observed.value)
    x_min <- min(combined_values) - 0.05 * abs(min(combined_values))
    x_max <- max(combined_values) + 0.05 * abs(max(combined_values))
    
    # Generate the histogram with ggplot2
    p <- ggplot(plot_data, aes(x = .data$permutation_record)) +  # Updated line
      geom_histogram(binwidth = (x_max - x_min) / 100, fill = "lightblue", color = "black") +
      geom_vline(xintercept = observed.value, color = "red", linewidth = 1.2) +
      labs(
        title = paste("Permutation Test for:", fe_name),
        x = expression(nu %*% hat(eta)),
        y = "Frequency"
      ) +
      coord_cartesian(xlim = c(x_min, x_max)) +
      theme_minimal()
    
    # Display the plot
    print(p)
    
    

    
    # Store the plot in the list
    plot_list[[j]] <- p
    
    # Progress Indicator
    cat("Completed permutation test and plot for fixed effect:", fe_name, "(", j, "of", n_tests, ")\n")
    
  }
  
  # Combine all plots into a single multi-panel figure using patchwork
  if(length(plot_list) > 0 && n_tests > 1){
    # Determine grid layout based on the number of plots
    n_cols_combined <- ceiling(sqrt(n_tests))
    n_rows_combined <- ceiling(n_tests / n_cols_combined)
    
    # Combine plots
    combined_plot <- wrap_plots(plot_list, ncol = n_cols_combined, nrow = n_rows_combined) + 
      plot_annotation(title = "Permutation Test Results for All Fixed Effects")
    
    # Display the combined plot
    print(combined_plot)
    

  }
  
  # Compile the results into a data frame
  permutation_results <- data.frame(
    Fixed_Effect = names(k_list),
    Nu_Prime_Eta = nu_prime_eta_rec,
    Permutation_P_Value = permutation_p_value_rec,
    stringsAsFactors = FALSE
  )
  
  # Arrange the columns
  permutation_results <- permutation_results[, c("Fixed_Effect", "Nu_Prime_Eta", "Permutation_P_Value")]
  
 
  print(permutation_results)
  

  
  # Return the results and plots
  return(list(permutation_results = permutation_results, plot_list = plot_list))
}
