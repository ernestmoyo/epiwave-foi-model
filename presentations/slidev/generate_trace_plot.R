# Generate MCMC trace plot from cached results
# Run from epiwave-foi-model root:
#   Rscript presentations/slidev/generate_trace_plot.R

library(ggplot2)

output_path <- "presentations/slidev/public/images/mcmc_trace_plots.png"

tryCatch({
  res <- readRDS("cache/mcmc_results.rds")
  cat("Loaded MCMC cache. Structure:\n")
  cat("  Class:", paste(class(res), collapse=", "), "\n")
  cat("  Names:", paste(names(res), collapse=", "), "\n")

  # Try to extract draws — greta mcmc returns mcmc.list (coda)
  # Common structures: res$draws_with, res$draws_without, or res itself is mcmc.list

  extract_draws <- function(obj, label) {
    if (inherits(obj, "mcmc.list")) {
      # coda mcmc.list — each element is an mcmc chain
      dfs <- lapply(seq_along(obj), function(i) {
        chain_df <- as.data.frame(obj[[i]])
        chain_df$chain <- paste("Chain", i)
        chain_df$iter  <- seq_len(nrow(chain_df))
        chain_df
      })
      df <- do.call(rbind, dfs)
      df$model <- label
      return(df)
    } else if (is.data.frame(obj) || is.matrix(obj)) {
      df <- as.data.frame(obj)
      df$chain <- "Chain 1"
      df$iter  <- seq_len(nrow(df))
      df$model <- label
      return(df)
    }
    return(NULL)
  }

  # Try common naming patterns
  draws_with <- NULL
  draws_without <- NULL

  if ("draws_with" %in% names(res)) {
    draws_with <- extract_draws(res$draws_with, "WITH offset")
  } else if ("with" %in% names(res)) {
    draws_with <- extract_draws(res[["with"]], "WITH offset")
  }

  if ("draws_without" %in% names(res)) {
    draws_without <- extract_draws(res$draws_without, "WITHOUT offset")
  } else if ("without" %in% names(res)) {
    draws_without <- extract_draws(res[["without"]], "WITHOUT offset")
  }

  # If res itself is mcmc.list (single model run)
  if (is.null(draws_with) && inherits(res, "mcmc.list")) {
    draws_with <- extract_draws(res, "WITH offset")
  }

  # Combine whatever we found
  all_draws <- rbind(draws_with, draws_without)

  if (!is.null(all_draws) && nrow(all_draws) > 0) {
    # Identify parameter columns (exclude metadata)
    param_cols <- setdiff(names(all_draws), c("chain", "iter", "model"))
    cat("  Parameters found:", paste(param_cols, collapse=", "), "\n")

    # Pivot longer
    plot_df <- tidyr::pivot_longer(all_draws,
                                    cols = all_of(param_cols),
                                    names_to = "parameter",
                                    values_to = "value")

    # Compute R-hat and ESS per model+parameter
    library(coda)
    add_diagnostics <- function(df) {
      labels <- c()
      for (mod in unique(df$model)) {
        for (par in unique(df$parameter)) {
          sub <- df[df$model == mod & df$parameter == par, ]
          chains <- split(sub$value, sub$chain)
          if (length(chains) >= 2) {
            mcmc_list <- mcmc.list(lapply(chains, mcmc))
            rhat <- tryCatch(gelman.diag(mcmc_list)$psrf[1,1], error = function(e) NA)
            ess  <- tryCatch(effectiveSize(mcmc_list)[1], error = function(e) NA)
          } else {
            rhat <- NA
            ess  <- tryCatch(effectiveSize(mcmc(chains[[1]])), error = function(e) NA)
          }
          rhat_str <- if (!is.na(rhat)) sprintf("R-hat=%.3f", rhat) else ""
          ess_str  <- if (!is.na(ess)) sprintf("ESS=%d", round(ess)) else ""
          flag <- if (!is.na(rhat) && rhat > 1.1) " \u26A0" else " \u2713"
          label <- paste0(mod, " \u2014 ", par, "\n", rhat_str, ", ", ess_str, flag)
          labels <- c(labels, setNames(label, paste(mod, par, sep="__")))
        }
      }
      labels
    }

    diag_labels <- add_diagnostics(plot_df)
    plot_df$facet <- diag_labels[paste(plot_df$model, plot_df$parameter, sep="__")]

    p3 <- ggplot(plot_df, aes(x = iter, y = value, colour = chain)) +
      geom_line(alpha = 0.6, linewidth = 0.3) +
      facet_wrap(~ facet, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = c("Chain 1" = "#F8766D", "Chain 2" = "#00BFC4",
                                     "Chain 3" = "#7CAE00", "Chain 4" = "#C77CFF")) +
      labs(title = "EpiWave MCMC Trace Plots",
           subtitle = "2000 samples, 1000 warmup, 2 chains",
           x = "iter", y = NULL, colour = "chain") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom",
            strip.text = element_text(size = 9),
            panel.grid.minor = element_blank())

    ggsave(output_path, p3, width = 10, height = 8, dpi = 200, bg = "white")
    cat("Trace plot saved to:", output_path, "\n")
  } else {
    cat("Could not extract draws from cache. Names found:", paste(names(res), collapse=", "), "\n")
    cat("Please copy your trace plot screenshot manually to:", output_path, "\n")
  }

}, error = function(e) {
  cat("Error reading MCMC cache:", conditionMessage(e), "\n")
  cat("This is expected if the cache contains greta/TF objects.\n")
  cat("Please copy your trace plot screenshot manually to:\n")
  cat("  ", output_path, "\n")
})
