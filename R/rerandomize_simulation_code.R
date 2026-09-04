
# Calibrated simulation code
#
# This file has code for running fingerprint and rerandomization
# simulations, which are simulations that target a specific dataset of
# set of estimates about a specific dataset.


# ****************************************************************************
# Fingerprint simulation code ----
# ****************************************************************************



#' Summarize set of compare_methods results across estimators
#'
#' @param rps A list of data frames, each containing the results of a
#'   compare_methods call.
#' @param summarize_results If "cross", summarize the results by
#'   calculating range statistics across the set, grouping estimators
#'   by estimand within each runID.  If "cross-agg", aggregate the
#'   sets of ranges after doing this and return summary of that.  If
#'   "method", summarize the results for each method. Otherwise just
#'   stack the initial set of results (if not already stacked)
#'
#' @return A data frame containing the summarized results.
#'
#' @export
summarize_simulation_results <- function( rps, summarize_results = "method" ) {

  if ( !is.data.frame(rps) ) {
    rps <- bind_rows( rps, .id = "runID" )
  }

  if ( !( "weight" %in% names(rps) ) ) {
    rps$weight = get_estimand( rps$method )
  }


  if ( summarize_results == "cross" || summarize_results == "cross-agg" ) {

    # Calculate range across estimators
    rr_overall <- rps %>%
      group_by( .data$runID ) %>%
      summarise( rangeATE = max(.data$ATE_hat, na.rm=TRUE) - min(.data$ATE_hat, na.rm=TRUE),
                 na_ATE = sum( is.na( .data$ATE_hat ) ),
                 ratioSE = max( .data$SE_hat, na.rm=TRUE ) / min( .data$SE_hat, na.rm=TRUE ),
                 na_SE = sum( is.na( .data$SE_hat ) ) )

    rr_estimand <- rps %>%
      # mutate( weight = forcats::fct_recode( weight,
      #                              C = "Cluster",
      #                              P = "Person" ) ) %>%
      group_by( .data$runID, .data$weight ) %>%
      summarise( rangeATE = max(.data$ATE_hat, na.rm=TRUE) - min(.data$ATE_hat, na.rm=TRUE),
                 na_ATE = sum( is.na( .data$ATE_hat ) ),
                 ratioSE = max( .data$SE_hat, na.rm=TRUE ) / min( .data$SE_hat, na.rm=TRUE ),
                 na_SE = sum( is.na( .data$SE_hat ) ),
                 .groups ="drop" )

    rr_overall$weight = "Overall"

    rr = bind_rows( rr_overall, rr_estimand )

    if ( summarize_results == "cross-agg" ) {
      ragg <- rr %>%
        group_by( .data$weight ) %>%
        rename( group = "weight" ) %>%
        summarise( mn_rngATE = mean( .data$rangeATE ),
                   sd_rngATE = sd( .data$rangeATE ),
                   prob_10p = mean( .data$rangeATE > 0.10 ),
                   mn_ratSE = mean( .data$ratioSE ),
                   prob_2x = mean( .data$ratioSE >= 2 ),
                   sd_ratSE = sd( .data$ratioSE ),
                   mn_naATE = mean( .data$na_ATE ),
                   sd_naATE = sd( .data$na_ATE ),
                   mn_naSE = mean( .data$na_SE ),
                   sd_naSE = sd( .data$na_SE ),
                   .groups = "drop" )

      stopifnot( nrow( ragg ) == 3 )

      ragg = ragg[ c( 1,3,2), ]

      ragg
    } else {
      rr
    }

  } else if ( summarize_results == "method" ) {
    # Calculate performance statistics for each method
    rps %>%
      group_by( .data$method, .data$weight ) %>%
      summarise( EATE = mean( .data$ATE_hat, na.rm=TRUE ),
                 SE = sd( .data$ATE_hat, na.rm=TRUE ),
                 p.noATE = mean( is.na( .data$ATE_hat ) ),
                 p.noSE = mean( is.na( .data$SE_hat ) ),
                 ESEhat = sqrt( mean( .data$SE_hat^2, na.rm=TRUE ) ),
                 sdSEhat = sd( .data$SE_hat, na.rm=TRUE ),
                 Q10 = quantile( .data$SE_hat, 0.1, na.rm=TRUE ),
                 Q90 = quantile( .data$SE_hat, 0.9, na.rm=TRUE ),
                 .groups = "drop" )
  } else {
    rps
  }

}


# ****************************************************************************
# Rerandomization simulation code ----
# ****************************************************************************




#' Run a finite sample simulation by permuting treatment assignment
#' labels within each block.
#'
#' Before starting simulation, put data in canonical form, imputing
#' missing values and so forth as needed. Call compare_methods
#' repeatidly on the resulting treatment-permuted data.
#'
#' @param formula A formula object specifying the model to estimate.
#' @param control_formula A formula object specifying the control
#'   variables; passed to compare_methods().
#' @param data A data frame containing the data to analyze.
#' @param R The number of simulations to run.
#' @param summarize_results If TRUE, summarize the results of the
#'   compare_methods call, FALSE return raw estimates.
#' @param parallel If TRUE, run simulations in parallel.
#' @param include_empirical If TRUE, include the empirical results as
#'   one of the "simulation replicates".  This will be included in any
#'   summarization, so be warned.
#' @param warn_missing If TRUE, warn when rows are dropped or
#'   covariates are imputed while canonicalizing/patching the data.
#' @param patch_data If TRUE, impute missing covariates via
#'   `patch_data_set()` before simulating.
#' @param ... Additional arguments to pass to compare_methods.
#'
#' @export
run_rerandomize_simulation <- function( formula,
                                        data,
                                        R = 100,
                                        control_formula = NULL,
                                        summarize_results = FALSE,
                                        include_empirical = TRUE,
                                        warn_missing = TRUE,
                                        patch_data = TRUE,
                                        parallel = FALSE, ... ) {

  # Clean up the dataset
  data = make_canonical_data( formula=formula, data=data,
                              control_formula = control_formula,
                              warn_missing = warn_missing,
                              patch_data = patch_data )
  if ( patch_data ) {
    control_formula = attr( data, "control_formula" )
  }



  if ( is.null( data$blockID ) ) {
    formula = Yobs ~ Z | clusterID
  } else {
    formula = Yobs ~ Z | clusterID | blockID

    # handle blocks with singletons
    data = patch_singleton_blocks( Yobs ~ Z | clusterID | blockID, data=data,
                                   drop_data = TRUE,
                                   warn_missing = warn_missing )

  }

  empirical = NULL
  if ( include_empirical ) {
    empirical <- compare_methods( formula = formula,
                                  control_formula=control_formula,
                                  data = data, ... )
  }

  rerandomize_and_analyze <- function( formula, data,
                                       ... ) {

    #form <- clusterRCT:::deconstruct_var_formula( formula, data = data )
    #data[form$Z] = rerandomize( data[[form$Z]], data[[form$clusterID]], data[[form$blockID]] )

    data$Z = rerandomize( data$Z, data$clusterID, data$blockID )
    compare_methods( formula = formula,
                     control_formula=control_formula,
                     data = data, ... )
  }

  if ( parallel ) {
    if ( !requireNamespace( "furrr", quietly = TRUE ) ) {
      stop( "parallel = TRUE requires the 'furrr' package. Install it with install.packages('furrr').", call. = FALSE )
    }

    if ( future::nbrOfWorkers() == 1 ) {
      warning( "You should set future::plan yourself -- using default multisession numCores - 1 plan" )

      oplan = future::plan(future::multisession, workers = parallel::detectCores() - 1 )

      on.exit(future::plan(oplan), add = TRUE)

    }
    rps = furrr::future_map( 1:R, \(.) { rerandomize_and_analyze( formula, data, ... ) },
                             .options = furrr::furrr_options( seed = TRUE )
                             )

  } else {
    rps = purrr::map( 1:R, \(.) { rerandomize_and_analyze( formula, data, ... ) },
                      .progress = TRUE )
  }

  names( rps ) = paste0( "R", 1:R )

  if ( include_empirical ) {
    rps = append(  list( empirical = empirical ), rps )
  }

  summarize_simulation_results( rps, summarize_results = summarize_results )

}






#' Return new treatment assignment based on passed one.
#'
#' Rerandomizes treatment at the cluster level (within block, if
#' `blockID` is given), preserving the original per-block (or
#' overall) proportion of clusters treated.
#'
#' @param Z Vector of original treatment assignments (1 = treated, 0 =
#'   control), one entry per row/unit.
#' @param clusterID Vector of cluster IDs, one entry per row/unit.
#' @param blockID Vector of block IDs, one entry per row/unit. If
#'   NULL, rerandomization ignores blocking.
#'
#' @return vector of 1s and 0s.
#'
#' @export
rerandomize <- function( Z, clusterID, blockID = NULL ) {

  stopifnot( length( Z ) == length( clusterID ) )
  stopifnot( is.null( blockID ) || length( Z ) == length( blockID ) )

  newZ = NA
  if ( !is.null( blockID ) ) {
    cnts <- tibble( clusterID = clusterID, blockID = blockID, Z = Z ) %>%
      group_by( blockID ) %>%
      mutate( p = mean( Z[ !duplicated( clusterID ) ] ) ) %>%
      ungroup() %>%
      mutate( Z = randomizr::block_and_cluster_ra( blocks=blockID,
                                                   clusters=clusterID,
                                                   prob_unit=p) )
    newZ = cnts$Z
  } else {
    # No blocking
    p = rep( mean( Z[ !duplicated( clusterID ) ] ), length( clusterID ) )
    newZ = randomizr::cluster_ra( clusters=clusterID,
                                  prob_unit=p )
  }

  return( newZ )

}

