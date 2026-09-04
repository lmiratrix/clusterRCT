
# MRStdCRT estimator wrapper




#' Estimate ATEs for a cluster RCT using the MRStdCRT pacakge
#'
#' NOT YET FUNCTIONAL -- future work.  Not exported from the package;
#' kept here as a skeleton to finish later.  See FUTURE_WORK.md.
#'
#' NOTE: You will have to manually install this package via GitHub:
#' Use:
#'
#' devtools::install_github("deckardt98/MRStdCRT")
#'
#' @param formula Formula for outcome and treatment and nesting.  If
#'   NULL, data is assumed to be in canonical form (see vignette for
#'   further discussion).
#'
#' @noRd
MRStdCRT_estimator <- function( formula,
                            data = NULL,
                            control_formula = NULL,
                            weight = c( "Person", "Cluster" ) ) {

    warning( "This method does not yet work due to difficulties mapping to the call" )

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }

    probs <- make_block_table( data ) %>%
        dplyr::select( blockID, p.tx )

    data <- left_join( data, probs, by = c("blockID" ) ) %>%
        mutate( assigned_value = ifelse( Z, p.tx, 1 - p.tx ) ) %>%
        group_by( blockID, clusterID ) %>%
        mutate( n = n() ) %>%
        ungroup()

    data

    form = Yobs ~ cluster(n)
    if ( !is.null( control_formula ) ) {
        form = update( control_formula, Yobs ~ . + cluster(n) )
    }

    data$clusterID = as.numeric( data$clusterID )
    data <- arrange( data, clusterID )
    res <- MRStdCRT_fit(
        formula = form,
        data = as.data.frame( data ),
        clusterID = "clusterID",
        trt = "Z",
        trtprob = data$assigned_value,
        method = "GEE",
        corstr = "independence",
        family = gaussian(link = "identity"),
        scale = "RD"
    )

    # TODO: mapping MRStdCRT_fit()'s return value (`res`) to this
    # package's canonical tibble( method, ATE_hat, weight, SE_hat,
    # p_value, df ) output format is the piece that isn't done yet.
    res
}




