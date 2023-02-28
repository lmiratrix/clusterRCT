#
# Design-based methods follow those found in Schochet JASA paper
# "Design-Based Ratio Estimators and Central Limit Theorems for
# Clustered, Blocked RCTs" Peter Z. Schochet, Nicole E. Pashley,
# Luke W. Miratrix, and Tim Kautz
#


get_overall_ATE <- function( mod, weights ) {
    cc = coef(mod)

    J = length( weights )

    ids <- grep( "Z:", names(cc) )
    stopifnot(length(ids) == J)
    ATE_hats <- cc[ids]

    weighted.mean( ATE_hats, w = weights )
}



#' Estimate ATE via design-based approach
#'
#' This follows "Design-Based Ratio Estimators and Central Limit
#' Theorems for Clustered, Blocked RCTs" by Peter Z. Schochet, Nicole
#' E. Pashley, Luke W. Miratrix, and Tim Kautz.
#'
#' It runs a weighted linear regression, as discussed in that paper.
#' In this implementation, we run the regression on the
#' cluster-aggregated data.
#'
#' @inheritParams linear_model_estimators
#' @param aggregated TRUE means data is already aggregated (and in
#'   canonical form).  FALSE means it is not.
#'
#' @return tibble of estimates using different varieties of the
#'   methods described in the paper.
#'
#' @export
design_based_estimators <- function( formula,
                                     data = NULL,
                                     control_formula = NULL,
                                     weight = c( "individual", "cluster" ),
                                     aggregated = FALSE ) {

    require( estimatr )

    # Determine which version of the estimator we are doing.
    weight <- match.arg(weight)

    if ( !aggregated ) {
        data = make_canonical_data( formula=formula, data=data,
                                control_formula = control_formula )
        # Collapse to clusters
        data_agg = aggregate_data( data, control_formula )
        control_formula = attr( data_agg, "control_formula" )
    } else {
        data_agg = data
    }
    has_site = "siteID" %in% names( data )


    # Make our weights variable
    suff = ""
    if (weight == "individual") {
        data_agg$.weight <- data_agg$n
        suff = "_indiv"
    } else {
        data_agg$.weight <- 1 / data_agg$n
        suff = "_clust"
    }

    # site level stuff: only if we have a siteID do we fit an
    # interacted model or allow for fixed effects for each site.
    ATE_int = NA
    if ( has_site ) {
        # Fully interacted model
        form = make_regression_formula( Yobs = "Ybar", interacted = TRUE,
                                        control_formula = control_formula )
        mod = lm_robust( form, data=data_agg, weights = .weight,
                         ci = FALSE, se_type="none" )
        stwt <- data_agg %>% group_by( siteID ) %>%
            summarise( wt = sum( .weight ) )
        ATE_int = get_overall_ATE(mod, stwt$wt )

        # fixed effects model
        form = make_regression_formula( Yobs = "Ybar", FE = TRUE,
                                        control_formula = control_formula )
        mod_FE = lm_robust( form, data=data_agg, weights = .weight,
                         ci = FALSE, se_type="none" )
        ATE_FE = coef(mod_FE)[ "Z" ]
        tibble(
            method = c( glue::glue( "DB{suff} (int)" ),
                        glue::glue( "DB{suff} (FE)" ) ),
            ATE_hat = c( ATE_int, ATE_FE ),
            SE_hat = c( NA, NA ),
            p_value = c( NA, NA )
        )
    } else {
        # We have to use simple Cluster randomized estimators.
        form = make_regression_formula( Yobs = "Ybar", FE = FALSE,
                                        control_formula = control_formula )
        mod = lm_robust( form, data=data_agg, weights = .weight,
                            ci = FALSE, se_type="none" )
        ATE_FE = coef(mod)[ "Z" ]
        tibble(
            method = glue::glue( "DB{suff}" ),
            ATE_hat = c( ATE_FE ),
            SE_hat = c( NA ),
            p_value = c( NA )
        )

    }
}






