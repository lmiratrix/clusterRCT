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
        datagg = aggregate_data( data, control_formula )
    } else {
        datagg = data
    }
    has_site = "siteID" %in% names( data )


    # Make our weights variable
    suff = ""
    if (weight == "individual") {
        datagg$.weight <- datagg$n
        suff = "_indiv"
    } else {
        datagg$.weight <- 1 / datagg$n
        suff = "_clust"
    }

    # site level stuff: only if we have a siteID do we fit an
    # interacted model or allow for fixed effects for each site.
    ATE_int = NA
    if ( has_site ) {
        # Fully interacted model
        form = make_regression_formula( Yobs = "Ybar", interacted = TRUE )
        mod = lm_robust( form, data=datagg, weights = .weight,
                         ci = FALSE, se_type="none" )
        stwt <- datagg %>% group_by( siteID ) %>%
            summarise( wt = sum( .weight ) )
        ATE_int = get_overall_ATE(mod, stwt$wt )

        # fixed effects model
        form = make_regression_formula( Yobs = "Ybar", FE = TRUE )
        mod_FE = lm_robust( form, data=datagg, weights = .weight,
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
        form = make_regression_formula( Yobs = "Ybar", FE = FALSE )
        mod = lm_robust( form, data=datagg, weights = .weight,
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






