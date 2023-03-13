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



# Notes on Shochet paper
#
# m_b^1 = Total number of clusters assigned to tx group in block b




# Implement formula (9) and following equations from Schochet et
# al paper.
#
# See page 2139-40 (section 4.2)
#
# @adt Aggregated (cluster-level) data.
# @param v Degrees of freedom.  Use v=# covariates.
schochet_variance_formula <- function( adt, v ) {

    stopifnot( !is.null( adt$resid ) )
    stopifnot( !is.null( adt$Z ) )

    calc_s2_partial <- function( adt ) {

        w_bar = mean( adt$.weight )
        tot = sum( ( (adt$.weight^2) / (w_bar^2) ) * adt$resid^2 )
        #mb = nrow( adt )

        #tot / (mb - v * p * q - 1)
        tot
    }

    adtw <- adt %>%
        group_by( siteID, Z ) %>%
        nest() %>%
        ungroup() %>%
        mutate( m = map_dbl( data, nrow ),
                totwt = map_dbl( data, ~ sum( .$.weight ) ),
                s2 = map_dbl( data, calc_s2_partial ) ) %>%
        group_by( siteID ) %>%
        mutate( p = m / sum(m) ) %>%
        ungroup() %>%
        mutate( q = totwt / sum(totwt) ) %>%
        group_by( Z ) %>%
        mutate( wt = totwt / sum(totwt),
                z = m - v * p * q - 1,
                s2 = s2 / z ) %>%
        ungroup()

    pt2 <- adtw %>%
        dplyr::select( -p, -q, -data ) %>%
        pivot_wider( names_from = "Z", values_from = c( totwt, s2, m, wt) ) %>%
        mutate( varD = s2_1 / m_1 + s2_0 / m_0,
                totwt = totwt_0 + totwt_1 )

    h = length( unique( adt$siteID ) )
    pt2$totwt = pt2$totwt / sum(pt2$totwt)

    SE_ATE <- weighted.mean( pt2$varD, w = pt2$totwt^2) / ( (h*mean(pt2$totwt) )^2 )

    list( SE_hat = sqrt( SE_ATE ),
          df = nrow(adt) - 2*h - 1 )
}




# From Schochet et al equation (15) for the restricted model with no
# interaction term for tx.
schochet_FE_variance_formula <- function( adt, v ) {

    stopifnot( !is.null( adt$resid ) )
    stopifnot( !is.null( adt$Z ) )

    adt <- adt %>%
        group_by( siteID ) %>%
        mutate( p = mean(Z),
                Ztilde = Z - p,
                m = n(),
                wbar = mean(.weight) ) %>%
        ungroup() %>%
        mutate( numers = .weight^2 * Ztilde^2 * resid^2,
                denoms = (m * p * (1-p) * wbar)^2 ) %>%
        ungroup()

    m = nrow(adt)
    h = length( unique(adt$siteID) )

    SE_ATE = sqrt( (m / (m - h - v - 1)) * sum( adt$numers ) / sum( adt$denoms ) )
    df = nrow(adt) - v - h - 1

    list( SE_hat = SE_ATE, df = df )
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

    v = number_controls(control_formula)

    # site level stuff: only if we have a siteID do we fit an
    # interacted model or allow for fixed effects for each site.
    ATE_int = NA
    if ( has_site ) {

        # Fully interacted model
        form = make_regression_formula( Yobs = "Ybar", interacted = TRUE,
                                        control_formula = control_formula )
        mod = lm( form, data=data_agg, weights = .weight )
        stwt <- data_agg %>% group_by( siteID ) %>%
            summarise( wt = sum( .weight ) )
        ATE_int = get_overall_ATE( mod, stwt$wt )

        data_agg$resid = residuals( mod )
        SE_int = schochet_variance_formula(data_agg, v = v )

        # fixed effects model
        form = make_regression_formula( Yobs = "Ybar", FE = TRUE,
                                        control_formula = control_formula )
        mod_FE = lm( form, data=data_agg, weights = .weight )
        ATE_FE = coef(mod_FE)[ "Z" ]
        data_agg$resid = residuals(mod_FE)
        SE_FE = schochet_FE_variance_formula( data_agg, v = v )

        tibble(
            method = c( glue::glue( "DB{suff} (int)" ),
                        glue::glue( "DB{suff} (FE)" ) ),
            ATE_hat = c( ATE_int, ATE_FE ),
            SE_hat = c( SE_int$SE_hat, SE_FE$SE_hat ),
            df = c( SE_int$df, SE_FE$df ),
            p_value = two_sided_p( ATE_hat, SE_hat, df )
        )
    } else {
        # We have to use simple Cluster randomized estimators.
        form = make_regression_formula( Yobs = "Ybar", FE = FALSE,
                                        control_formula = control_formula )
        mod = lm( form, data=data_agg, weights = .weight )
        ATE_FE = coef(mod)[ "Z" ]

        data_agg$siteID = "Ind"
        data_agg$resid = residuals( mod )
        SE = schochet_variance_formula(data_agg, v = v )

        tibble(
            method = glue::glue( "DB{suff}" ),
            ATE_hat = c( ATE_FE ),
            SE_hat = c( SE$SE_hat ),
            df = c( SE$df ),
            p_value = two_sided_p(ATE_hat, SE_hat, df )
        )

    }
}









#' Estimate ATE via Middleton & Aronow unbiased approach
#'
#' This follows Middleton, J. A., & Aronow, P. M. (2015). Unbiased
#' Estimation of the Average Treatment Effect in Cluster-Randomized
#' Experiments. 6(1–2), 312. doi: 10.1515/spp-2013-0002
#'
#' They introduce a Raj Difference estimator that is an extension of a
#' Horvitz-Thompson estimator where, in effect, we divide the sum of
#' the outcomes by a fixed constant rather than the realized sample
#' size to avoid biasing our estimate.
#'
#' @inheritParams linear_model_estimators
#'
#' @return tibble of estimates using different varieties of the
#'   methods described in the paper.
#'
#' @export
middleton_aronow_estimator <- function( formula,
                                     data = NULL,
                                     control_formula = NULL,
                                     weight = c( "individual", "cluster" ),
                                     aggregated = FALSE ) {


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

    if ( !has_site ) {
        data_agg$siteID = "Uni"
    }
    n_site = length( unique( data_agg$siteID ) )

    M = nrow( data_agg )
    m_t = sum( data_agg$Z )

    tots <- data_agg %>% group_by( siteID ) %>%
        mutate( Ytot = Ybar * n ) %>%
        summarise( YT_1 = mean( Ytot[Z==1] ),
                   YT_0 = mean( Ytot[Z==0] ),
                   M = n(),
                   N = sum( n ),
                   ATE_hat = (M/N) * (YT_1 - YT_0) )
    stopifnot( nrow( tots ) == n_site )

    #ATE_hat <- (M / N) * ( tots$YT[[2]] - tots$YT[[1]] )

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












