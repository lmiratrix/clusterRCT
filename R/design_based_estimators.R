#
# Design-based methods follow those found in Schochet JASA paper
# "Design-Based Ratio Estimators and Central Limit Theorems for
# Clustered, Blocked RCTs" Peter Z. Schochet, Nicole E. Pashley, Luke
# W. Miratrix, and Tim Kautz
#
#
# Also implements the Middleton and Aronow design based Raj estimator.



get_overall_ATE <- function( mod, weights ) {
    cc = coef(mod)

    J = length( weights )

    ids <- grep( "Z:", names(cc) )
    stopifnot(length(ids) == J)
    ATE_hats <- cc[ids]

    weighted.mean( ATE_hats, w = weights )
}



# Notes on Schochet paper notation
#
# m_b^1 = Total number of clusters assigned to tx group in block b




# Implement formula (9) and following equations from Schochet et
# al paper.
#
# See page 2139-40 (section 4.2)
#
# @param adt Aggregated (cluster-level) data.
# @param v Degrees of freedom.  Use v=# covariates.
schochet_variance_formula_block <- function( adt, v ) {

    require( tidyr )
    require( dplyr )

    stopifnot( !is.null( adt$resid ) )
    stopifnot( !is.null( adt$Z ) )

    calc_s2_partial <- function( .weight, resid ) {

        w_bar = mean( .weight )
        tot = sum( ( (.weight^2) / (w_bar^2) ) * resid^2 )
        #mb = nrow( adt )

        #tot / (mb - v * p * q - 1)
        tot
    }

    adtw <- adt %>%
        group_by( blockID, Z ) %>%
        summarize( m = n(),
                   n = sum( n ),
                   totwt = sum( .weight ),
                   s2 = calc_s2_partial( .weight, resid ),
                   .groups = "drop" ) %>%
        dplyr::group_by( blockID ) %>%
        dplyr::mutate( p = m / sum(m) ) %>%
        ungroup() %>%
        dplyr::mutate( q = totwt / sum(totwt) ) %>%
        dplyr::group_by( Z ) %>%
        dplyr::mutate( wt = totwt / sum(totwt),
                       s2 = s2 / (m - v * p * q - 1) ) %>%
        ungroup()

    pt2 <- adtw %>%
        dplyr::select( -p, -q ) %>%
        tidyr::pivot_wider( names_from = "Z",
                            values_from = c( totwt, s2, m, n, wt) ) %>%
        dplyr::mutate( varD = s2_1 / m_1 + s2_0 / m_0,
                       SE_hat = sqrt( varD ),
                       n = n_0 + n_1,
                       m = m_0 + m_1,
                       weight = totwt_0 + totwt_1 )

    pt2$weight = pt2$weight / sum(pt2$weight)   # normalize weights to avoid overflow

    pt2
}



#' Given block specific estimates, calculate the ATE, associated
#' standard error, and design-based degrees of freedom.
#'
#' @param block_estimates Result of the
#'   schochet_variance_formula_block() function call.
#' @return List of ATE, SE, and degrees of freedom.
schochet_overall_estimator <- function( block_estimates ) {

    h = length( unique( adt$blockID ) )

    SE_ATE <- with(block_estimates, sum(varD*totwt^2) / (h*mean(totwt))^2)

    list( SE_hat = sqrt( SE_ATE ),
          df = nrow(adt) - 2*h - 1 )

}



# From Schochet et al equation (15) for the restricted model with no
# interaction term for tx and block.
schochet_FE_variance_formula <- function( adt, v ) {

    stopifnot( !is.null( adt$resid ) )
    stopifnot( !is.null( adt$Z ) )

    m = nrow(adt)
    h = length( unique(adt$blockID) )

    SE_ATE <- adt %>%
        group_by( blockID ) %>%
        summarize( p_b = mean(Z),
                   m_b = n(),
                   wbar_b = mean(.weight),
                   num = sum(.weight^2 * (Z-p_b)^2 * resid^2)) %>%
        summarize(SE_ATE = sqrt( m / (m-h-v-1) * sum(num) /
                                     sum(m_b*p_b*(1-p_b)*wbar_b)^2 )) %>%
        pull(SE_ATE)

    df = m - v - h - 1

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
                                     weight = c( "Person", "Cluster" ),
                                     aggregated = FALSE,
                                     include_block_estimates = FALSE ) {

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
    has_block = "blockID" %in% names( data )


    # Make our weights variable
    suff = ""
    if (weight == "Person") {
        data_agg$.weight <- data_agg$n
        suff = "Person"
    } else {
        # data_agg$.weight <- 1 / data_agg$n
        data_agg$.weight <- 1
        suff = "Cluster"
    }

    v = number_controls(control_formula)
    J = length( unique( data$clusterID ) )

    # block level stuff: only if we have a blockID do we fit an
    # interacted model or allow for fixed effects for each block.
    ATE_int = NA
    if ( has_block ) {
        K = length( unique( data$blockID ) )

        # Fully interacted model
        form = make_regression_formula( Yobs = "Ybar", interacted = TRUE,
                                        control_formula = control_formula )
        mod = lm( form, data=data_agg, weights = .weight )

        data_agg$resid = residuals(mod)
        block_tab = schochet_variance_formula_block( data_agg, v = v )

        # Translate Schochet notation to our notation.
        block_tab <- rename( block_tab,
                             J = m )

        ests = generate_all_interacted_estimates( mod, data = data_agg,
                                                  use_full_vcov = FALSE,
                                                  SE_table = block_tab,
                                                  method = glue::glue("DB_FI_{suff}"),
                                                  aggregated = TRUE,
                                                  include_block_estimates = include_block_estimates )

        # df = m - 2h - v* = #clusters - 2 #blocks - #covariates
        ests$df = J - 2*K - v

        # fixed effects model
        form = make_regression_formula( Yobs = "Ybar", FE = TRUE,
                                        control_formula = control_formula )
        mod_FE = lm( form, data=data_agg, weights = .weight )
        ATE_FE = coef(mod_FE)[ "Z" ]
        data_agg$resid = residuals(mod_FE)
        SE_FE = schochet_FE_variance_formula( data_agg, v = v )

        ests$weight = paste( suff, ests$weight, sep="_" )

        FE <- tibble(
            method = c( glue::glue( "DB_FE_{suff}" ) ),
            weight = suff,
            ATE_hat = c( ATE_FE ),
            SE_hat = c( SE_FE$SE_hat ),
            df = c( SE_FE$df ),
            p_value = two_sided_p( ATE_hat, SE_hat, df )
        )

        # Pack up results
        ests <- bind_rows( ests, FE )

        if ( include_block_estimates ) {
            blocks <- attr( ests, "blocks" )
            blocks = left_join( blocks, block_tab, by="blockID" )
            attr( ests, "blocks" ) <- blocks
        }

        ests

    } else {
        # No Blocks, so we use simple Cluster randomized estimators.
        form = make_regression_formula( Yobs = "Ybar", FE = FALSE,
                                        control_formula = control_formula )
        mod = lm( form, data=data_agg, weights = .weight )
        ATE_FE = coef(mod)[ "Z" ]

        data_agg$blockID = "Ind"
        data_agg$resid = residuals( mod )
        SE = schochet_variance_formula_block( data_agg, v = v )

        df = nrow(data_agg) - 2 - v

        tibble(
            method = glue::glue( "DB_{suff}" ),
            ATE_hat = c( ATE_FE ),
            SE_hat = c( SE$SE_hat ),
            df = df,
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
    has_block = "blockID" %in% names( data )

    if ( !has_block ) {
        data_agg$blockID = "Uni"
    }
    n_block = length( unique( data_agg$blockID ) )

    M = nrow( data_agg )
    m_t = sum( data_agg$Z )

    tots <- data_agg %>% group_by( blockID ) %>%
        mutate( Ytot = Ybar * n ) %>%
        summarise( YT_1 = mean( Ytot[Z==1] ),
                   YT_0 = mean( Ytot[Z==0] ),
                   M = n(),
                   N = sum( n ),
                   ATE_hat = (M/N) * (YT_1 - YT_0) )
    stopifnot( nrow( tots ) == n_block )

    #ATE_hat <- (M / N) * ( tots$YT[[2]] - tots$YT[[1]] )

    # Make our weights variable
    suff = ""
    if (weight == "individual") {
        data_agg$.weight <- data_agg$n
        suff = "Person"
    } else {
        data_agg$.weight <- 1 / data_agg$n
        suff = "Cluster"
    }

    # block level stuff: only if we have a blockID do we fit an
    # interacted model or allow for fixed effects for each block.
    ATE_int = NA
    if ( has_block ) {
        # Fully interacted model
        form = make_regression_formula( Yobs = "Ybar", interacted = TRUE,
                                        control_formula = control_formula )
        mod = lm_robust( form, data=data_agg, weights = .weight,
                         ci = FALSE, se_type="none" )
        stwt <- data_agg %>% group_by( blockID ) %>%
            summarise( wt = sum( .weight ) )
        ATE_int = get_overall_ATE(mod, stwt$wt )

        # fixed effects model
        form = make_regression_formula( Yobs = "Ybar", FE = TRUE,
                                        control_formula = control_formula )
        mod_FE = lm_robust( form, data=data_agg, weights = .weight,
                            ci = FALSE, se_type="none" )
        ATE_FE = coef(mod_FE)[ "Z" ]
        tibble(
            method = c( glue::glue( "DB_{suff} (int)" ),
                        glue::glue( "DB_{suff} (FE)" ) ),
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
            method = glue::glue( "DB_{suff}" ),
            ATE_hat = c( ATE_FE ),
            SE_hat = c( NA ),
            p_value = c( NA )
        )

    }
}



#### Testing #####

if ( FALSE ) {


    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 30                          # number of schools
        , K = 10                          # number of districts
        , nbar = 10                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of district assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.1              # percent of district variation
        , ICC.3 = 0.2 # district intraclass correlation
        , omega.3 = 0.1           # ratio of district effect size variability
        , R2.2 = 0.1              # percent of school variation
        , ICC.2 = 0.2             # school intraclass correlation
        , omega.2 = 0.1           # ratio of school effect size variability
        , R2.1 = 0.1    # percent of indiv variation explained
    )


    sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
    rm( model.params.list )

    head( sim.data )
    table( sim.data$D.id)
    aa = design_based_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )
    aa

    # No block-level.
    sim.data$D.id = NULL
    bb = design_based_estimators( Yobs ~ T.x | S.id, data=sim.data )
    bb


}











