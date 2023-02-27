



#' Cluster RCT with multi-level modeling
#'
#' This uses the lmerTest package to get the p-values, and fits a
#' model with assumed homoskedasticity, etc.  I.e., this is the
#' vanilla MLM that one would typically fit.
#'
#' @inheritParams linear_model_estimators
#'
#' @export
MLM_estimators <- function( formula,
                                     data = NULL,
                                     control_formula = NULL ) {

    require( lme4 )
    require( lmerTest)

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }


    form = Y ~ 1 + Z + (1|clusterID)
    if ( !is.null( data$siteID ) ) {
        form = Y ~ 0 + Z + siteID + (1|clusterID)
    }

    M1 <- lmerTest::lmer( form, data=data )

    est1 <- fixef( M1 )[["Z"]]
    se1 <- arm::se.fixef( M1 )[["Z"]]
    pv1 <- summary(M1)$coefficients["Z",5]



    # Compile our results
    tibble(
        method = c( "MLM" ),
        ATE_hat = c( est1 ),
        SE_hat = c( se1 ),
        p_value = c( pv1 )
    )

}


