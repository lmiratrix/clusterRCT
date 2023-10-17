
#' Estimate ATEs for cluster RCT using data aggregation.
#'
#' This method aggregates data at the cluster level, and analyzes
#' using linear regression. Can optionally pass pre-aggregated data if
#' desired.
#'
#' @inheritParams linear_model_estimators
#'
#' @export
aggregation_estimators <- function( formula,
                                    data = NULL,
                                    control_formula = NULL,
                                    aggregated = FALSE ) {

    if ( !is.null( formula ) && !aggregated ) {
        data = make_canonical_data( formula=formula, data=data,
                                    control_formula=control_formula )
    }

    datagg = data
    if ( !aggregated ) {
        datagg = aggregate_data(data, control_formula)
        control_formula = attr( datagg, "control_formula" )
    }

    form = make_regression_formula( Yobs = "Ybar",
                                    FE = "blockID" %in% names( data ),
                                    control_formula = control_formula )

    M3 <- lm_robust( form, data=datagg, se_type = "HC2" )
    est3 <- M3$coefficients[["Z"]]
    se3 <- M3$std.error[["Z"]]
    pv3 <- M3$p.value[["Z"]]


    M4 <- lm_robust( form, data=datagg, weights = n, se_type = "HC2" )
    est4 <- M4$coefficients[["Z"]]
    se4 <- M4$std.error[["Z"]]
    pv4 <- M4$p.value[["Z"]]

    # Compile our results
    tibble(
        method = c( "LR (agg)", "LR (agg, wt)" ),
        ATE_hat = c( est3, est4 ),
        SE_hat = c( se3, se4 ),
        p_value = c(pv3, pv4 )
    )


}


