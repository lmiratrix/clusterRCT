
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


    get_agg_ests <- function( M1, name = "MLM" ) {
        est1 <- M1$coefficients[["Z"]]
        se1 <- M1$std.error[["Z"]]
        pv1 <- M1$p.value[["Z"]]

        # Compile our results
        tibble(
            method = c( name ),
            ATE_hat = c( est1 ),
            SE_hat = c( se1 ),
            p_value = c( pv1 )
        )
    }


    if ( !is.null( formula ) && !aggregated ) {
        data = make_canonical_data( formula=formula, data=data,
                                    control_formula=control_formula )
    }

    datagg = data
    if ( !aggregated ) {
        datagg = aggregate_data(data, control_formula)
        control_formula = attr( datagg, "control_formula" )
    }

    needFE = "blockID" %in% names( data )

    form = make_regression_formula( Yobs = "Ybar",
                                    FE = needFE,
                                    control_formula = control_formula )

    M3 <- lm_robust( form, data=datagg, se_type = "HC2" )
    Agg_FE_cluster = get_agg_ests( M3, ifelse( needFE, "Agg_FE_Cluster", "Agg_noFE_Cluster" ) )


    M4 <- lm_robust( form, data=datagg, weights = n, se_type = "HC2" )
    Agg_FE_person = get_agg_ests( M4, ifelse( needFE, "Agg_FE_Person", "Agg_noFE_Person" ) )


    if ( !needFE ) {
        return( bind_rows( Agg_FE_cluster, Agg_FE_person ) )
    }

    formI = make_regression_formula( Yobs = "Ybar",
                                    FE = needFE, interacted = TRUE,
                                    control_formula = control_formula )
    M5 <- lm_robust( formI, data=datagg, se_type = "HC2" )
    Agg_FI = generate_all_interacted_estimates( M5, data,
                                                use_full_vcov=TRUE,
                                                method = "Agg_FI" )

    M6 <- lm_robust( formI, data=datagg, se_type = "HC2", weights = n )
    Agg_wFI = generate_all_interacted_estimates( M6, data,
                                                 use_full_vcov=TRUE,
                                                 method = "Agg_wFI" )

    # Compile our results
    res <- bind_rows( Agg_FE_cluster,
               Agg_FE_person,
               Agg_FI,
               Agg_wFI )

    # convert NaNs to NAs
    res$SE_hat[ is.nan(res$SE_hat) ] = NA
    res$p_value[ is.na(res$SE_hat) ] = NA

    res
}



if ( FALSE ) {

    data( fakeCRT )

    fakeCRT

    formula =  Yobs ~ T.x | S.id | D.id
    control_formula = ~ X.jk + C.ijk
    data = clusterRCT:::make_canonical_data(formula=formula, data=fakeCRT,
                                            control_formula = control_formula)

    aggregated = FALSE

    aggregation_estimators(formula, control_formula = control_formula,
                   data = fakeCRT )

    aggregation_estimators( Yobs ~ T.x | S.id | D.id,
                    data = fakeCRT )

    aggregation_estimators( Yobs ~ T.x | S.id,
                    data = fakeCRT )

}



