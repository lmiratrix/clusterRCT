

# GEE estimators





#' Estimate ATEs for a cluster RCT using a GEE approach
#'
#' @param formula Formula for outcome and treatment and nesting.  If
#'   NULL, data is assumed to be in canonical form (see vignette for
#'   further discussion).
#'
#' @importFrom geepack geeglm
#'
#' @export
gee_estimators <- function( formula,
                                     data = NULL,
                                     control_formula = NULL,
                                     weight = c( "Person", "Cluster" ) ) {

    weight = match.arg(weight)
    est_method = ifelse( weight == "Person", "GEE", "GEEcw" )

    require( estimatr )

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }


    # Make our weights variable
    if (weight == "Person") {
        data$.weight <- 1
    } else {
        stop( "weighting not yet implemented for GEE approach" )
        data <- data %>%
            group_by( clusterID ) %>%
            mutate( .weight = 1 / n() ) %>%
            ungroup()
    }


    # Linear regression (possibly with fixed effects if block exists)
    needFE = "blockID" %in% names(data)
    form = make_regression_formula( FE = needFE,
                                    control_formula = control_formula )


    # Fix the factors
    data$clusterID = as.numeric( as.factor( data$clusterID ) )
    data <- arrange( data, clusterID )

    gee_model <- geepack::geeglm(
        formula = form,
        id = clusterID,
        data = data,
        family = gaussian,          # or binomial, etc.
        corstr = "exchangeable"     # assumes common intra-cluster correlation
    )
    gee_model

    npar = length( gee_model$coefficients )
    tt = broom::tidy( gee_model ) %>%
        filter( term == "Z")

    # Compile our results
    nm = paste0( est_method,
                 ifelse( needFE, "-FE", "" ) )
    df2 = n_distinct( data$clusterID ) - npar

    tibble(
        method = c( nm ),
        ATE_hat = c( tt$estimate ),
        weight = weight,
        SE_hat = c( tt$std.error ),
        p_value = c( tt$p.value ),
        df = c( df2 )
    )

}




