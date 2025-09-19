
# MRStdCRT estimator wrapper




#' Estimate ATEs for a cluster RCT using the MRStdCRT pacakge
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
#' @export
MRStdCRT_estimator <- function( formula,
                            data = NULL,
                            control_formula = NULL,
                            weight = c( "Person", "Cluster" ) ) {

    warning( "This method does not yet work due to difficulties mapping to the call" )
    require( MRStdCRT )

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

    res
    names( example )

    example$estimate

    example$m


    library( clusterRCT )
    describe_clusterRCT( PEGS ~ INTERVENTION | CLUST, data = ppact )
    compare_methods( PEGS ~ INTERVENTION | CLUST, data = ppact )

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




