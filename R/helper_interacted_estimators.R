

calc_agg_estimate <- function( wts, ATE_hats, SE_hats ) {

    # normalize weights
    wts = wts / sum(wts)

    ATE_hat <- weighted.mean(ATE_hats, wts)
    SE_hat = NA
    if ( is.matrix(SE_hats) ) {
        wts = as.matrix(wts, ncol=1)
        SE_hat = as.numeric( sqrt( t(wts) %*% SE_hats %*% wts ) )
    } else {
        SE_hat = sqrt( sum(wts ^ 2 * SE_hats) )
    }
    tibble( ATE_hat = ATE_hat,
            SE_hat = SE_hat,
            p_value = 2*pnorm( abs(ATE_hat) / SE_hat, lower.tail=FALSE ) )
}



calculate_avg_impacts <- function( ATE_hats, SE_hats,
                                   sizes = NULL,
                                   clusterID = NULL,
                                   blockID = NULL ) {
    if ( is.null( sizes ) ) {
        stopifnot( !is.null( clusterID ) && !is.null( blockID ) )
    } else {
        stopifnot( is.null( clusterID ) && is.null( blockID ) )
    }

    if ( is.null( sizes ) ) {
        sizes <- count_block_sizes( clusterID, blockID )
    }

    stopifnot( all( sizes$blockID == names(ATE_hats) ) )

    K = nrow( sizes )

    # equal weighting
    ATE_eq <- calc_agg_estimate( rep(1,K), ATE_hats, SE_hats )

    # weight by number of clusters
    ATE_cluster <- calc_agg_estimate(sizes$J, ATE_hats, SE_hats)

    # weight by number of students
    ATE_indiv <- calc_agg_estimate(sizes$n, ATE_hats, SE_hats)

    rsp <- bind_rows( Block = ATE_eq,
                      Cluster = ATE_cluster,
                      Person = ATE_indiv, .id="weight" )
    if ( "weight" %in% names(sizes) ) {
        # weight by arbitrary weight, if given
        ATE_weight = calc_agg_estimate(sizes$weight, ATE_hats, SE_hats)
        ATE_weight$weight = "Weight"
        rsp = bind_rows( rsp, ATE_weight )
    }

    return( rsp )
}


#' Average block-level estimates to get overall ATE
#'
#' Also calculate associated standard error.
#'
#' @inheritParams interacted_linear_model_estimators
#' @param fitModel The model that has the interacted estimates in it.
#' @param method Prefix of the method.  Will add the weighting to
#'   stem.
#' @param aggregated TRUE means data is summarized at the Block level.
#'   False means it is not, and needs to be aggregated.
#' @export
generate_all_interacted_estimates <- function( fitModel, data,
                                               use_full_vcov = FALSE,
                                               SE_table = NULL,
                                               method = "LR_FI_CRVE",
                                               aggregated = FALSE,
                                               include_block_estimates = FALSE ) {
    J = length( unique( data$blockID ) )

    cc = NA
    if ( !( class(fitModel) %in% c( "lm", "lm_robust" ) ) ) {
        # Multilevel model!
        cc = fixef( fitModel )
    } else {
        cc = coef( fitModel )
    }
    ids <- grep( "Z:", names(cc) )
    stopifnot(length(ids) == J)

    VC <- vcov(fitModel)

    ATE_hats <- cc[ids]
    names(ATE_hats) = stringr::str_remove(names(ATE_hats), "Z:blockID")

    sizes = data
    if ( !aggregated ) {
        sizes <- count_block_sizes( data$clusterID, data$blockID )
    }

    ests = NA
    SE_hats = NA
    if ( !is.null( SE_table ) ) {
        stopifnot( all( names(ATE_hats) == SE_table$blockID ) )
        SE_hats = SE_table$SE_hat
        ests <- calculate_avg_impacts( ATE_hats,
                                       SE_hats,
                                       sizes = SE_table )
    } else if ( use_full_vcov ) {
        # This is the cautious way we don't need since we have 0s in
        # the off diagonal
        SE_hats = diag( as.matrix( VC ) )[ids] # for table below
        ests <- calculate_avg_impacts( ATE_hats, VC[ids,ids],
                                       sizes = sizes )
    } else {
        SE_hats <- diag( as.matrix( VC ) )[ids]
        ests <- calculate_avg_impacts( ATE_hats, SE_hats,
                                       sizes = sizes )
    }


    ests <- ests %>%
        mutate( method = paste0( method, "_", weight ) ) %>%
        relocate( method )

    if ( include_block_estimates ) {
        tb = tibble( blockID = names(ATE_hats),
                     ATE_hat = ATE_hats,
                     SE_hat = SE_hats )
        attr( ests, "blocks" ) <- tb
    }
    ests
}
