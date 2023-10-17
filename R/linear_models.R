


#### Linear model estimators ####





#' Estimate ATEs for a cluster RCT using linear model
#'
#' @param formula Formula for outcome and treatment and nesting.  If
#'   NULL, data is assumed to be in canonical form (see vignette for
#'   further discussion).
#'
#' @importFrom estimatr lm_robust
#'
#' @export
linear_model_estimators <- function( formula,
                                     data = NULL,
                                     control_formula = NULL ) {

    require( estimatr )

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }

    # Linear regression (possibly with fixed effects if block exists)
    needFE = "blockID" %in% names(data)
    form = make_regression_formula( FE = needFE,
                                    control_formula = control_formula )

    M2 <- estimatr::lm_robust( form, data=data, clusters=clusterID )
    est2 <- M2$coefficients[["Z"]]
    se2  <- M2$std.error[["Z"]]
    pv2 <- M2$p.value[["Z"]]

    # Compile our results
    tibble(
        method = c( "LR w/ CRSE" ),
        ATE_hat = c( est2 ),
        SE_hat = c( se2 ),
        p_value = c( pv2 )
    )

}


calc_agg_estimate <- function( wts, ATE_hats, SE_hats ) {

    # normalize weights
    wts = wts / sum(wts)

    ATE_hat <- weighted.mean(ATE_hats, wts)
    SE_hat = NA
    if ( is.matrix(SE_hats) ) {
        wts = as.matrix(wts, ncol=1)
        SE_hat = sqrt( t(wts) %*% SE_hats %*% wts )
    } else {
        SE_hat = sqrt( sum(wts ^ 2 * SE_hats) )
    }
    tibble( ATE_hat = ATE_hat,
            SE_hat = SE_hat )
}


calculate_avg_impacts <- function( ATE_hats, SE_hats, clusterID, blockID ) {

    sizes <- count_block_sizes( clusterID, blockID )
    sizes
    stopifnot( all( sizes$bid == names(ATE_hats) ) )

    K = nrow( sizes )

    # equal weighting
    ATE_eq = calc_agg_estimate( rep(1,K), ATE_hats, SE_hats )

    # weight by number of clusters
    ATE_cluster <-calc_agg_estimate(sizes$J, ATE_hats, SE_hats)

    # weight by number of students
    ATE_indiv = calc_agg_estimate(sizes$n_stud, ATE_hats, SE_hats)

    bind_rows( equal = ATE_eq,
               cluster = ATE_cluster,
               indiv = ATE_indiv, .id="weight" )
}


#' Estimate ATEs for a cluster RCT using linear model
#'
#' Estimate district-average effects using an interacted FE model and
#' then aggregate.
#'
#' @inheritParams linear_model_estimators
#' @param use_full_vcov TRUE/FALSE. When calculating standard errors,
#'   should the full variance-covariance matrix of the block-level
#'   estimates be used, or just the diagonal of standard errors?
#' @importFrom estimatr lm_robust
#'
#' @export
interacted_linear_model_estimators <- function( formula,
                                     data = NULL,
                                     control_formula = NULL,
                                     use_full_vcov = FALSE ) {

    require( estimatr )

    if ( !is.null( formula ) ) {
        data = clusterRCT:::make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }

    # Linear regression (possibly with fixed effects if block exists)
    needFE = "blockID" %in% names(data)
    if ( !needFE ) {
        return(     tibble(
            method = c(),
            ATE_hat = c(  ),
            SE_hat = c(  ),
            p_value = c(  )
        ) )
    }
    form = make_regression_formula( FE = TRUE, interacted = TRUE,
                                    control_formula = control_formula )

    M0.int <- estimatr::lm_robust( form, data=data, clusters=clusterID )

    J = length( unique( data$blockID ) )
    ids <- grep( "Z:", names(coef(M0.int)))
    stopifnot(length(ids) == J)

    VC <- vcov(M0.int)

    ATE_hats <- coef(M0.int)[ids]
    names(ATE_hats) = str_remove(names(ATE_hats), "Z:blockID")

    if ( use_full_vcov ) {
        # This is the cautious way we don't need since we have 0s in the off diagonal
        ests <- calculate_avg_impacts( ATE_hats, VC[ids,ids],
                           data$clusterID, data$blockID )
    } else {
        SE_hats <- diag(VC)[ids]
        ests <- calculate_avg_impacts( ATE_hats, SE_hats,
                                       data$clusterID, data$blockID )

    }

    ests <- ests %>%
        mutate( method = c( "LR-FI-CRVE-Block", "LR-FI-CRVE-Cluster", "LR-FI-CRVE-Person" ) ) %>%
        relocate( method )

    if (!is.null(control_formula)) {
        ests$method <- paste0(ests$method, "-adj")
    }

    ests
}





if ( FALSE ) {

    data( fakeCRT )
    fakeCRT

    interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id, control_formula = ~ X.jk + C.ijk,
                             data = fakeCRT )

    interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                                        data = fakeCRT )

    interacted_linear_model_estimators( Yobs ~ T.x | S.id,
                                        data = fakeCRT )

}






#' Interacted linear regression models
#'
#' These linear models have blockID by treatment interaction terms.
#' The final ATE estimates are then weighted average of the block
#' (block) specific ATE estimates.
#'
#' SEs come from the overall variance-covariance matrix.
#'
#' @inheritParams linear_model_estimators
#'
#' @return Dataframe of the different weighting versions of this
#'   estimator (by number of schools (clusters) or simple average)
#' @export
#'
interacted_aggregated_estimators <- function( formula, data = NULL,
                                              control_formula = NULL) {

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula, data, control_formula )
    }

    data$B <- droplevels(as.factor(data$B))
    J <- length(unique(data$B))
    nj <- table(data$B)
    n <- nrow(data)

    formula <- make_FE_int_formula("Yobs", "Z", "B", control_formula, data)
    M0.int <- lm(formula, data = data)
    ids <- grep( "Z:", names(coef(M0.int)))
    stopifnot(length(ids) == J)
    VC <- vcov(M0.int)
    ATE_hats <- coef(M0.int)[ids]

    # block weighting
    if (!is.null( blockID)) {
        # aggregate!
        wts <- data %>% group_by(B, blockID) %>%
            dplyr::summarise(n = n()) %>% group_by(blockID) %>% mutate(wts = n / sum(n))

        # some checks to make sure we are matching RA blocks and blocks to the
        # right things
        stopifnot( nrow( wts ) == J )
        nms <- gsub( "Z:B", "", names(ATE_hats))
        stopifnot(all(nms == wts$B))
        wts <- wts$wts / sum(wts$wts)
    } else {
        wts <- rep(1 / J, J)
    }

    # the block SEs from our linear model
    SE_hat <- diag(VC)[ids]
    ATE_hat_block <- weighted.mean(ATE_hats, wts)
    # Calculate SE for ATE_hat_block
    SE_block <- sqrt(sum(wts ^ 2 * SE_hat))


    wts.indiv <- nj / n
    ATE_hat_indiv <- weighted.mean(ATE_hats, wts.indiv)
    SE_indiv <- sqrt(sum(wts.indiv ^ 2 * SE_hat))

    # This is the cautious way we don't need since we have 0s in the off diagonal
    # SE_block = sqrt( t(wts) %*% VC %*% wts )
    # faster way---this should work easily.
    # sqrt( t(wts.indiv) %*% VC %*% wts.indiv )

    interactModels <- data.frame(method = c("FE-Int-blocks", "FE-Int-Persons"),
                                 ATE_hat = c(ATE_hat_block, ATE_hat_indiv),
                                 SE = c(SE_block, SE_indiv), stringsAsFactors = FALSE)
    if (!is.null(control_formula)) {
        interactModels$method <- paste0(interactModels$method, "-adj")
    }
    interactModels
}







