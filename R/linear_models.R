


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

    form = Y ~ 1 + Z
    if ( !is.null( data$siteID ) ) {
        form = Y ~ 0 + Z + siteID
    }

    # Linear regression
    M2 <- estimatr::lm_robust( form, data=data, clusters=clusterID )
    est2 <- M2$coefficients[["Z"]]
    se2  <- M2$std.error[["Z"]]
    pv2 <- M2$p.value[["Z"]]


    # If blocked, average two ways
    if ( !is.null( data$siteID ) ) {


    }


    # Compile our results
    tibble(
        method = c( "LR w/ CRSE" ),
        ATE_hat = c( est2 ),
        SE_hat = c( se2 ),
        p_value = c( pv2 )
    )

}











#' Interacted linear regression models
#'
#' These linear models have siteID by treatment interaction terms.  The final ATE
#' estimates are then weighted average of the block (site) specific ATE
#' estimates.
#'
#'#' If siteID passed, it will weight the RA blocks within site and then average
#' these site estimates.
#'
#' SEs come from the overall variance-covariance matrix.
#'
#' @inheritParams linear_model_estimators
#'
#' @return Dataframe of the different versions of this estimator (person and
#'   site weighted)
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

    # site weighting
    if (!is.null( siteID)) {
        # aggregate!
        wts <- data %>% group_by(B, siteID) %>%
            dplyr::summarise(n = n()) %>% group_by(siteID) %>% mutate(wts = n / sum(n))

        # some checks to make sure we are matching RA blocks and sites to the
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
    ATE_hat_site <- weighted.mean(ATE_hats, wts)
    # Calculate SE for ATE_hat_site
    SE_site <- sqrt(sum(wts ^ 2 * SE_hat))


    wts.indiv <- nj / n
    ATE_hat_indiv <- weighted.mean(ATE_hats, wts.indiv)
    SE_indiv <- sqrt(sum(wts.indiv ^ 2 * SE_hat))

    # This is the cautious way we don't need since we have 0s in the off diagonal
    # SE_site = sqrt( t(wts) %*% VC %*% wts )
    # faster way---this should work easily.
    # sqrt( t(wts.indiv) %*% VC %*% wts.indiv )

    interactModels <- data.frame(method = c("FE-Int-Sites", "FE-Int-Persons"),
                                 ATE_hat = c(ATE_hat_site, ATE_hat_indiv),
                                 SE = c(SE_site, SE_indiv), stringsAsFactors = FALSE)
    if (!is.null(control_formula)) {
        interactModels$method <- paste0(interactModels$method, "-adj")
    }
    interactModels
}







