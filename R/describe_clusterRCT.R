


#' Describe a cluster RCT dataset
#'
#' Calculate various statistics describing total clusters, etc.
#'
#' @param formula
#'
#' @export
describe_clusterRCT <- function( formula = NULL,
                                 data = NULL,
                                 control_formula = NULL ) {

    data = make_canonical_data( formula=formula, data=data,
                                control_formula = control_formula )

    cnames = NULL
    if ( !is.null( control_formula ) ) {
        cnames = setdiff( colnames(data), c("Yobs", "Z", "clusterID", "siteID" ) )
    }

    K = length( unique( data$siteID ) )

    sizes = data %>%
        group_by( siteID, clusterID, Z ) %>%
        summarise( n = n(), .groups = "drop" )

    cstat = sizes %>%
        summarise( J = n(),
                   nbar = mean(n),
                   ncv = sd(n) / nbar,
                   n.25 = quantile( n, 0.25 ),
                   n.75 = quantile( n, 0.75 ),
                   n.IQR = n.75 - n.25,
                   p.tx = mean(Z) )

    sstat = sizes %>% group_by( siteID ) %>%
        summarise( J = n(),
                   n = sum(n),
                   p.tx = mean(Z) ) %>%
        summarise( K = n(),
                   Jbar = mean(J),
                   Jcv = sd(J) / Jbar,
                   J.25 = quantile( J, 0.25 ),
                   J.75 = quantile( J, 0.75 ),
                   J.IQR = J.75 - J.25,
                   n_site = mean(n),
                   n_site_cv = sd( n ) / n_site,
                   tx.avg = mean( p.tx ),
                   tx.cv = sd( p.tx ) / tx.avg,
                   tx.25 = quantile( p.tx, 0.25 ),
                   tx.75 = quantile( p.tx, 0.75 ),
                   tx.IQR = tx.75 - tx.25 )


    library( lme4 )
    form = Yobs ~ 1 + Z + (1 | siteID ) + (1 | clusterID )
    #if ( !is.null( control_formula ) ) {
    #    form = update( control_formula, Yobs ~ . + 1 + Z + (1 | siteID ) + (1 | clusterID ) )
    #}
    M = lmer( form, data=data )
    #arm::display(M)

    a = VarCorr( M )

    C.ICC = as.numeric( a$clusterID )
    S.ICC = as.numeric( a$siteID )
    sigma = sigma( M )
    tvar = (C.ICC + S.ICC + sigma)
    C.ICC = C.ICC / tvar
    S.ICC = S.ICC / tvar

    stats <- tibble( n = nrow(data) )
    stats = cbind( stats, cstat, sstat )
    stats$cluster_ICC = C.ICC
    stats$site_ICC = S.ICC
    stats$sdY0 = sd( data$Yobs[ data$Z == 0] )

    class( stats ) <- c( "clusterRCTstats", class( stats ) )


    # R2 values
    if ( !is.null( cnames ) ) {
        stats <- bind_cols( stats, calc_covariate_R2s(data) )
    }

    stats
}


#' Calculate R2 of level 1, level 2, and district.
#'
#' District is always the fixed effects for district.  No covariates
#' allowed at that level.  In all models district fixed effects are
#' included, so level 2 covarites R2 are additional explanatory power
#' beyond their representing district differences (e.g., they are
#' considered district centered).
#'
#' @importFrom purrr map_dbl
#' @param data in canonical form.
#' @param pooled Calculate pooled Tx and Co groups, attempting to
#'   remove systematic shift of average treatment effect.  If FALSE,
#'   subset to control observations only.
#'
calc_covariate_R2s <- function( data, pooled = FALSE ) {

    # Covert covariates to numerical, if they are not already.
    dm = model.matrix( Yobs ~ . - clusterID - siteID, data=data ) %>%
        as.data.frame() %>%
        dplyr::select( -`(Intercept)` )

    cnames = colnames(dm)

    dm$clusterID = data$clusterID
    dm$siteID = data$siteID

    # Add group means and group mean-centered versions of all
    # covariates
    result <- dm %>%
        group_by(siteID, clusterID) %>%
        mutate(
            across( everything(),
                    list(mn = ~ mean(.),
                         cent = ~. - mean(.)),
                    .names = "{.col}_{.fn}")
        ) %>%
        ungroup() %>%
        dplyr::select( -all_of( cnames ) )

    # Drop centered covariates with no variation (i.e., level 2
    # covariates get dropped from level 1 consideration)
    cents = which( endsWith( colnames(result), "_cent" ) )
    sds = purrr::map_dbl( cents, ~ sd( result[[ . ]] ) )
    ncov.1 = sum( sds != 0 )
    result = result[ -c( cents[sds == 0] ) ]
    ncov.2 = sum( endsWith( colnames(result), "_mn" ) ) - 1

    result$Yobs = data$Yobs

    # Two ways of calculating R2 values.
    if ( pooled ) {
        Msite = lm( Yobs ~ Z_mn*siteID, data=result )
        R2.site = summary(Msite)$adj.r.squared

        cents = result %>% dplyr::select( !ends_with( "_mn" ) )
        cents$Z_mn = result$Z_mn
        M = lm( Yobs ~ . + Z_mn*siteID - clusterID, data=cents )
        R2.cent = summary(M)$adj.r.squared

        mns = result %>% dplyr::select( !ends_with("_cent" ) )
        Mmn = lm( Yobs ~ . + Z_mn*siteID - clusterID, data=mns )
        R2.mn = summary( Mmn )$adj.r.squared

        Mtx = lm( Yobs ~ Z_mn, data=result )
        sTx <- summary(Mtx)
        tx.R2 = sTx$r.squared
        rat = 1 / (1 - sTx$r.squared)

        tibble( R2.1 = (R2.cent - R2.site)*rat,
                ncov.1 = ncov.1,
                R2.2 = (R2.mn - R2.site)*rat,
                ncov.2 = ncov.2,
                R2.3 = (R2.site - tx.R2)*rat,
                R2.adj.rat = rat )
    } else {
        resCo = dplyr::filter( result, Z_mn == 0 )
        resCo$Z_mn = NULL

        Msite = lm( Yobs ~ siteID, data=resCo )
        R2.site = summary(Msite)$adj.r.squared

        cents = resCo %>% dplyr::select( !ends_with( "_mn" ) )
        M = lm( Yobs ~ . - clusterID, data=cents )
        summary( M )
        R2.cent = summary(M)$adj.r.squared

        mns = resCo %>% dplyr::select( !ends_with("_cent" ) )
        Mmn = lm( Yobs ~ . - clusterID, data=mns )
        summary( Mmn )
        R2.mn = summary( Mmn )$adj.r.squared

        tibble( R2.1 = R2.cent - R2.site,
                ncov.1 = ncov.1,
                R2.2 = R2.mn - R2.site,
                ncov.2 = ncov.2,
                R2.3 = R2.site )
    }
}



#' Make table of statistics by district
#'
#' Given individual/school/district data of a blocked, cluster RCT,
#' calculate statistics for each block (district).
#'
#' @param check_data_integrity TRUE means runs some checks and give
#'   errors if data fails them (e.g., incorrectly processed treatment
#'   vector.). FALSE means calculate statistics without these checks.
#' @export
make_site_table <- function(  formula = NULL,
                              data = NULL,
                              control_formula = NULL,
                              check_data_integrity = FALSE ) {

    data = make_canonical_data( formula=formula, data=data,
                                control_formula = control_formula )

    if ( check_data_integrity ) {
        check_data_integrity( data )
    }

    n <- nrow( data )

    K = length( unique( data$siteID ) )

    sizes = data %>%
        group_by( siteID, clusterID, Z ) %>%
        summarise( n = n(), .groups = "drop" )

    sstat = sizes %>% group_by( siteID ) %>%
        summarise( J = n(),
                   nbar = mean(n),
                   ncv = sd(n) / nbar,
                   n.25 = quantile( n, 0.25 ),
                   n.75 = quantile( n, 0.75 ),
                   n.IQR = n.75 - n.25,
                   p.tx = mean(Z) )
    sstat
}





#' Pretty print result from describe cluster RCT call
#'
#' @export
#' @param x A clusterRCTstats object.
#' @param ... No extra options passed.
#' @family clusterRCTstats
#'
print.clusterRCTstats <- function( x, ... ) {
    stopifnot( is.clusterRCTstats(x) )

    if (  is.null( x$K ) || x$K == 1 ) {
        scat( "Cluster RCT: %d units in %d clusters\n", x$n, x$J )
    } else {
        scat( "Cluster RCT: %d units in %d clusters across %d sites\n",
              x$n, x$J, x$K )
    }


    if ( !is.null( x$K ) && x$K != 1 ) {
        scat( "Site Statistics:\n\tAvg units/site: %.2f (coef var %.2f)\n\tAvg clusters/site: %.2f (coef var %.2f)\n\t25-75 Quantiles: %.1f -- %.1f w/ IQR = %.1f\n\tICC: %.2f\n\tAvg tx: %.2f (coef var %.2f)\n\t25-75 Quantile tx: %.2f -- %.2f w/ IQR %.2f\n",
              x$n_site, x$n_site_cv,
              x$Jbar, x$Jcv,
              x$J.25, x$J.75, x$J.IQR,
              x$site_ICC,
              x$tx.avg, x$tx.cv,
              x$tx.25, x$tx.75, x$tx.IQR)
        if ( !is.null( x$R2.3 ) ) {
            scat( "\tR2.3: %.2f (%d fixed effects)\n", x$R2.3, x$K)
        }
    }

    scat( "Cluster Statistics:\n\tAvg units: %.2f (coef var %.2f)\n\t25-75 Quantiles: %.1f -- %.1f w/ IQR = %.1f\n\tICC: %.2f\n",
          x$nbar, x$ncv,
          x$n.25, x$n.75, x$n.IQR,
          x$cluster_ICC )
    if ( !is.null( x$R2.2 ) ) {
        scat( "\tR2.2: %.2f (%d covariates)\n", x$R2.2, x$ncov.2)
    }


    scat( "Unit Statistics:\n\tprop clusters tx: %.2f\n\tstddev( Y0 ): %.2f\n",
          x$p.tx, x$sdY0 )
    if ( !is.null( x$R2.1 ) ) {
        scat( "\tR2.1: %.2f (%d covariates)\n", x$R2.1, x$ncov.1)
    }


}



#' Is object a clusterRCTstats object?
#'
#' @export
#' @aliases clusterRCTstats
#' @param x the object to check.
#' @family clusterRCTstats
is.clusterRCTstats = function( x ) {
    inherits(x, "clusterRCTstats")
}




#' Cast cluster RCT info result to data.frame
#'
#' @export
#' @aliases clusterRCTstats
#' @param x the clusterRCTstats object to covert
#' @family clusterRCTstats
as.data.frame.clusterRCTstats = function( x ) {
    class(x) = "list"
    as.data.frame( x )
}




#### Testing code ####

if ( FALSE ) {



    library( tidyverse )
    library( clusterRCT )

    set.seed( 1039 )


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
        , ICC.3 = 0.2             # district intraclass correlation
        , omega.3 = 0.2           # ratio of district effect size variability
        , R2.2 = 0.1              # percent of school variation
        , ICC.2 = 0.2             # school intraclass correlation
        , omega.2 = 0.0          # ratio of school effect size variability
        , R2.1 = 0.1    # percent of indiv variation explained
    )

    data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )

    data = slice_sample( data, n = nrow(data) / 2 )
    dd = sample( unique(data$S.id), 10 )
    data = filter( data, !( S.id %in% dd ) )
    data <- clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=data )
    head( data )
    describe_clusterRCT( data )

    data( "fakeCRT" )
    head( fakeCRT )
    dsc = describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT )
    class( dsc )
    print( dsc )

    # Check R2 calculations
    head( fakeCRT )
    dsc = describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                               control_formula = ~ V.k + X.jk + C.ijk )
    as.data.frame(dsc)
    print( dsc )



    dd = data.frame( cid = rep( 1:10, 1:10 ) )
    dd = mutate( dd,
                 Ybar = cid,
                 Z = as.numeric( cid %% 3 == 1 ) )
    dd2 = bind_rows( dd, dd )
    dd$cid = dd$cid + 100
    dd3 = bind_rows( dd, dd2 )
    dd$cid = dd$cid - 100
    dd2$Ybar = dd2$Ybar + 1
    dd3$Ybar = dd3$Ybar + 2
    dd = bind_rows( A = dd, B = dd2, C = dd3, .id = "sid" )
    dd$Y = dd$Ybar + rnorm( nrow(dd) )
    dd <- mutate( dd,
                  Z = ifelse( cid <= as.numeric(as.factor(sid))*2, 1, 0 ) )
    head(dd)
    skimr::skim(dd)

    make_site_table(  Y ~ Z | cid | sid, data=dd )
    describe_clusterRCT( Y ~ Z | cid | sid, data=dd)


    dd <- dd %>%
        group_by( sid ) %>%
        mutate(   X1 = Y + rnorm(n(), sd=3),
                  X2 = cid,
                  X3 = mean(Y),
                  X4 = sample( c( "A", "B", "C"), n(), replace=TRUE ),
                  Y = Y + 1 * Z ) %>%
        ungroup()

    data = clusterRCT:::make_canonical_data( Y ~ Z | cid | sid, data=dd,
                                            control_formula = ~ X1 + X2 + X3 + X4 )
    head( data )


    calc_covariate_R2s( data, pooled=TRUE )
    calc_covariate_R2s( data, pooled=FALSE )

}
