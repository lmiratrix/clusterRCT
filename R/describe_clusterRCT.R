


#' Describe a cluster RCT dataset
#'
#' Calculate various statistics describing total clusters, etc.  If
#' formula is NULL, data assumed to be in canonical form.
#'
#' Will patch data via patch_data() to handle missing values.
#'
#' @inheritParams compare_methods
#' @seealso [patch_data_set()]
#'
#' @export
describe_clusterRCT <- function( formula = NULL,
                                 data = NULL,
                                 control_formula = NULL,
                                 warn_missing = TRUE ) {

    if ( !is.null( formula ) ) {
        # Deal with multiple outcomes
        lhs_vars <- formula.tools::lhs.vars(formula)
        if (length(lhs_vars) > 1) {
            res <- purrr::map(
                lhs_vars,
                function(v) {
                    # see https://stackoverflow.com/questions/37824625/how-to-modify-the-left-side-of-a-formula
                    outcome <- as.name(v)
                    formula1 <- as.formula(bquote(.(outcome) ~ .(formula.tools::rhs(formula))))
                    res <- describe_clusterRCT(
                        formula1,
                        data = data,
                        control_formula = control_formula ) %>%
                        dplyr::mutate(outcome = v, .before=n)
                }) %>%
                purrr::list_rbind()
            return(res)
        } else {
            data = make_canonical_data( formula=formula, data=data,
                                        control_formula = control_formula,
                                        drop_missing = FALSE,
                                        warn_missing = FALSE )
        }
    }

    cnames = NULL
    if ( !is.null( control_formula ) ) {
        cnames = setdiff( colnames(data), c("Yobs", "Z", "clusterID", "blockID" ) )
    }

    # count up missing data and report, then patch dataset.
    miss = colSums(is.na(data))

    data = patch_data_set( NULL, data=data, control_formula = control_formula,
                           warn_missing = warn_missing )
    data = patch_singleton_blocks(NULL, data=data,
                                  warn_missing = warn_missing )

    control_formula = attr( data, "control_formula" )

    K = length( unique( data$blockID ) )

    # Get cluster sizes
    sizes = data %>%
        group_by( blockID, clusterID, Z ) %>%
        summarise( n = n(), .groups = "drop" )

    cstat = sizes %>%
        summarise( J = n(),
                   nbar = mean(n),
                   ncv = sd(n) / nbar,
                   n.25 = quantile( n, 0.25 ),
                   n.75 = quantile( n, 0.75 ),
                   n.IQR = n.75 - n.25,
                   p.tx = mean(Z) )

    sstat = sizes %>%
        group_by( blockID ) %>%
        summarise( J = n(),
                   n = sum(n),
                   p.tx = mean(Z) ) %>%
        summarise( K = n(),
                   Jbar = mean(J),
                   Jcv = sd(J) / Jbar,
                   J.25 = quantile( J, 0.25 ),
                   J.75 = quantile( J, 0.75 ),
                   J.IQR = J.75 - J.25,
                   n_block = mean(n),
                   n_block_cv = sd( n ) / n_block,
                   tx.avg = mean( p.tx ),
                   tx.cv = sd( p.tx ) / tx.avg,
                   tx.25 = quantile( p.tx, 0.25 ),
                   tx.75 = quantile( p.tx, 0.75 ),
                   tx.IQR = tx.75 - tx.25 )


    ICCs = calc_ICCs(data)

    stats <- tibble( n = nrow(data) )
    stats = bind_cols( stats, cstat, sstat )
    stats$cluster_ICC = ICCs$C.ICC
    stats$block_ICC = ICCs$S.ICC
    stats$sdY0 = sd( data$Yobs[ data$Z == 0] )

    # Calculate some notes about block distribution.
    notes = c()
    tb = table( sizes$blockID, sizes$Z )

    if ( any( tb == 1 ) ) {
        if ( all( tb <= 1 ) ) {
            notes = c( notes, "matched pairs design" )
        } else {
            nblk = apply( tb == 1, 1, max )
            notes = c( notes, glue::glue( "{sum(nblk)} blocks have at least one singleton treated or control clusters" ) )
            nblk = apply( tb == 1, 1, min )
            notes = c( notes, glue::glue( "{sum(nblk)} blocks have exactly 2 clusters, one treated and one control" ) )
        }
        stats$has_singletons = TRUE
    } else {
        stats$has_singletons = FALSE
    }

    if ( any( tb == 2 ) ) {
        notes = c( notes, glue::glue( "{sum(tb==2)} within-block treatment arms have exactly 2 clusters." ) )
    }
    attr( stats, "notes" ) <- notes

    class( stats ) <- c( "clusterRCTstats", class( stats ) )

    # R2 values
    if ( !is.null( cnames ) ) {
        stats <- bind_cols( stats, calc_covariate_R2s(data, ICCs) )
    }

    if ( sum( miss ) > 0 ) {
        stats$.missing = list( miss )
    }

    stats
}



#' Calculate ICCs
#'
#' Use a multilevel model to calculate ICCs
#'
calc_ICCs <- function( data ) {

    require( lme4 )

    form = Yobs ~ 1 + Z + (1 | blockID ) + (1 | clusterID )
    #if ( !is.null( control_formula ) ) {
    #    form = update( control_formula, Yobs ~ . + 1 + Z + (1 | blockID ) + (1 | clusterID ) )
    #}
    M = lmer( form, data=data )
    #arm::display(M)

    a = VarCorr( M )

    C.ICC = as.numeric( a$clusterID )
    S.ICC = as.numeric( a$blockID )
    sigma2 = sigma( M )^2
    tvar = (C.ICC + S.ICC + sigma2)
    C.ICC = C.ICC / tvar
    S.ICC = S.ICC / tvar

    list( C.ICC = C.ICC, S.ICC = S.ICC )
}


#' Calculate R2 of level 1, level 2, and block.
#'
#' block is always the fixed effects for block.  No covariates
#' allowed at that level.  In all models block fixed effects are
#' included, so level 2 covarites R2 are additional explanatory power
#' beyond their representing block differences (e.g., they are
#' considered block centered).
#'
#' @importFrom purrr map_dbl
#' @param data in canonical form.
#' @param ICCs List of ICC values, level 2 and level 3 in order.
#' @param pooled Calculate pooled Tx and Co groups, attempting to
#'   remove systematic shift of average treatment effect.  If FALSE,
#'   subset to control observations only.
#'
calc_covariate_R2s <- function( data, ICCs = NULL, pooled = FALSE ) {

    if ( is.null( ICCs ) ) {
        ICCs = calc_ICCs(data)
    }

    # Covert covariates to numerical, if they are not already.
    dm = model.matrix( Yobs ~ . - clusterID - blockID, data=data ) %>%
        as.data.frame() %>%
        dplyr::select( -`(Intercept)` )

    cnames = colnames(dm)

    dm$clusterID = data$clusterID
    dm$blockID = data$blockID

    # Add group means and group mean-centered versions of all
    # covariates
    result <- dm %>%
        group_by(blockID, clusterID) %>%
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
    R2s = NA
    if ( pooled ) {
        # This way is not preferred--I forget why.
        Mblock = lm( Yobs ~ Z_mn*blockID, data=result )
        R2.block = summary(Mblock)$adj.r.squared

        cents = result %>% dplyr::select( !ends_with( "_mn" ) )
        cents$Z_mn = result$Z_mn
        M = lm( Yobs ~ . + Z_mn*blockID - clusterID, data=cents )
        R2.cent = summary(M)$adj.r.squared

        mns = result %>% dplyr::select( !ends_with("_cent" ) )
        Mmn = lm( Yobs ~ . + Z_mn*blockID - clusterID, data=mns )
        R2.mn = summary( Mmn )$adj.r.squared

        Mtx = lm( Yobs ~ Z_mn, data=result )
        sTx <- summary(Mtx)
        tx.R2 = sTx$r.squared
        rat = 1 / (1 - sTx$r.squared)

        R2s <- tibble( R2.1 = (R2.cent - R2.block)*rat,
                       ncov.1 = ncov.1,
                       R2.2 = (R2.mn - R2.block)*rat,
                       ncov.2 = ncov.2 )
    } else {
        # This is the way we default to!

        # Look at control observations only.
        resCo = dplyr::filter( result, Z_mn == 0 )
        resCo$Z_mn = NULL

        Mblock = lm( Yobs ~ blockID, data=resCo )
        R2.block = summary(Mblock)$adj.r.squared

        cents = resCo %>%
            dplyr::select( !ends_with( "_mn" ) )
        M = lm( Yobs ~ . - clusterID, data=cents )
        summary( M )
        R2.cent = summary(M)$adj.r.squared

        mns = resCo %>% dplyr::select( !ends_with("_cent" ) )
        Mmn = lm( Yobs ~ . - clusterID, data=mns )
        summary( Mmn )
        R2.mn = summary( Mmn )$adj.r.squared

        R2s <- tibble( R2.1 = R2.cent - R2.block,
                       ncov.1 = ncov.1,
                       R2.2 = R2.mn - R2.block,
                       ncov.2 = ncov.2 )
    }

    R2s$R2.1 = R2s$R2.1 / (1 - ICCs[[1]] - ICCs[[2]])
    R2s$R2.2 = R2s$R2.2 / ICCs[[1]]

    R2s
}



#' Make table of statistics by block
#'
#' Given individual/school/block data of a blocked, cluster RCT,
#' calculate statistics for each block (block).
#'
#' @param check_data_integrity TRUE means runs some checks and give
#'   errors if data fails them (e.g., incorrectly processed treatment
#'   vector.). FALSE means calculate statistics without these checks.
#' @export
make_block_table <- function(  formula = NULL,
                               data = NULL,
                               control_formula = NULL,
                               check_data_integrity = FALSE ) {

    data = make_canonical_data( formula=formula, data=data,
                                control_formula = control_formula )

    if ( check_data_integrity ) {
        check_data_integrity( data )
    }

    n <- nrow( data )

    K = length( unique( data$blockID ) )

    sizes = data %>%
        group_by( blockID, clusterID, Z ) %>%
        summarise( n = n(), .groups = "drop" )

    sstat = sizes %>% group_by( blockID ) %>%
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

    if ( nrow( x ) > 1 ) {
        rr <- print.data.frame( x, ... )
        return( rr )
    }

    if (  is.null( x$K ) || x$K == 1 ) {
        scat( "Cluster RCT: %d units in %d clusters\n", x$n, x$J )
    } else {
        scat( "Cluster RCT: %d units in %d clusters across %d blocks\n",
              x$n, x$J, x$K )
    }


    if ( !is.null( x$K ) && x$K != 1 ) {
        scat( "Block Statistics:\n\tAvg units/block: %.2f (coef var %.2f)\n\tAvg clusters/block: %.2f (coef var %.2f)\n\t25-75 Quantiles: %.1f-%.1f w/ IQR = %.1f\n\tICC: %.2f\n\tAvg tx: %.2f (coef var %.2f)\n\t25-75 Quantile tx: %.2f-%.2f w/ IQR %.2f\n",
              x$n_block, x$n_block_cv,
              x$Jbar, x$Jcv,
              x$J.25, x$J.75, x$J.IQR,
              x$block_ICC,
              x$tx.avg, x$tx.cv,
              x$tx.25, x$tx.75, x$tx.IQR)
        if ( hasName( x, "R2.3" ) ) {
            scat( "\tR2.3: %.2f (%d fixed effects)\n", x$R2.3, x$K)
        }
    }

    scat( "Cluster Statistics:\n\tAvg size: %.2f (coef var %.2f)\n\t25-75 Quantiles: %.1f-%.1f w/ IQR = %.1f\n\tICC: %.2f\n\tProp treated: %.2f\n",
          x$nbar, x$ncv,
          x$n.25, x$n.75, x$n.IQR,
          x$cluster_ICC,
          x$p.tx )
    if ( hasName( x, "R2.2" ) ) {
        scat( "\tR2.2: %.2f (%d covariates)\n", x$R2.2, x$ncov.2)
    }


    scat( "Unit Statistics:\n\tstddev( Y0 ): %.2f\n",
          x$sdY0 )
    if ( hasName( x, "R2.1" ) ) {
        scat( "\tR2.1: %.2f (%d covariates)\n", x$R2.1, x$ncov.1)
    }

    if ( ".missing" %in% names(x) ) {
        scat( "missing data counts:\n" )
        print( x$.missing[[1]] )
    }

    notes = attr( x, "notes" )
    if ( !is.null(notes) && ( length( notes ) > 0 ) ) {
        cat( "Notes:\n" )
        for ( n in notes ) {
            cat( "\t", n, "\n" )
        }
    }

    invisible( x )

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
    x$.missing = NULL
    class(x) = "list"
    x <- as.data.frame( x )
    rownames(x) = 1:nrow(x)

    x
}





count_block_sizes <- function( clusterID, blockID ) {
    tt = tibble( clusterID = clusterID, blockID = blockID )
    t2 <- tt %>%
        group_by( blockID ) %>%
        summarize( n = n(),
                   J = length( unique( clusterID ) ) )

    t2
}




#### Testing code ####

if ( FALSE ) {

    library( tidyverse )
    library( clusterRCT )

    set.seed( 1039 )


    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 30                          # number of clusters
        , K = 10                          # number of blocks
        , nbar = 10                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of block assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.1              # percent of block variation
        , ICC.3 = 0.2             # block intraclass correlation
        , omega.3 = 0.2           # ratio of block effect size variability
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
    dd <- describe_clusterRCT( data )
    class( dd )
    dd2 = dd

    b = bind_rows( dd, dd2 )
    class(b)
    as.data.frame(b)

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

    make_block_table(  Y ~ Z | cid | sid, data=dd )
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
