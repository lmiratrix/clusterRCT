


# Simple simulation as sanity check of SEs etc.
#
# This sim has two sizes of clusters, with big clusters having larger
# treatment impacts.
#
# Inspired by Mike Weiss questions on MLM
#


set.seed( 1039 )
library( tidyverse )
library( clusterRCT )


ATE = 0.2
tx_het = FALSE

model.params.list <- list(
    M = 1                            # number of outcomes
    , J = 10                          # number of schools
    , K = 4                          # number of districts
    , nbar = 10                       # number of individuals per school
    , S.id = NULL                     # N-length vector of school assignments
    , D.id = NULL                     # N-length vector of district assignments
    , Xi0 = 0                         # scalar grand mean outcome under no treatment
    , MDES = ATE            # minimum detectable effect size
    , R2.3 = 0.1              # percent of district variation
    , ICC.3 = 0.2             # district intraclass correlation
    , omega.3 = 0.2*tx_het           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.0*tx_het          # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)



#### Functions for the simulations in this file ####


#' @param blocks TRUE means multiple districts (blocks)
#' @param het_block TRUE means blocks are different shapes and styles
#'
make_two_tier_data <- function( n1, ATE1, n2 = 5*n1, ATE2 = 3*ATE1,
                                blocks = TRUE, params = model.params.list,
                                het_block = FALSE ) {

    design = "d3.2_m3ff2rc"
    if ( het_block && !blocks ) {
        stop( "Invalid--can't have het blocks with no blocks" )
    }

    if ( !blocks ) {
        params$ICC.3 = 0
        params$R2.3 = 0
        params$omega.3 = 0
        design = "d2.2"
        params$J = params$J * params$K
        params$K = NULL
    }

    # Make small blocks across the districts
    params$nbar = n1
    params$MDES = ATE1
    sim.data.small <- PUMP::gen_sim_data( d_m = design, params, Tbar = 0.30, include_POs = TRUE )

    # Make large blocks across the districts
    params$nbar = n2
    params$MDES = ATE2
    sim.data <- PUMP::gen_sim_data( d_m = design, params, Tbar = 0.30,
                                    include_POs = TRUE )
    sim.data$S.id = paste0( "b-", sim.data$S.id )

    # Bind so each district has small and large blocks
    sim.data = bind_rows( sim.data.small, sim.data )

    if ( het_block ) {
        # Add a third size of blocks in District 1, and also change
        # the treatment proportion
        params$K = 1
        params$J = round( 1.5*params$J )
        params$MDES = ATE1 + ATE2
        params$nbar = n1 + 2
        sim.data3 <- PUMP::gen_sim_data( d_m = design, params,
                                         Tbar = 0.70, include_POs = TRUE )

        sim.data3$D.id = sim.data.small$D.id[[1]]
        sim.data3$S.id = paste0( "h-", sim.data3$S.id )
        sim.data = bind_rows( sim.data, sim.data3 )

        sim.data <- filter( sim.data, !( S.id %in% c( paste0( "b-", 11:16 ) ) ) )

        # browser()
        if ( FALSE ) {
            ss <- sim.data %>% group_by( D.id, S.id ) %>%
                summarise( n = n(), .groups = "drop" )
            table( ss$n, ss$D.id )
            sim.data %>% group_by( S.id, D.id, Tx ) %>%
                summarise( n = n(), .groups = "drop",
                           Tx = mean( Tx ) ) %>%
                group_by( D.id ) %>%
                summarise( p_tx = mean( Tx ),
                           sizes = paste0( unique( n ), collapse=", " ),
                           nT = sum( n ),
                           J = n(), .groups = "drop" )
        }

        # Rerandomize units within the blocks
        sim.data$Tx = randomizr::block_and_cluster_ra( blocks = sim.data$D.id,
                                                       clusters = sim.data$S.id,
                                                       block_m = c( 10, 6, 10, 5 ) )
        sim.data <- mutate( sim.data,
                            Yobs = ifelse( Tx == 1, Y1, Y0 ) )
    }

    sdY0 = sd( sim.data$Y0 )

    # Standardize
    sim.data <- mutate( sim.data,
                        Yobs = Yobs / sdY0,
                        Y0 = Y0 / sdY0,
                        Y1 = Y1 / sdY0 )
    as_tibble( sim.data )
}





if ( FALSE ) {
    dat = make_two_tier_data( 10, 0.2, het_block = TRUE)
    slice_sample( dat, n=10 )
    table( dat$D.id )
    table( table( dat$S.id ) )
    sd( dat$Y0 )

    describe_clusterRCT( Yobs ~ Tx | S.id | D.id, data=dat )
    compare_methods(  Yobs ~ Tx | S.id | D.id, data=dat ) %>%
        arrange( ATE_hat )

    # Person average ATE
    mean( dat$Y1 - dat$Y0 )
    table( round( dat$Y1 - dat$Y0, digits=5 ) )
    dat %>%
        group_by( D.id, S.id ) %>%
        summarise( tauk = mean( Y1 - Y0 ) ) %>%
        pull( tauk ) %>% mean()



    # More varied data to tease apart methods
    dat = make_two_tier_data( 10, 0.2, het_block = TRUE)
    dat
    table( dat$D.id )
    table( table( dat$S.id ) )

    ad <- dat %>% group_by( D.id, S.id ) %>%
        summarise( n = n(),
                   Z = mean( Tx ) )
    table( ad$D.id, ad$n )
    sd( dat$Y0 )

    describe_clusterRCT( Yobs ~ Tx | S.id | D.id, data=dat )
    compare_methods(  Yobs ~ Tx | S.id | D.id, data=dat ) %>%
        arrange( ATE_hat )




    # Person average ATE
    mean( dat$Y1 - dat$Y0 )
    table( round( dat$Y1 - dat$Y0, digits=5 ) )
    dat %>%
        group_by( D.id, S.id ) %>%
        summarise( tauk = mean( Y1 - Y0 ) ) %>%
        pull( tauk ) %>% mean()




    dat = make_two_tier_data( 10, 0.2, blocks=FALSE )
    dat
    dat$D.id
    table( table( dat$S.id ) )
    sd( dat$Y0 )
}





#' Make data with very different block structure and cluster mixtures
#' that correlate with treatment effects.
#'
#' Each block will be a mix of clusters of specified sizes, with
#' specified number of clusters of each size.  The clusters will each
#' have a mean outcome for that cluster, and a cluster-specific
#' treatment effect.  I.e., blocks are a mixture of a discrete number
#' of cluster types.
#'
#' "Randomizations" will be perfectly balanced to focus on bias, not
#' estimation error.
#'
#' @param n List lists of cluster sizes within each block
#'
#'
make_block_data <- function( n, J, Y0 = 0, ATE, ptx = 0.5 ) {

    statblock = NA
    if ( is.data.frame( n ) ) {
        statblock = n

        if ( !( "ptx" %in% colnames( statblock ) ) ) {
            statblock$ptx = 0.5
        }
    } else {

        if ( length(ptx) == 1 ) {
            ptx = rep( ptx, length(J) )
        }
        if ( length(Y0) == 1 ) {
            Y0 = rep( Y0, length(J) )
        }

        stopifnot( length(n) == length(J) )
        stopifnot( length(n) == length(Y0) )
        stopifnot( length(n) == length(ATE) )
        stopifnot( length(n) == length(ptx) )

        statblock = list( n, J, Y0, ATE, ptx )

    }

    # Make J clusters of size n with an ATE of given ATE.
    make_blocks <- function( n, J, Y0, ATE, ptx = 0.5, prefix ) {
        dat = tibble( S.id = rep( 1:J, each=n ),
                      Y0 = Y0 + rnorm( n*J, sd=1 ) )
        dat$Y1 = dat$Y0 + ATE
        dat$Tx = 0 + (dat$S.id %% (1/ptx) == 0)
        dat$Yobs = ifelse( dat$Tx, dat$Y1, dat$Y0 )
        dat$S.id = paste0( prefix, "-", dat$S.id )
        dat
    }

    # Given vectors of n, J, Y0, and ATE, make batches of clusters of
    # given sizes via pmap
    make_blocks_v <- function( n, J, Y0, ATE, ptx=0.5 ) {
        stopifnot( length(n) == length(J) )
        stopifnot( length(n) == length(Y0) )
        stopifnot( length(n) == length(ATE) )
        if ( length(ptx) == 1 ) {
            ptx = rep( ptx, length(J) )
        }
        stopifnot( length(n) == length(ptx) )
        if ( length( Y0 ) == 1 ) {
            Y0 = rep( Y0, length(n) )
        }

        res <- pmap_df( list( n=n, J=J, Y0=Y0, ATE=ATE, ptx=ptx,
                              prefix=letters[ 1:length(n) ] ),
                        make_blocks )
        res
    }

    dats = pmap( statblock, make_blocks_v )
    names(dats) = LETTERS[ 1:length(dats) ]
    dat = bind_rows( dats, .id="D.id" )

    # Standardize
    sdY0 = sd( dat$Y0 )
    sim.data <- mutate( dat,
                        Yobs = Yobs / sdY0,
                        Y0 = Y0 / sdY0,
                        Y1 = Y1 / sdY0 )

    # Add some covariates
    sim.data$X.jk = sim.data$Y0 + rnorm( nrow(sim.data), sd=1)
    sim.data <- sim.data %>%
        group_by( S.id ) %>%
        mutate( X.jk = mean( X.jk ) ) %>%
        ungroup()
    sim.data$C.ijk = sim.data$Y0 + rnorm( nrow(sim.data), sd=1)

    as_tibble( sim.data )
}



if ( FALSE ) {

    dd <- tibble(
        n = list( c( 200, 400 ), 300 ),
        J = list( c( 4, 2 ), 2 ),
        Y0 = list( c( 0, 0 ), 0 ),
        ATE = list( c( 1, 2 ), 3 ),
        ptx = c( 0.5, 0.5 ) )
    dd

    dd %>%
        make_block_data( ) %>%
        group_by( D.id, S.id ) %>%
        summarise( Tx = mean( Tx ),
                   Yobs = mean( Yobs ),
                   n = n(), .groups = "drop" ) %>%
        print( n = 100 )
}







#' Make data with very different block structure and cluster mixtures
#' that correlate with treatment effects.
#'
#' "Randomizations" will be perfectly balanced to focus on bias, not
#' estimation error.
#'
#' @param blocks TRUE means multiple districts (blocks)
#' @param het_block TRUE means blocks are different shapes and styles
#'
make_simple_block_example <- function( version = "A", tx_scale = 1, var_tx = FALSE ) {

    ptx = rep( 0.5, 3 )
    if ( var_tx ) {
        ptx = c( 1/2, 1/6, 1/2 )
    }

    if ( version == "A" ) {
        blks <- tribble( ~n, ~ J, ~ Y0, ~ ATE,
                         c( 20, 80 ), c( 2, 4 ), c( 1, 0.5 ), c( 1.4, 0.2 ),

                         c( 20, 80 ), c( 6, 12 ), c( 0, 1.5 ), c( 1.2, 2.4 ),
                         c( 5, 100 ), c( 20, 4 ), c( 4, 0 ), c(  4, -4 ) )
    } else {
        blks <- tribble( ~n, ~ J, ~ Y0, ~ ATE,
                         c( 20, 80 ), c( 2, 4 ), c( 1, 0 ),  c( 0, -2 ),
                         c( 10, 30 ), c( 6, 6 ), c( 0, 2 ), c( 3, 1 ),
                         c( 5, 10 ), c( 10, 10 ), c( 2, 0 ), c(  4, 2 ) )
    }
    blks$ptx = ptx
    blks$ATE = map( blks$ATE, ~ .x * tx_scale )

    make_block_data( blks )
}








##### Test make_simple_block_example() #####

if ( FALSE ) {
    dat = make_simple_block_example( version = "B", var_tx = TRUE )
    dat

    ds = dat %>%
        group_by( D.id, S.id, Tx ) %>%
        summarise( tauk = mean( Y1 - Y0 ),
                   n = n(), .groups = "drop" )

    ff <- filter( ds, D.id == "A" )
    ff
    mean( ff$tauk )
    weighted.mean( ff$tauk , ff$n )

    # Make table of block-average and person-average impacts by block
    ds <- ds %>%
        group_by( D.id ) %>%
        summarise( tau_c = mean( tauk ),
                   tau_p = weighted.mean( x=tauk, w=n ),
                   #       tau_pp = mean( tauk * n ) / mean(n),
                   p_tx = mean( Tx ),
                   J = n(),
                   Jwt = J * p_tx * (1-p_tx),
                   n = sum( n ),
                   nwt = n * p_tx * (1-p_tx ) )
    ds %>% knitr::kable( digits = 3 )

    ds %>%
        summarise( tau_cb = mean( tau_c ),
                   tau_pb = mean( tau_p ),
                   tau_p_FE = weighted.mean( tau_p, w = n*p_tx*(1-p_tx) ),
                   tau_c_FE = weighted.mean( tau_c, w = n*p_tx*(1-p_tx) ),
                   tau_c = weighted.mean( tau_c, w=J ),
                   tau_p = weighted.mean( tau_p, n ),

                   J = sum( J ),
                   n = sum( n ) )
}






one_run <- function( blocks = TRUE, het_block = FALSE, simple = FALSE, ICC.2 = NULL,
                     include_covariates = FALSE,
                     model.params = model.params.list ) {

    p = model.params
    if ( !is.null( ICC.2 ) ) {
        p$ICC.2 = ICC.2
    }

    if ( simple != FALSE ) {
        sim.data = make_simple_block_example( version = simple, var_tx = TRUE )
    } else {
        sim.data = make_two_tier_data( n1 = 10, ATE1 = 0,
                                       n2 = 40, ATE2 = 0.4,
                                       params = p,
                                       blocks = blocks, het_block = het_block )
    }

    form = Yobs ~ Tx | S.id | D.id
    if ( !blocks ) {
        form = Yobs ~ Tx | S.id
    }
    control_form = NULL
    if ( include_covariates ) {
        if ( blocks ) {
            control_form = ~ X.jk + C.ijk
        } else {
            control_form = ~ C.ijk
        }
    }

    c1 <- clusterRCT::compare_methods( form, data=sim.data,
                                       control_formula = control_form,
                                       include_DB = TRUE,
                                       include_MLM = TRUE,
                                       include_LM = TRUE,
                                       include_agg = FALSE,
                                       include_dumb = FALSE,
                                       include_method_characteristics = FALSE )

    c2 <- clusterRCT::canonical_weight_estimators( form, sim.data )

    c1 = bind_rows( c1, c2 )

    if ( !blocks ) {
        sim.data$D.id = "sing"
    }

    ds = sim.data %>%
        group_by( D.id, S.id ) %>%
        summarise( tauk = mean( Y1 - Y0 ),
                   n = n(), .groups = "drop" ) %>%
        group_by( D.id ) %>%
        summarise( tau_c = mean( tauk ),
                   tau_p = weighted.mean( tauk, w=n ),
                   J = n(),
                   n = sum( n ) ) %>%
        summarise( tau_cb = mean( tau_c ),
                   tau_c = weighted.mean( tau_c, w=J ),
                   tau_pb = mean( tau_p ),
                   tau_p = weighted.mean( tau_p, w=n ),
                   J = sum( J ),
                   n = sum( n ) )

    c1$tau_c = ds$tau_c
    c1$tau_p = ds$tau_p
    c1$tau_cb = ds$tau_cb
    c1$tau_pb = ds$tau_pb
    c1$J = ds$J
    c1$n = ds$n
    c1$sdY0 = sd( sim.data$Y0 )

    c1 %>% dplyr::select( -any_of( c( "weight", "df" ) ) )
}



if ( FALSE ) {
    # testing
    one_run( simple=TRUE )

    one_run( blocks = TRUE, het_block = TRUE )


    sim.data = make_two_tier_data( n1 = 10, ATE1 = 0,
                                   n2 = 40, ATE2 = 0.4,
                                   params = model.params.list,
                                   blocks = TRUE, het_block = TRUE )

    c1 <- clusterRCT::compare_methods( form = Yobs ~ Tx | S.id | D.id,
                                       data=sim.data,
                                       include_DB = TRUE,
                                       include_MLM = TRUE,
                                       include_LM = TRUE,
                                       include_agg = TRUE,
                                       include_method_characteristics = FALSE,
                                       include_dumb = TRUE )

    c1
}


run_simple_simulation <- function( ICC.2 = 0.4, omega.2 = 1,
                                   simple = FALSE, R = 1000 ) {

    model.params.list$ICC.2 = ICC.2
    model.params.list$omega.2 = omega.2

    if ( TRUE ) {
        library( furrr )
        plan(multisession, workers = parallel::detectCores() - 1 )
        tictoc::tic()
        rps = future_map_dfr( 1:R, ~ one_run( blocks = TRUE,
                                              het_block = TRUE,
                                              simple = simple,
                                              model.params = model.params.list ),
                              .id = "runID", .progress=TRUE,
                              .options = furrr_options( seed=TRUE ) )
        tictoc::tic()
    } else {
        rps = map_df( 1:R, ~ one_run( blocks = TRUE,
                                      het_block = TRUE,
                                      model.params = model.params.list ),
                      .id = "runID", .progress=TRUE )
    }


    head( rps )
    sd( rps$tau_c)

    res <- rps %>%
        group_by( method ) %>%
        summarise( EATE = mean( ATE_hat ),
                   tau_p = mean( tau_p ),
                   tau_c = mean( tau_c ),
                   tau_cb = mean( tau_cb ),
                   tau_pb = mean( tau_pb ),
                   SE = sd( ATE_hat ),
                   SE_E = sd( ATE_hat ) / sqrt(n()),
                   ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
        mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )

    res
}


if ( FALSE) {

    debug( one_run )
    one_run( simple = TRUE )
}


#### Simple simulation ####

# This simulation generates data and compares the estimators to the 4
# different estimands to see which are biased.
#
# The third DGP ("simple") is a special DGP that hand-builds blocks of
# different sizes to ensure that there are correlations with impact
# and size of block and number of blocks in a district, and so forth.
#
# The other DGPs are more standard two-tier DGPs with different
# intraclass correlations and effect size variability.

if ( TRUE ) {

    R = 100

    #    res1 <- run_simple_simulation( ICC.2 = 0.4, omega.2 = 1, R = R )
    #    res2 <- run_simple_simulation( ICC.2 = 0, omega.2 = 0, R = R )


    res3 <- run_simple_simulation( simple = "B", R = R )
    res <- res3 # hack to skip res1 and res2
    res$sim = "simple"

    #res <- bind_rows( varied=res1, not_varied=res2, simple=res3, .id="sim" )

    res %>%
        arrange( EATE ) %>%
        knitr::kable(digits=2)

    res$ESE_hat[ is.na( res$ESE_hat ) ] = 0
    res <- mutate( res, method = fct_reorder2( method, ESE_hat, EATE ) )
    res
    #res <- mutate( res, method = factor( method,
    #                                     levels = rev( sort( levels(method) ) ) ) )

    # Look at third simulation only
    res %>%
        filter( sim == "simple" ) %>%
        ggplot( aes( method, EATE ) ) +
        geom_hline( aes( yintercept = mean( tau_c ),  lty="clean", col="Cluster" ) ) +
        geom_hline( aes( yintercept = mean( tau_p ), lty="clean", col="Person" ) ) +
        geom_hline( aes( yintercept = mean( tau_cb ),  lty = "block", col="Cluster" ) ) +
        geom_hline( aes( yintercept = mean( tau_pb ), lty="block", col="Person" ) ) +
        geom_point(size=3) +
        #        geom_errorbar( aes( ymin = EATE - ESE_hat, ymax = EATE + ESE_hat ), width=0) +
        geom_errorbar( aes( ymin = EATE - SE_E, ymax = EATE + SE_E ), linewidth = 1, width=0) +
        theme_minimal() +
        coord_flip()


}


#### Looking at what estimators align with other estimators ####

# This runs a few iterations
if ( FALSE ) {

    R = 20

    #rps = map_df( 1:R, ~ one_run( blocks = TRUE, het_block = TRUE ),
    #              .id = "runID", .progress=TRUE )

    rps = map_df( 1:R, ~ one_run( blocks = TRUE, het_block = FALSE, simple = TRUE,
                                  include_covariates = TRUE ),
                  .id = "runID", .progress=TRUE )
    head( rps )

    source( here::here( "simulations/simulation_helper_functions.R" ) )
    groups = group_estimators( rps )
    groups


    groups = group_estimators( rps, tolerance = 0.05, return_matrix = TRUE )
    summary( as.numeric( groups ) )

    # List the groups
    map_chr( groups, paste, collapse = ", " ) %>%
        as.data.frame() %>%
        mutate( grp = 1:n() ) %>%
        relocate( grp ) %>%
        knitr::kable()

    if ( FALSE ) {
        # Make list of estimators with group IDs (optional)
        group_df <- map_df(seq_along(groups), ~tibble(estimator = groups[[.x]], group = .x))

        print(group_df)
        table( group_df$group )

    }
}




#### Mike's ICC Simulation ####

# Mike wanted to see how the MLM estimators targeted different
# estimands based on the ICC

if ( FALSE ) {

    run_sim <- function( ICC, R = 100 ) {
        cat( "ICC ", ICC, "\n" )

        rps = map_df( 1:R, ~ one_run( ICC.2 = ICC, blocks = FALSE ),
                      .id = "runID", .progress=TRUE )

        head( rps )


        res <- rps %>%
            group_by( method ) %>%
            summarise( EATE = mean( ATE_hat ),
                       tau_p = mean( tau_p ),
                       tau_c = mean( tau_c ),
                       SE = sd( ATE_hat ),
                       ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
            mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )


        res
    }

    #one_run( blocks=FALSE, ICC.2 = 0.8 )

    ICCs = c( 0, 0.2, 0.4, 0.6, 0.8 )
    res <- ICCs %>%
        set_names() %>%
        map_dfr( run_sim, R = 100, .id = "ICC" )

    res

    # Make families of estimators for grouping on the plot
    res$meth_group = str_extract( res$method, "^[^_]+" )

    # Note how all methods are level w.r.t. the ICC, except MLM which
    # goes from person-weighted to cluster-weighted as ICC increases.
    ggplot( res, aes( ICC, EATE ) ) +
        facet_wrap( ~ meth_group, nrow=1 ) +
        scale_y_continuous( limits = c(0, 0.5 ) ) +
        geom_line( aes( col = method, group=method )) +
        geom_hline( aes( yintercept = mean( tau_c ), lty="Cluster" ), col="black" ) +
        geom_hline( aes( yintercept = mean( tau_p ), lty="Person" ), col="black" )


}


# Checking ICC is working
if ( FALSE ) {

    rr = rerun( 100, {

        p = model.params.list
        p$ICC.2 = 0.8
        dat = make_two_tier_data( 10, 0.2, params=p, blocks = FALSE )
        m <- lmer( Yobs ~ 1 + Tx + (1|S.id), data=dat )
        # arm::display(m)
        ICC1 = as.numeric( VarCorr(m)$S.id / (VarCorr(m)$S.id + sigma(m)^2 ) )

        p = model.params.list
        p$ICC.2 = 0
        dat = make_two_tier_data( 10, 0.2, params=p, blocks = FALSE )
        m2 <- lmer( Yobs ~ 1 + Tx + (1|S.id), data=dat )
        # arm::display(m2)
        ICC2 = as.numeric( VarCorr(m2)$S.id / (VarCorr(m2)$S.id + sigma(m2)^2 ) )

        res <- bind_rows( `0.8`=fixef( m ), `0`=fixef( m2 ), .id = "mod" )
        res$ICC = c( ICC1, ICC2 )
        res$tau = mean(dat$Y1 - dat$Y0 )
        res
    } )
    rr = bind_rows( rr )
    rr
    rr %>% group_by( mod ) %>%
        summarise( mean_tx = mean( Tx ),
                   SE = sd( Tx ),
                   SE_Etx = sd( Tx ) / sqrt(n()),
                   mean_ICC = mean(ICC),
                   tau_p = mean(tau),
                   sd_t = sd(tau) )
}



