


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
    sim.data <- PUMP::gen_sim_data( d_m = design, params, Tbar = 0.30, include_POs = TRUE )
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
        sim.data3 <- PUMP::gen_sim_data( d_m = design, params, Tbar = 0.70, include_POs = TRUE )

        sim.data3$D.id = sim.data.small$D.id[[1]]
        sim.data3$S.id = paste0( "h-", sim.data3$S.id )
        sim.data = bind_rows( sim.data, sim.data3 )

        sim.data <- filter( sim.data, !( S.id %in% c( paste0( "b-", 11:16 ) ) ) )

        # browser()
        if ( FALSE ) {
            ss <- sim.data %>% group_by( D.id, S.id ) %>%
                summarise( n = n(), .groups = "drop" )
            table( ss$n, ss$D.id )
            sim.data %>% group_by( S.id, D.id, T.x ) %>%
                summarise( n = n(), .groups = "drop",
                           Tx = mean( T.x ) ) %>%
                group_by( D.id ) %>%
                summarise( p_tx = mean( Tx ),
                           sizes = paste0( unique( n ), collapse=", " ),
                           nT = sum( n ),
                           J = n(), .groups = "drop" )
        }

        # Rerandomize units within the blocks
        sim.data$T.x = randomizr::block_and_cluster_ra( blocks = sim.data$D.id,
                                                        clusters = sim.data$S.id,
                                                        block_m = c( 10, 6, 10, 5 ) )
        sim.data <- mutate( sim.data,
                            Yobs = ifelse( T.x == 1, Y1, Y0 ) )
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

    describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=dat )
    compare_methods(  Yobs ~ T.x | S.id | D.id, data=dat ) %>%
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
                   Z = mean( T.x ) )
    table( ad$D.id, ad$n )
    sd( dat$Y0 )

    describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=dat )
    compare_methods(  Yobs ~ T.x | S.id | D.id, data=dat ) %>%
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






one_run <- function( blocks = TRUE, het_block = FALSE, ICC.2 = NULL, include_covariates = FALSE,
                     model.params = model.params.list ) {

    p = model.params
    if ( !is.null( ICC.2 ) ) {
        p$ICC.2 = ICC.2
    }

    sim.data = make_two_tier_data( n1 = 10, ATE1 = 0,
                                   n2 = 40, ATE2 = 0.4,
                                   params = p,
                                   blocks = blocks, het_block = het_block )

    form = Yobs ~ T.x | S.id | D.id
    if ( !blocks ) {
        form = Yobs ~ T.x | S.id
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
                                       include_agg = TRUE,
                                       include_dumb = TRUE,
                                       include_method_characteristics = FALSE )

    c1$tau_p =  mean( sim.data$Y1 - sim.data$Y0 )
    sim.data$D.id = "sing"
    c1$tau_c = sim.data %>%
        group_by( D.id, S.id ) %>%
        summarise( tauk = mean( Y1 - Y0 ), .groups = "drop" ) %>%
        pull( tauk ) %>% mean()

    c1$sdY0 = sd( sim.data$Y0 )

    c1 %>% dplyr::select( -any_of( c( "weight", "df" ) ) )
}



if ( FALSE ) {
    # testing
    one_run()

    one_run( blocks = TRUE, het_block = TRUE )


    sim.data = make_two_tier_data( n1 = 10, ATE1 = 0,
                                   n2 = 40, ATE2 = 0.4,
                                   params = model.params.list,
                                   blocks = TRUE, het_block = TRUE )

    c1 <- clusterRCT::compare_methods( form = Yobs ~ T.x | S.id | D.id,
                                       data=sim.data,
                                       include_DB = TRUE,
                                       include_MLM = TRUE,
                                       include_LM = TRUE,
                                       include_agg = TRUE,
                                       include_method_characteristics = FALSE,
                                       include_dumb = TRUE )

    c1
}


run_simple_simulation <- function( ICC.2 = 0.4, omega.2 = 1, R = 1000 ) {

    model.params.list$ICC.2 = ICC.2
    model.params.list$omega.2 = omega.2

    if ( TRUE ) {
        library( furrr )
        plan(multisession, workers = parallel::detectCores() - 1 )
        tictoc::tic()
        rps = future_map_dfr( 1:R, ~ one_run( blocks = TRUE,
                                              het_block = TRUE,
                                              model.params = model.params.list ),
                              .id = "runID", .progress=TRUE, .options = furrr_options( seed=TRUE ) )
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
                   SE = sd( ATE_hat ),
                   SE_E = sd( ATE_hat ) / sqrt(n()),
                   ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
        mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )

    res
}

#### Simple simulation ####

if ( FALSE ) {

    res <- run_simple_simulation( ICC.2 = 0.4, omega.2 = 1, R = 250 )
    res2 <- run_simple_simulation( ICC.2 = 0, omega.2 = 0, R = 250 )

    res <- bind_rows( varied=res, not_varied=res2, .id="sim" )

    res %>%
        arrange( EATE ) %>%
        knitr::kable(digits=2)

    res <- mutate( res, method = fct_reorder2( method, ESE_hat, EATE ) )

    res <- mutate( res, method = factor( method, levels = rev( sort( levels(method) ) ) ) )

    ggplot( res, aes( method, EATE ) ) +
        facet_wrap( ~ sim, nrow=1 ) +
        geom_hline( aes( yintercept = mean( tau_c ),  col="Cluster" ) ) +
        geom_hline( aes( yintercept = mean( tau_p ), col="Person" ) ) +
        geom_point(size=3) +
#        geom_errorbar( aes( ymin = EATE - ESE_hat, ymax = EATE + ESE_hat ), width=0) +
        geom_errorbar( aes( ymin = EATE - SE_E, ymax = EATE + SE_E ), linewidth = 1, width=0) +
        theme_minimal() +
        coord_flip()

}


#### Looking at what estimators align with other estimators ####


if ( FALSE ) {

    R = 10

    #rps = map_df( 1:R, ~ one_run( blocks = TRUE, het_block = TRUE ),
    #              .id = "runID", .progress=TRUE )

    rps = map_df( 1:R, ~ one_run( blocks = TRUE, het_block = FALSE,
                                  include_covariates = TRUE ),
                  .id = "runID", .progress=TRUE )
    head( rps )

    groups = group_estimators( rps )
    groups

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
        map_dfr( run_sim, R = 400, .id = "ICC" )

    res
    res$meth_group = str_extract( res$method, "^[^_]+" )
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
        m <- lmer( Yobs ~ 1 + T.x + (1|S.id), data=dat )
        # arm::display(m)
        ICC1 = as.numeric( VarCorr(m)$S.id / (VarCorr(m)$S.id + sigma(m)^2 ) )

        p = model.params.list
        p$ICC.2 = 0
        dat = make_two_tier_data( 10, 0.2, params=p, blocks = FALSE )
        m2 <- lmer( Yobs ~ 1 + T.x + (1|S.id), data=dat )
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
        summarise( mean_tx = mean( T.x ),
                   SE = sd( T.x ),
                   SE_Etx = sd( T.x ) / sqrt(n()),
                   mean_ICC = mean(ICC),
                   tau_p = mean(tau),
                   sd_t = sd(tau) )
}



