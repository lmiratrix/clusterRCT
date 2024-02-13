
# Test designs with singleton clusters treated, or doubleton, etc.
#
# generally we should always get an ATE, but we might not get a SE,
# unless we are doing covariate adjustment, which can then bring down
# some of the estimators.
#
# Describe data should note these as unusual designs in the printout.


test_that( "matched pairs designs get handled right", {

    blk = blkvar::generate_blocked_data( n_k = c( 5, 10, 20, 10, 5, 10 ) )
    head( blk )
    blks = tibble( B = paste0( "B", 1:6 ),
                   D = rep( 1:3, each=2 ),
                   Z = rep( c(0,1), 3 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    blk$X = rnorm( nrow(blk) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk, control_formula = ~ X )
    dsc
    expect_true( any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )

    expect_warning( expect_warning( comp <- compare_methods( Yobs ~ Z | B | D, data=blk ) ) )

    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )


    head(blk)
    suppressWarnings( comp <- compare_methods( Yobs ~ Z | B | D, data=blk, control_formula = ~ X ) )
    comp
    expect_true( any( is.na( comp$ATE_hat ) ) )

})


test_that( "singleton tx or co designs get handled right", {

    blk = blkvar::generate_blocked_data( n_k = c( 5, 10, 20, 10, 11, 13,
                                                  10, 10, 20, 10,
                                                  30, 10, 5, 10, 30, 15 ) )
    head( blk )
    blks = tibble( B = paste0( "B", 1:16 ),
                   D = rep( 1:3, c(6,4,6) ),
                   Z = randomizr::block_ra(blocks=D, block_m = c(3,1,3) ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    dsc
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )

    comp = compare_methods( Yobs ~ Z | B | D, data=blk )
    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )



    # One block with 2 clusters
    n_k = c( 5,  10, 10, 20, 30, 40,
             20, 10,
             10, 20, 30, 10, 40, 10,
             40, 10, 5, 10, 20, 30 )
    blk = blkvar::generate_blocked_data( n_k = n_k )
    head( blk )
    blks = tibble( B = paste0( "B", 1:length(n_k) ),
                   D = rep( c( 1, 2, 3, 4 ), c(6,2,6,6) ),
                   Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    dsc
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )

    comp2 = compare_methods( Yobs ~ Z | B | D, data=blk )
    comp2
    expect_true( all( !is.na( comp2$ATE_hat ) ) )

    left_join( comp[ c(1:3) ], comp2[ c(1:3) ], by="method" )
    expect_equal( is.na( comp$SE_hat ), is.na( comp2$SE_hat ) )
})




test_that( "doubleton tx or co designs get handled right", {

    blk = blkvar::generate_blocked_data( n_k = c( 5,  10, 10, 20, 30, 40,
                                                  20, 10,  5, 10, 10,
                                                  10, 20, 30, 10, 40, 10 ) )
    head( blk )
    blks = tibble( B = paste0( "B", 1:17 ),
                   D = rep( c( 1, 2, 3 ), c(6,5,6) ),
                   Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    dsc
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )

    comp = compare_methods( Yobs ~ Z | B | D, data=blk )
    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )
    expect_true( all( !is.na( comp$SE_hat ) ) )


    # All blocks have  2 tx and 2 co units in it.
    n_k = c( 5,  10, 30, 40,
             20, 10,  5, 10,
             10, 20, 40, 10 )
    blk = blkvar::generate_blocked_data( n_k = n_k )
    head( blk )
    blks = tibble( B = paste0( "B", 1:length(n_k) ),
                   D = rep( c( 1, 2, 3 ), c(4,4,4) ),
                   Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    dsc
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )

    comp = compare_methods( Yobs ~ Z | B | D, data=blk )
    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )
    expect_true( all( !is.na( comp$SE_hat ) ) )



})

