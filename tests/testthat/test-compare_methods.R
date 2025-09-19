

if ( FALSE ) {
    library( testthat )

}

data( "fakeCRT" )

test_that("compare methods gives multiple estimates as expected", {

    mtab <- compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                             include_method_characteristics = TRUE,
                             include_disfavored = TRUE )
    mtab
    expect_true( is.data.frame(mtab) )
    expect_true( sum(is.na(mtab)) == 0 )
    nrow( mtab )

    mtab


    mc = method_characteristics( include_weight = TRUE )
    mc
    expect_true( is.data.frame(mc) )

    mtab$mm = get_estimand( mtab$method, simple = FALSE )
    expect_true( all( !is.na( mtab$mm ) ) )

    mm = left_join( mtab, mc, by = "method" )
    if ( FALSE) {
        table( mm$weight.x )
        mm %>%
            dplyr::select( method, weight.x, weight.y, mm ) %>%
            filter( weight.x != weight.y | weight.x != mm )
        names(mm)
    }
    expect_true( all( mm$weight.x == mm$weight.y ) )
    expect_true( all( mm$mm == mm$weight.x ) )


    mtab$m2 = get_estimand( mtab$method, simple = TRUE )
    expect_equal( names( table( mtab$m2 ) ), c( "Cluster", "Person" ) )

    mtab_cov <-  compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                  control_formula = ~ X.jk + C.ijk,
                                  include_method_characteristics = FALSE,
                                  include_disfavored = TRUE )

    expect_equal( nrow( mtab ), nrow( mtab_cov ) )


    # Ensure covariate correction changes the estimates!  (I.e.,
    # control formula is getting propagated)
    mtab$dels = mtab$ATE_hat - mtab_cov$ATE_hat
    mtab

    # HT methods ignore covariates, so 2 should be the same
    expect_true( sum( mtab$dels == 0 ) == 2 )


    # Now without blocking
    mtab <- compare_methods( Yobs ~ T.x | S.id, data=fakeCRT,
                             include_method_characteristics = TRUE,
                             include_disfavored = TRUE )
    mtab
    expect_true( is.data.frame(mtab) )
    expect_true( sum( is.na(mtab) ) == 0 )
    mtab

    mtab_cov <- compare_methods( Yobs ~ T.x | S.id, data=fakeCRT,
                                 include_method_characteristics = FALSE,
                                 include_disfavored = TRUE,
                                 control_formula = ~ X.jk + C.ijk )

    mtab_cov
    mtab$dels = mtab$ATE_hat - mtab_cov$ATE_hat
    expect_true( sum( mtab$dels == 0 ) == 2 )

    mm = left_join( mtab, mc, by = "method" )
    if ( FALSE) {
        mm %>%
            dplyr::select( method, weight.x, weight.y ) %>%
            filter( weight.x != weight.y )
        names(mm)
    }
    expect_true( all( mm$weight.x == mm$weight.y ) )



})




test_that("categorical covariates handled", {

    data( fakeCRT )
    cc <- compare_methods( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                           patch_data = FALSE, include_MLM = FALSE,
                           control_formula = ~ X )
    expect_true( is.data.frame(cc) )

    cc <- compare_methods( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                           patch_data = TRUE, include_MLM = FALSE,
                           control_formula = ~ X )
    expect_true( is.data.frame(cc) )

})



test_that("missing data handled as desired", {
    data( fakeCRT )
    fake2 = fakeCRT
    head( fake2 )
    fake2$Yobs[1:10] = NA
    fake2$S.id[5:20] = NA
    fake2$X.jk[(1:30)*2] = NA
    fake2$T.x[c(1, 5, 17, 61)] = NA
    fake2$D.id[c( 5, 66 )] = NA

    nrow(fake2)

    expect_warning( mtab <- compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                             patch_data = FALSE,
                                             include_method_characteristics = FALSE) )

    expect_true( is.data.frame(mtab) )

    expect_warning( mtab_cov <-  compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                                  control_formula = ~ X.jk + C.ijk,
                                                  patch_data = FALSE,
                                                  include_method_characteristics = FALSE ) )
    expect_true( is.data.frame(mtab_cov) )

    mtab_cov

    mtab_cov2 <-  compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                   control_formula = ~ X.jk + C.ijk,
                                   include_method_characteristics = FALSE,
                                   warn_missing = FALSE )

    expect_true( is.data.frame(mtab_cov2) )

    # Kept more data, lower standard errors!
    expect_true( mean( mtab_cov2$SE_hat ) < mean( mtab_cov$SE_hat ) )

    # Check that if missingness drops all tx in a district, we stop.
    fake2$T.x[ is.na( fake2$D.id ) | is.na( fake2$T.x ) | (fake2$D.id == 1 & fake2$T.x == 0) ] = NA

    expect_warning( tt <- make_block_table( Yobs ~ T.x | S.id | D.id, data=fake2,
                                            control_formula = ~ X.jk + C.ijk ) )
    expect_true( tt$p.tx[[1]] == 1 )

    # Will throw warning about dropping all-tx blocks when trying to analyze.
    expect_warning( expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                                     control_formula = ~ X.jk + C.ijk,
                                                     patch_data = FALSE,
                                                     include_method_characteristics = FALSE ) ) )

    # Patching data will not help with all tx or all co blocks.
    expect_warning( expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                                     control_formula = ~ X.jk + C.ijk,
                                                     patch_data = TRUE, warn_missing = TRUE,
                                                     include_method_characteristics = FALSE ) ) )
    #expect_true( is.data.frame(mtab_cov) )
})





test_that("covariate adjustment works", {

    data( fakeCRT )

    c1 <- compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT )


    set.seed( 40440 )
    fakeCRT$X = fakeCRT$Yobs + rnorm( nrow(fakeCRT), sd=0.25 )
    cor( fakeCRT$X, fakeCRT$Yobs )

    c2 <- compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                           control_formula = ~ X )


    expect_true( min( c1$ATE_hat - c2$ATE_hat  ) > 0.1 )
    expect_true( all( c2$SE_hat / c1$SE_hat < 0.2 ) )



    # Check for matched pairs experiment
    mp <- fakeCRT %>%
        group_by( T.x, D.id ) %>%
        mutate( blockID = cur_group_id() ) %>%
        ungroup()

    # individual level covariate
    mp$X = mp$Yobs + rnorm( nrow(mp), sd=0.75 )
    mp
    str( mp$X )

    mp <- mp %>%
        group_by( blockID ) %>%
        mutate( W = mean(X) )

    describe_clusterRCT( Yobs ~ T.x | blockID | D.id, data = mp,
                         control_formula = ~ X )
    describe_clusterRCT( Yobs ~ T.x | blockID | D.id, data = mp,
                         control_formula = ~ W )

    expect_message( expect_message (
        c1 <- compare_methods( Yobs ~ T.x | blockID | D.id, data=mp,
                           include_method_characteristics = FALSE )
    ))
    expect_message( expect_message (
        c2 <- compare_methods( Yobs ~ T.x | blockID | D.id, data=mp,
                           control_formula = ~ X,
                           include_method_characteristics = FALSE )
    ))
    expect_message( expect_message (
        c3 <- compare_methods( Yobs ~ T.x | blockID | D.id, data=mp,
                           control_formula = ~ W,
                           include_method_characteristics = FALSE )

    ))
    c1
    c2
    c3

    cc = left_join( c1[1:3], c2[1:3], by = "method", suffix = c("", ".X") ) %>%
        left_join( c3[1:3], by = "method", suffix = c( "", ".W") ) %>%
        relocate( method, ATE_hat, ATE_hat.X, ATE_hat.W,
                  SE_hat, SE_hat.X, SE_hat.W )
    cc

    ccc = filter( cc, !is.na( SE_hat.W ) )
    if ( FALSE ) {
        sort( ccc$SE_hat.X / ccc$SE_hat )
        sort( ccc$SE_hat.W / ccc$SE_hat )
    }
    expect_true( min( c1$ATE_hat - c2$ATE_hat, na.rm=TRUE ) > 0.1 )
    expect_true( min( c1$ATE_hat - c3$ATE_hat, na.rm=TRUE ) > 0.1 )
    #expect_true( all( ccc$SE_hat.W / ccc$SEc2$SE_hat / c1$SE_hat < 0.2 ) )


})








