


data( "fakeCRT" )

test_that("compare methods aggregates as expected", {

    mtab <- compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                             include_method_characteristics = FALSE)

    expect_true( is.data.frame(mtab) )

    mtab_cov <-  compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                  control_formula = ~ X.jk + C.ijk, include_method_characteristics = FALSE )

    # Ensure covariate correction changes the estimates!  (I.e.,
    # control formula is getting propagated)
    mtab$dels = mtab$ATE_hat - mtab_cov$ATE_hat
    mtab
    expect_true( all( mtab$dels != 0 ) )

    mtab <- compare_methods( Yobs ~ T.x | S.id , data=fakeCRT,
                             include_method_characteristics = FALSE)

    expect_true( is.data.frame(mtab) )
    mtab

})





test_that("missing data does not crash", {

    fake2 = fakeCRT
    head( fake2 )
    fake2$Yobs[1:10] = NA
    fake2$S.id[5:20] = NA
    fake2$X.jk[(1:30)*2] = NA
    fake2$T.x[c(1, 5, 17, 61)] = NA
    fake2$D.id[c( 5, 66 )] = NA

    nrow(fake2)

    expect_warning( mtab <- compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                             include_method_characteristics = FALSE) )

    expect_true( is.data.frame(mtab) )

    expect_warning( mtab_cov <-  compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                  control_formula = ~ X.jk + C.ijk,
                                  include_method_characteristics = FALSE ) )
    expect_true( is.data.frame(mtab_cov) )


    # Check if missingness drops all tx in a district, we stop.
    fake2$T.x[ is.na( fake2$D.id ) | is.na( fake2$T.x ) | (fake2$D.id == 1 & fake2$T.x == 0) ] = NA

    expect_warning( tt <- make_site_table( Yobs ~ T.x | S.id | D.id, data=fake2,
                     control_formula = ~ X.jk + C.ijk ) )
    expect_true( tt$p.tx[[1]] == 1 )

    # Will crash when trying to analyze.
    expect_warning( expect_error( compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                  control_formula = ~ X.jk + C.ijk,
                                  include_method_characteristics = FALSE ) ) )
    #expect_true( is.data.frame(mtab_cov) )
})
