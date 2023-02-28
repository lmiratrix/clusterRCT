


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
