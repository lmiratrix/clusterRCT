

test_that("canonical weight estimators run", {


    data( fakeCRT )

    formula =  Yobs ~ T.x | S.id | D.id
    control_formula = ~ X.jk + C.ijk

    res = canonical_weight_estimators(formula,
                                      data=fakeCRT)
    res

    res2 = canonical_weight_estimators( Yobs ~ T.x | S.id,
                                        data=fakeCRT)
    res2

    expect_true( is.data.frame(res2) )
})
