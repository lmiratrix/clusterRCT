

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


test_that("canonical (non-MLM) weight point estimates match their corresponding linear-model estimators", {

    data( fakeCRT )

    formula = Yobs ~ T.x | S.id | D.id

    cw <- canonical_weight_estimators( formula, data = fakeCRT )
    cw_val <- function( m ) unname( cw$ATE_hat[ cw$method == m ] )

    cm <- compare_methods( formula, data = fakeCRT,
                           include_gee = FALSE, include_dumb = TRUE,
                           include_disfavored = TRUE )
    cm_val <- function( m ) unname( cm$ATE_hat[ cm$method == m ] )

    # The non-MLM canonical weight combinations are, per the docstring,
    # meant to exactly reproduce the corresponding linear-model
    # estimator's point estimate -- this is the actual "canonical
    # weight formula matches the regression" validation these
    # estimators exist for.
    expect_equal( cw_val( "WT_Person-Person" ), cm_val( "LR-FIpw-crve" ), tolerance = 1e-6 )
    expect_equal( cw_val( "WT_Person-FE" ),     cm_val( "LR-FE-crve" ),   tolerance = 1e-6 )
    expect_equal( cw_val( "WT_Person-Block" ),  cm_val( "LR-FIbw-crve" ), tolerance = 1e-6 )
    expect_equal( cw_val( "WT_Cluster-Cluster" ), cm_val( "LRcw-FIcw-crve" ), tolerance = 1e-6 )
    expect_equal( cw_val( "WT_Cluster-FE" ),      cm_val( "LRcw-FE-crve" ),   tolerance = 1e-6 )
    expect_equal( cw_val( "WT_Cluster-Block" ),   cm_val( "LRcw-FIbw-crve" ), tolerance = 1e-6 )
})
