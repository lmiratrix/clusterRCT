

test_that("check data integrity works", {

    data( fakeCRT )
    head( fakeCRT )

    expect_true( is_nested( fakeCRT$S.id, fakeCRT$D.id ) )
    expect_true( !is_nested( fakeCRT$D.id, fakeCRT$S.id ) )

    dat1 = clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
    rs <- clusterRCT:::check_data_integrity( dat1 )
    expect_true( rs )

    dat1$Z[ 1:3 ] = 1
    expect_error( clusterRCT:::check_data_integrity( dat1 ) )

    dat1$Z[1:3] = 2
    expect_error( clusterRCT:::check_data_integrity( dat1 ) )

    d2 = dat1
    d2$Z[ d2$SiteID == d2$SiteID[[1]] ] = 0
    expect_error( clusterRCT:::check_data_integrity( d2 ) )

    dat1 = dat1[ -c(1:3), ]

    dat1$clusterID[1] = dat1$clusterID[2]
    expect_error( clusterRCT:::check_data_integrity( dat1 ) )
    dat1 = dat1[ -c(1:3), ]

    dat2 = clusterRCT:::make_canonical_data( Yobs ~ T.x | D.id | S.id, data=fakeCRT )
    expect_error( clusterRCT:::check_data_integrity( dat2 ) )
})




test_that("check missing data doesn't crash", {

    data( fakeCRT )
    fakeCRT$Yobs[1:10] = NA
    expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT ) )


    fakeCRT$T.x[5:20] = NA
    expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT ) )

    data(fakeCRT)
    fakeCRT$X.jk[30:40] = NA
    expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                     control_formula = ~ X.jk ) )

    data(fakeCRT)
    fakeCRT$S.id[ 100:150 ] = NA
    expect_warning( describe_clusterRCT( Yobs ~ T.x | S.id | D.id,
                                         data=fakeCRT, control_formula = ~ X.jk ) )

    fakeCRT$D.id[ 70:150 ] = NA
    expect_warning( clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id,
                                                      data=fakeCRT, control_formula = ~ X.jk ) )

})
