library( testthat )


test_that( "aggregation works", {
    dd = data.frame( Y = 1:20,
                     clusterID = sample(LETTERS[1:5], 20, replace=TRUE),
                     siteID = sample( LETTERS[1:2], 20, replace=TRUE ),
                     x = 1:20,
                     x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    dd$Z = randomizr::cluster_ra( clusters = dd$clusterID )
    dd = clusterRCT:::make_canonical_data( Y ~ Z | clusterID | siteID, data=dd,
                                           control_formula = ~ x + x2 )
    dd$Z = randomizr::cluster_ra( clusters = dd$clusterID )
    aa <- clusterRCT:::aggregate_data( dd )
    expect_equal( names(aa), c("siteID", "clusterID","Z", "Ybar", "n"))

    bb <-     clusterRCT:::aggregate_data( dd, control_formula = ~ x + x2 )
    bb
    expect_true( all( c( "x", "x2B", "x2C" ) %in% names(bb) ) )
    expect_true( nrow(bb) == nrow( aa ) )
    expect_true( all( bb$Ybar == aa$Ybar ) )


    dd$siteID = NULL
    a2 = clusterRCT:::aggregate_data( dd, control_formula = ~ x + x2 )
    a2

    a4 = clusterRCT:::aggregate_data( dd, control_formula = ~ x )
    a4
    expect_true( all( a2$x == a4$x ) )
})




test_that("formula generator works", {
    dd = data.frame( Y = 1:10, tx = rep(2,10),
                     c = sample(c("A","B"), 10, replace=TRUE),
                       SITE = sample( LETTERS[1:4], 10, replace=TRUE ),
                       x = 1:10,
                       x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    dd

    rr = make_regression_formula()

    expect_equal( as.character(rr), "Yobs ~ 1 + Z" )

    rr = make_regression_formula( siteID = "SITE", FE = TRUE )
    expect_equal( as.character(rr), "Yobs ~ 0 + Z + SITE" )

    rr = make_regression_formula( siteID = "SITE", interacted = TRUE )
    expect_equal( as.character(rr), "Yobs ~ 0 + Z:SITE + SITE" )

    rr = make_regression_formula( siteID = "SITE", FE = TRUE )
    expect_equal( as.character(rr), "Yobs ~ 0 + Z + SITE" )

    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  siteID = "SITE", interacted = TRUE,
                                  control_formula = ~ x + x2 )
    expect_equal( as.character(rr), "Y ~ 0 + Z:SITE + SITE + x + x2" )

    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  siteID = "SITE", FE = TRUE,
                                  control_formula = ~ x + x2 )
    expect_equal( as.character(rr), "Y ~ 0 + Z + SITE + x + x2" )

    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  siteID = "SITE", FE = FALSE,
                                  control_formula = ~ x + x2 )
    expect_equal( as.character(rr), "Y ~ 1 + Z + x + x2" )


})






test_that("formula generator random effects", {

    dd = data.frame( Y = 1:10, tx = rep(2,10),
                     c = sample(c("A","B"), 10, replace=TRUE),
                     SITE = sample( LETTERS[1:4], 10, replace=TRUE ),
                     x = 1:10,
                     x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    dd

    rr = make_regression_formula( cluster_RE = TRUE )
    a <- as.character(rr)
    a
    expect_true( a == "Yobs ~ 1 + Z + (1 | clusterID)" )


    rr = make_regression_formula( siteID = "SITE", FE = TRUE, cluster_RE = TRUE )
    a <- as.character(rr)
    a
    expect_equal( a, "Yobs ~ 0 + Z + SITE + (1 | clusterID)" )


})







test_that("making canonical data works", {

    odat = data.frame( Y = 1:10, tx = rep(2,10), c = sample(c("A","B"), 10, replace=TRUE),
                      s = sample( c("X","Y"), 10, replace=TRUE ),
                      x = 1:10,
                      x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    odat
    odat$tx[1:5] = 1

    # dat = make_canonical_data(formula=Y ~ tx | c | s, data=odat)

    dat = clusterRCT:::make_canonical_data( formula = Y ~ tx | c | s, data=odat)

    dat
    expect_true( ncol(dat) == 4 )
    expect_equal( names(dat), c( "Yobs", "Z", "clusterID", "siteID" ) )


    dat2 = clusterRCT:::make_canonical_data(formula=Y ~ tx | c | s, data=odat, control_formula = ~ x + x2 )
    head( dat2 )
    expect_true( ncol(dat2) == 6 )

    expect_equal( names(dat2), c( "Yobs", "Z", "clusterID", "siteID", "x", "x2" ) )

})
