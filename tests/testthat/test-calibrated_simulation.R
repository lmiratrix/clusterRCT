

test_that("calibrated simulation works and summary methods work", {

    data( fakeCRT )

    rs <- run_rerandomize_simulation( Yobs ~ T.x | S.id | D.id,
                                      data = fakeCRT,
                                      R = 3,
                                      summarize_results = FALSE )

    expect_true( is.data.frame(rs) )

    rss = summarize_simulation_results( rs, summarize_results = "method" )
    rss
    expect_true( is.data.frame(rss) )
    expect_true( "method" %in% names(rss) )

    rs2 = summarize_simulation_results( rs, summarize_results = "cross" )
    rs2
    expect_true( "runID" %in% names(rs2) )

    rs3 = summarize_simulation_results( rs, summarize_results = "cross-agg" )
    rs3
    expect_equal( nrow(rs3), 3 )
    expect_equal( rs3$group, c( "Cluster", "Person", "Overall" ) )
})
