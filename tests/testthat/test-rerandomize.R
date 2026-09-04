

test_that("rerandomize preserves cluster structure and block treatment proportions (blocked)", {

    set.seed( 40404 )

    clusterID <- rep( 1:12, times = c( 4,5,6,3,4,5,6,3,4,5,6,3 ) )
    blockID <- rep( c( "A", "B", "C" ), each = 4 )[ as.numeric( factor( clusterID ) ) ]
    # Build clusterID/blockID/Z at the cluster level, then expand to the person level.
    cluster_tbl <- tibble( clusterID = 1:12,
                           blockID = rep( c("A","B","C"), each = 4 ),
                           Z = rep( c(1,1,0,0), 3 ) )
    n_per_cluster <- c( 4,5,6,3,4,5,6,3,4,5,6,3 )
    person_tbl <- cluster_tbl[ rep( 1:12, times = n_per_cluster ), ]

    newZ <- rerandomize( person_tbl$Z, person_tbl$clusterID, person_tbl$blockID )

    expect_equal( length( newZ ), nrow( person_tbl ) )

    # Treatment status must be constant within each cluster.
    per_cluster_levels <- tapply( newZ, person_tbl$clusterID, function(z) length(unique(z)) )
    expect_true( all( per_cluster_levels == 1 ) )

    # Per-block proportion of *clusters* treated should match the original design (2 of 4 per block).
    new_cluster_Z <- tapply( newZ, person_tbl$clusterID, unique )
    new_block_id <- tapply( person_tbl$blockID, person_tbl$clusterID, unique )
    new_prop_by_block <- tapply( new_cluster_Z, new_block_id, mean )
    expect_equal( as.vector( new_prop_by_block ), rep( 0.5, 3 ) )
})


test_that("rerandomize preserves cluster structure and treatment proportion (unblocked)", {

    set.seed( 40405 )

    clusterID <- rep( 1:10, times = 3:12 )
    Z <- rep( c(1,1,1,1,0,0,0,0,0,0), times = 3:12 )  # 4 of 10 clusters treated

    newZ <- rerandomize( Z, clusterID, blockID = NULL )

    expect_equal( length( newZ ), length( Z ) )

    per_cluster_levels <- tapply( newZ, clusterID, function(z) length(unique(z)) )
    expect_true( all( per_cluster_levels == 1 ) )

    new_cluster_Z <- tapply( newZ, clusterID, unique )
    expect_equal( mean( new_cluster_Z ), 0.4 )
})


test_that("rerandomize's implied treatment probability matches the target proportion on average", {

    set.seed( 40406 )

    clusterID <- rep( 1:8, each = 5 )
    blockID <- rep( c("A","B"), each = 20 )
    # 3 of 4 clusters treated in block A, 1 of 4 in block B.
    Z <- rep( rep( c(1,1,1,0), 5 )[1:20], 1 )
    Z <- c( rep( c(1,1,1,0), 5 ), rep( c(1,0,0,0), 5 ) )

    R <- 200
    sims <- replicate( R, rerandomize( Z, clusterID, blockID ), simplify = FALSE )

    cluster_treated_rate <- sapply( sims, function(z) tapply( z, clusterID, unique ) )
    # rows are clusters 1-8; average treatment rate per cluster should track block-level proportion.
    avg_rate_by_cluster <- rowMeans( cluster_treated_rate )
    block_of_cluster <- tapply( blockID, clusterID, unique )

    expect_equal( mean( avg_rate_by_cluster[ block_of_cluster == "A" ] ), 0.75, tolerance = 0.1 )
    expect_equal( mean( avg_rate_by_cluster[ block_of_cluster == "B" ] ), 0.25, tolerance = 0.1 )
})
