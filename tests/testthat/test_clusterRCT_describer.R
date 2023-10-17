

library( tidyverse )
library( clusterRCT )

set.seed( 1039 )


model.params.list <- list(
    M = 1                            # number of outcomes
    , J = 30                          # number of schools
    , K = 10                          # number of districts
    , nbar = 10                       # number of individuals per school
    , S.id = NULL                     # N-length vector of school assignments
    , D.id = NULL                     # N-length vector of district assignments
    , Xi0 = 0                         # scalar grand mean outcome under no treatment
    , MDES = 0.125            # minimum detectable effect size
    , R2.3 = 0.1              # percent of district variation
    , ICC.3 = 0.2             # district intraclass correlation
    , omega.3 = 0.2           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.0          # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)

data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )

data <- clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=data )
head( data )




test_that( "describer works", {


    data( fakeCRT )

    d <- describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                              control_formula = ~ V.k + X.jk + C.ijk )

    expect_true( is.clusterRCTstats(d) )

    dd = as.data.frame(d)
    expect_true( nrow( dd ) == 1 )

    # We can check nesting of School in District
    expect_true( is_nested( fakeCRT$S.id, fakeCRT$D.id ) )
    expect_true( !is_nested( fakeCRT$D.id, fakeCRT$S.id ) )

    # We can calculate statistics for each district
    st <- make_block_table( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                     control_formula = ~ V.k + X.jk + C.ijk )
    expect_true( nrow( st ) == 10 )


    # It will trap errors!
    fakeCRT$T.x[1:4] = 2
    expect_error( describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                         control_formula = ~ V.k + X.jk + C.ijk ) )


} )

