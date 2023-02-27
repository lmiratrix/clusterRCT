

library( testthat )

library( tidyverse )
library( PUMP )

library( clusterRCT )

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
    , ICC.3 = 0.2 # district intraclass correlation
    , omega.3 = 0.1           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.1           # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)



sim.data <- gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )




test_that("linear models works", {
    head( sim.data )
    lme <- linear_model_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )
    lme
    expect_true( nrow( lme) == 1 )


    agg = aggregation_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )
    agg
    expect_true( !is.null( agg ) )

})
