

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



sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
rm( model.params.list )




test_that("DB estimators work", {
    head( sim.data )
    table( sim.data$D.id)
    aa = design_based_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )
    aa
    expect_true( nrow( aa ) == 2 )

    # No block-level.
    sim.data$D.id = NULL
    bb = design_based_estimators( Yobs ~ T.x | S.id, data=sim.data )
    bb
    expect_true( nrow(bb) == 1 )

})
