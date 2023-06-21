


# Simple simulation as sanity check of SEs etc

set.seed( 1039 )
library( tidyverse )
# library( PUMP )
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
    , ICC.3 = 0.2             # district intraclass correlation
    , omega.3 = 0.2           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.0          # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)



# Hack--- use binom dist to get number of big clusters (needed?)
# generate and then stack the data from little and big cluster estimates



one_run <- function() {

    sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )

    c1 <- clusterRCT::compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data, include_method_characteristics = FALSE )

    c1
}

one_run()


rps = map_df( 1:50, ~ one_run(), .id = "runID" )

head( rps )


rps %>% group_by( method ) %>%
    summarise( EATE = mean( ATE_hat ),
               SE = sd( ATE_hat ),
               ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
    mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )





# Profile the above code
