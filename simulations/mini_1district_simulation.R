


# Simple simulation as sanity check of SEs etc. This simulation
# generates code for a single district CRT (i.e., we just have
# clusters).
#
# We can thus use the blkvar package and then re-assign tx with
# randomizr, so we can have clusters of unequal size, etc.




set.seed( 1039 )
library( tidyverse )
library( PUMP )
library( blkvar )
library( clusterRCT )



one_run <- function() {

    sim.data <- generate_multilevel_data( n.bar = 50, J = 30, p = 0.3,
                                          size.impact.correlate = 1,
                                          variable.n = TRUE )
    sim.data <- mutate( sim.data,
                        Z = randomizr::cluster_ra(clusters = sid),
                        Yobs = ifelse( Z == 1, Y1, Y0 ) )

    c1 <- compare_methods( Yobs ~ Z | sid, data=sim.data,
                           include_method_characteristics = FALSE )

    c1$tau = mean( sim.data$Y1 - sim.data$Y0 )

    c1
}

one_run()


rps = map_df( 1:20, ~ one_run(), .id = "runID" )

head( rps )


rps %>% group_by( method ) %>%
    summarise( tau = mean( tau ),
               EATE = mean( ATE_hat ),
               SE = sd( ATE_hat ),
               ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
    mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )





# Profile the above code
