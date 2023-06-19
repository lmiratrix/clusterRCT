


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

    c1 <- clusterRCT::compare_methods( Yobs ~ Z | sid, data=sim.data,
                                       include_method_characteristics = FALSE )

    c1$tau_indiv = mean( sim.data$Y1 - sim.data$Y0 )
    c1$tau_clust = sim.data %>%
        group_by(sid) %>%
        summarize(tau = mean(Y1-Y0)) %>%
        summarize(tau = mean(tau)) %>%
        pull(tau)

    c1
}

one_run()


rps = map_df( 1:100, ~ one_run(), .id = "runID", .progress=TRUE )

head( rps )


rps %>% group_by( method ) %>%
    summarise( tau_indiv = mean( tau_indiv ),
               sdtau_clust = sd( tau_clust ),
               tau_clust = mean( tau_clust ),
               EATE = mean( ATE_hat ),
               SE = sd( ATE_hat ),
               ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
    mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )


# RESULTS:
# # A tibble: 6 x 7
#   method       tau_indiv tau_clust  EATE    SE ESE_hat calib
#   <glue>           <dbl>     <dbl> <dbl> <dbl>   <dbl> <dbl>
# 1 DB_clust         0.373     0.219 0.196 0.308   0.301 0.977
# 2 DB_indiv         0.373     0.219 0.368 0.319   0.316 0.989
# 3 LR (agg)         0.373     0.219 0.196 0.308   0.301 0.977
# 4 LR (agg, wt)     0.373     0.219 0.368 0.319   0.321 1.01
# 5 LR w/ CRSE       0.373     0.219 0.368 0.319   0.321 1.01
# 6 MLM              0.373     0.219 0.200 0.307   0.300 0.978




# Profile the above code
