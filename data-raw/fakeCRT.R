## code to prepare `fakeCRT` dataset goes here


set.seed( 1039 )
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
    , ICC.3 = 0.2             # district intraclass correlation
    , omega.3 = 0.1           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.1           # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)


``
sim.data <- gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
fakeCRT = slice_sample(sim.data, n = nrow(sim.data)*0.5 )


library( lme4 )
M <- lmer( Yobs ~ 1 + T.x + (1|S.id) + (1|D.id), data=fakeCRT )
arm::display(M)


v <- (1.06^2 + 0.56^2+ 0.66^2)
v
0.19 / v

0.66^2 / (1.06^2 + 0.56^2+ 0.66^2)
0.56^2 / (1.06^2 + 0.56^2+ 0.66^2)


usethis::use_data(fakeCRT, overwrite = TRUE)
