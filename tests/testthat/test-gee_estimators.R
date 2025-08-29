

test_that("gee estimation works", {

        data( fakeCRT )
        sim.data = bind_rows( fakeCRT, fakeCRT[1:100, ] )

        J = length( unique( sim.data$S.id ) )

        lme <- gee_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )

        expect_true( nrow( lme) == 1 )
        expect_equal( lme$df, J - 10 - 1 )

        expect_true( all( !is.na( lme ) ) )

})
