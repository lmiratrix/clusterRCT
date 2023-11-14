

# Patch a dataset with mean imputation of the covariates


#' Fill in missing values and convert to canonical form
#'
#' This method will fill in means for all missing covariate values and
#' return a dataset where rows with missing outcome, cluster or block
#' membership, and treatment assignment are all dropped and everything
#' else is filled in.  To do this, it will convert categorical
#' covariates to a set of dummy variables (dropping the reference
#' group).
#'
#' It will drop missing data indicators that are colinear with prior
#' missing data indicators.
#'
#' If there are categorical covariates, then it will convert those to
#' dummy variables regardless of missingness.  The `control_formula`
#' attribute handed back will give the updated control_formula in this
#' case.
#'
#' This attribute will also have all missing value indicators added in
#' as covariates.
#'
#' @export
patch_data_set <- function( formula = NULL,
                            data = NULL,
                            control_formula = NULL) {

    if ( !is.null(formula) ) {
        data = make_canonical_data( formula=formula, data=data,
                                control_formula = control_formula,
                                drop_missing = FALSE, warn_missing = FALSE )
    }

    mod_data = NULL
    if ( !is.null( control_formula ) ) {
        mod_data <- clusterRCT:::expand_control_variables( data, control_formula )
        data = mod_data$data
        control_formula = mod_data$control_formula
    }

    data = filter( data,
                   !is.na( Z ), !is.na( Yobs ),
                   !is.na( clusterID ) )
    if ( "blockID" %in% colnames(data) ) {
        data <- filter( data, !is.na( blockID ) )
    }

    if ( is.null( control_formula ) ) {
        return( data )
    }

    # Locate missing data and do imputation on the columns with
    # missingness.  Also generate a dummy var.  Make sure they are not
    # colinear.
    miss = colSums(is.na(data))

    which_miss = which( miss > 0 )
    if ( length( which_miss ) > 0 ) {
        missInd = 0 + is.na( data[ , which_miss, drop=FALSE ] )
        colnames(missInd) = paste0( colnames(missInd) , "_mi" )
        qr = qr( missInd )
        missInd <- missInd[ , qr$pivot[ 1:qr$rank ], drop=FALSE ] %>%
            as.data.frame()

        data <- data %>%
            mutate( across(all_of(which_miss), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

        data <- bind_cols(data, missInd)

        mod_data$control_formula <- update( mod_data$control_formula,
                                            as.formula( paste0( ". ~ . + ", paste( colnames(missInd), collapse = " + " ) ) ) )
    }

    if ( !is.null( mod_data ) ) {
        attr( data, "control_formula" ) <- mod_data$control_formula
    } else {
        attr( data, "control_formula" ) <- NULL
    }

    data
}



#### For debugging and demonstration ####

if ( FALSE ) {

    data <- tribble( ~ X, ~ R, ~ Q, ~ W,
                   1,  2, 3, 4,
                   NA, 2, 4, 5,
                   NA, NA,4, 6,
                   2,  3, NA, 7)

    data$T.x = c( 1, 1, 0, 0 )
    data$S.id = c( 10, 5, 10, 5 )
    data = bind_rows( data, data, data, data, data)
    data$D.id = rep( c("A","B"), c(8, 12) )
    data$Y = rnorm( nrow(data), sd=10 )
    data$X[1:5] = 1:5
    data$W[5:10] = 4
    data$Y[12] = NA
    data$T.x[9] = NA
    data$D.id[11] = NA
    data$S.id[2] = NA
    data$RR = sample( c("A","B","C"), nrow(data), replace=TRUE )
    data$RR[17] = NA
    data

    ptch <- patch_data_set( Y ~ T.x | S.id | D.id, data=data,
                    control_formula = ~ X + R + Q + W + RR )
    ptch




    #### Code for debugging ####
    data <- clusterRCT:::make_canonical_data( Y ~ T.x | S.id | D.id, data=data,
                                              control_formula = ~ X + R + Q + W + RR,
                                              drop_missing = FALSE )
    data

    formula = Y ~ T.x | S.id | D.id
    control_formula = ~ X + R + Q + W + RR




    library(tidyverse)

    # Replace missing values with the mean of each column
    data <- data %>%
        mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))
    data

    # Now, 'data' contains mean-imputed values for missing data


    #data[is.na(data)] <- sapply(data, function(x) mean(x, na.rm = TRUE))
    data

    data( "fakeCRT" )
    data <- clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
    head( data )

}
