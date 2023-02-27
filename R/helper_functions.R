##
## Utility functions for the rest of the package
##




library( Formula )
library( tidyverse )



scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}



#### Expand control variables to dummies as needed ####



#' Take data with possibly categorical control variables and expand to
#' dummy variables.
#'
#' @return List of three objects: the revised data frame, revised
#'   control_formula, and final list of names of control variables.
expand_control_variables <- function( data, control_formula ) {
    # Get control variables and expand control matrix with dummy variables, etc.
    controls = model.matrix( control_formula, data )
    stopifnot( colnames(controls)[1] == "(Intercept)" )
    controls = controls[,-1, drop=FALSE]

    c.names = make.names( colnames(controls), unique=TRUE )
    colnames(controls) = c.names

    # Update dataframe to expanded controls
    data = data %>%
        dplyr::select( -all_of( formula.tools::rhs.vars(control_formula) ) ) %>%
        cbind( as.data.frame( controls ) )

    control_formula = as.formula( paste0( "~ ", paste( c.names, collapse = " + " ) ) )

    return( list( data = data,
            control_formula = control_formula,
            c.names = c.names ) )
}





#### Aggregate data to clusters ####


#' Aggregate data to the cluster level.
#'
#' Also implements checks to ensure that all clusters are entirely
#' treated or not treated.
#'
#' If covariates specified, categorical covariates will be converted
#' to dummy variables and then averaged (even if they are only level-2
#' or level-3 covariates)
#'
#' @param data  Dataframe with 'clusterID' and 'Z' as columns.
#'   'siteID' optional column.
#'
#' @return tibble of cluster-aggregated data, including Ybar, n,
#'   siteID, clusterID, and Z
aggregate_data <- function( data, control_formula = NULL ) {
    # Aggregation to the cluster level

    datagg = NA
    stopifnot( "clusterID" %in% names(data) )
    stopifnot( "Z" %in% names(data) )


    if ( !is.null( control_formula ) ) {

        dd = expand_control_variables( data, control_formula )
        data = dd$data

        my.vars = dd$c.names
        datagg <-
            datagg <- data %>%
            group_by( across( any_of( c( "siteID", "clusterID", "Z" ) ) ) ) %>%
            summarise( Ybar = mean( Y ),
                       n = n(),
                       across( all_of( my.vars ), mean ),
                       .groups = "drop" )
    } else {
        datagg <-
            datagg <- data %>%
            group_by( across( any_of( c( "siteID", "clusterID", "Z" ) ) ) ) %>%
            summarise( Ybar = mean( Y ),
                       n = n(),
                       .groups = "drop" )
    }

    J = nrow( datagg )
    stopifnot( J == length(unique(data$clusterID) ) )

    return( datagg )
}





#### Functions to help process control functions  ####


#' Make a canonical regression formula (no FE), possibly with control
#' variables.
#'
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z Name of treatment variable (assumed to exist in data)
#' @param clusterID Name of cluster ID variable (assumed to exist in
#'   data)
#' @param siteID Name of site ID variable (assumed to exist in data,
#'   if this is not null).
#' @param control_formula What variables to control for, in the form
#'   of "~ X1 + X2".
#' @param data Dataframe holding all variables to be used in formula.
#' @param interacted  TRUE means include treatment by site
#'   interactions.  Will override FE flag and set FE to true.
#' @param FE TRUE means include site dummy variables.
#'
#' @return Something like "Yobs ~ 1 + Z" or "Yobs ~ 1 + Z + X1 + X2"
#'
#' @export
make_regression_formula = function( Yobs = "Yobs", Z = "Z",
                                    clusterID = "clusterID",
                                    siteID = "siteID",
                                    control_formula = NULL,
                                    interacted = FALSE,
                                    FE = FALSE,
                              data = NULL) {

    if ( interacted ) {
        FE = TRUE
    }

    # Formula types
    # Yobs ~ 1 + Z + X
    # Yobs ~ 0 + siteID + Z + X
    # Yobs ~ 0 + siteID + Z:siteID + X

    stem = "1"
    if ( FE ) {
        stem = "0"
    }

    txS = "Z"
    if ( interacted ) {
        txS = paste0( Z, ":", siteID )
    }

    if ( FE || interacted ) {
        txS = paste0( txS, " + ", siteID )
    }

    if( !is.null( control_formula ) ) {
        if (length(formula.tools::lhs.vars(control_formula)) != 0 |
            length(formula.tools::rhs.vars(control_formula)) < 1) {
            stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
        }

        if (!is.null(data)) {
            control.vars <- formula.tools::get.vars(control_formula, data = data)
            if (any(!(control.vars %in% colnames(data)))) {
                stop("Some variables in control_formula are not present in your data.")
            }
        }

        c.names <- formula.tools::rhs.vars(control_formula)
        c.names <- paste( c.names, collapse =" + " )

        new.form <- sprintf( "%s ~ %s + %s + %s", Yobs, stem, txS, c.names)
        return( as.formula( new.form ) )
    } else {
        new.form <- sprintf( "%s ~ %s + %s", Yobs, stem, txS )
        return( as.formula( new.form ) )
    }
}








#' Make canoncial data for analysis
#'
#' Given a formula, assemble dataframe of outcome, treatment, and any
#' clustering variables all with canonical names (Y, Z, clusterID,
#' siteID)
#'
#' @param formula Notation for Y ~ Z | clusterID | siteID (| siteID is
#'   optional).
#' @param data Dataframe to pull data from.
#' @return canonical dataframe.
#'
#' @return Dataset with now-canonical variable names, and no
#'   extraneous variables.
#'
#' @noRd
make_canonical_data <- function(formula, control_formula = NULL, data,
                                give_default_site = FALSE ) {

    formula <- as.formula(formula)
    rhs <- formula.tools::rhs.vars( formula )
    parts <- strsplit( rhs, split=" \\| " )[[1]]
    n_part <- length( parts )
    outcome <- formula.tools::lhs.vars(formula)
    if ( (  n_part != 2 && n_part != 3 ) ||
         (length(outcome) != 1 || length(rhs) != 1) ) {
        stop("The formula argument must be of the form outcome ~ treatment | clusterID or outcome ~ treatment | clusterID | siteID.")
        #stop( "Need cluster ID (and possibly site ID) in formula, e.g. , e.g., Y ~ Z | clusterID or Y ~ Z | clusterID | siteID" )
    }

    clusterID = parts[[2]]
    siteID = NULL
    if ( n_part == 3 ) {
        siteID = parts[[3]]
    }
    main.vars = c( outcome, parts )
    if (any(!(main.vars %in% colnames(data)))) {
        stop("Outcome, treatment variables, or cluster variables in formula are not present in your data.")
    }

    Y = data[[ outcome ]] # model.part( formula, data=data, lhs=1, drop = TRUE)
    Z = data[[ parts[[1]] ]] # model.part( formula, data=data, rhs=1, drop = TRUE)
    clusterID = data[[ clusterID ]] # model.part( formula, data=data, rhs=2, drop=TRUE )
    clusterID = as.factor( clusterID )

    new_dat = data.frame( Y = Y, Z = Z, clusterID = clusterID )
    if ( n_part == 3 ) {
        siteID = data[[ siteID ]]
        new_dat$clusterID = interaction(siteID, clusterID, drop=TRUE)
        new_dat$siteID = as.factor( siteID )
    }

    if ( !is.null( control_formula ) ) {
        if(length(formula.tools::lhs.vars(control_formula)) != 0 |
           length(formula.tools::rhs.vars(control_formula)) < 1) {
            stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
        }
        control.vars <- formula.tools::get.vars(control_formula, data = data)
        if (any(!(control.vars %in% colnames(data)))) {
            stop("Some variables in control_formula are not present in your data.")
        }

        if ( any( control.vars %in% colnames(new_dat) ) ) {
            stop( "Control variables with canonical names of Y, Z, clusterID, or siteID are not allowed." )
        }
        xd = data %>% dplyr::select( dplyr::all_of( control.vars ) )
        new_dat = bind_cols( new_dat, xd )
    }


    # Check for valid treatment variable
    n_levels = length( unique( new_dat$Z ) )
    if ( n_levels == 1 ) {
        warning( glue::glue( "Only single level in treatment variable {parts[[1]]}" ) )
    }
    if ( n_levels > 2 ) {
        stop( sprintf( "Identified treatment variable '%s' has more than two values. Did you swap treatment and block?",
                       parts[[1]] ) )
    }

    if ( give_default_site && is.null( new_dat$siteID ) ) {
        new_dat$siteID = "Sall"
    }

    return( new_dat )
}





# # #
# if ( FALSE ) {

# make_base_formula()
# make_base_formula( control_formula = ~ X1 + X2 + X3 )
# make_FE_formula(  )
# make_FE_formula( control_formula = ~ X1 + X2 + X3 )

# make_FE_int_formula( )
# make_FE_int_formula( control_formula = ~ X1 + X2 + X3 )

# }



