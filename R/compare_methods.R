


method_characteristics <- function() {

    mc <- tribble( ~ fullname,              ~ method, ~ weight, ~ population, ~ biased,
                   "DB (cluster-finite)",   ~ "DB",   "cluster",  "finite",    "0" )

    return( mc )
}



#' Compare different estimators for a (possibly blocked) cluster
#' randomized trial.
#'
#' This function calculates the point estimates and SE estimates for a
#' variety of the estimators implemented by this package.
#'
#' @param Yobs vector observed outcomes (or column name in data)
#' @param Z vector of assignment indicators (1==treated) (or column
#'   name in data)
#' @param B vector of block ids (or column name in data)
#' @param siteID site ids (variable name as string if data frame
#'   passed) (if randomization blocks are nested in site).
#' @param data frame holding Y, Z, B and (possibly a column with name
#'   specified by siteID).
#' @param include_MLM Include MLM estimators
#' @param include_DB Include Design-Based estimators (taken from
#'   RCTYes documentation and prior literature).
#' @param include_LM Include Linear Model-based estimators (including
#'   Huber-White SEs, etc.)
#' @param include_agg Include estimators applied to aggregated data,
#'   aggregated at the cluster level.
#' @param control_formula What variables to control for, in the form
#'   of "~ X1 + X2".
#' @param include_method_characteristics Include details of the
#'   methods (target estimands and sampling framework assumed) in the
#'   return value.
#'
#' @return Dataframe of point estimates and standard errors for each
#'   method considered. If \code{include_method_characteristics=TRUE}
#'   also add some features of the methods as additional columns.
#' @export
compare_methods <- function(formula,
                            data = NULL,
                            control_formula = NULL,
                            weight = c( "individual", "cluster" ),
                            include_MLM = TRUE,
                            include_DB = TRUE,
                            include_LM = TRUE,
                            include_agg = TRUE,
                            include_method_characteristics = TRUE ) {

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula=formula, data=data,
                                    control_formula = control_formula )
    }

    n <- nrow( data )
    #Quick check that input is correct
    if (is.numeric(data$Z) == FALSE ) {
        stop("Treatment indicator should be vector of ones and zeros")
    }
    if ((sum(data$Z == 1) + sum(data$Z == 0)) != n) {
        stop("Treatment indicator should be vector of ones and zeros")
    }

    # aggregate data once
    aggdat = aggregate_data(data, control_formula)
    control_formula_agg = attr( aggdat, "control_formula" )


    summary_table <- data.frame()

    if (include_DB) {
        db_res_i <- design_based_estimators(formula = NULL, data = aggdat,
                                            control_formula = control_formula_agg,
                                            weight = "individual",
                                            aggregated = TRUE)
        db_res_c <- design_based_estimators(formula = NULL, data = aggdat,
                                            control_formula = control_formula_agg,
                                            weight = "cluster",
                                            aggregated = TRUE)

        summary_table = dplyr::bind_rows( summary_table, db_res_i, db_res_c )
    }


    if (include_LM) {
        lms <- linear_model_estimators(formula = NULL, data = data,
                                       control_formula = control_formula )
        summary_table = dplyr::bind_rows(summary_table, lms)
    }


    if (include_agg) {
        lms <- aggregation_estimators( formula = NULL, data = aggdat,
                                       control_formula = control_formula_agg,
                                       aggregated = TRUE )
        summary_table = dplyr::bind_rows(summary_table, lms)
    }

    if (include_MLM) {
        mlms <- MLM_estimators(formula = NULL, data = data,
                               control_formula = control_formula )
        summary_table <- dplyr::bind_rows( summary_table, mlms )
    }

    # Add info on the methods (e.g., what estimand they are targeting)
    if (include_method_characteristics) {
        mc <- method_characteristics()

        if (!is.null(control_formula)) {
            mc$method <- paste0(mc$method, "-adj")
        }
        summary_table <- merge( summary_table, mc, by = "method",
                                all.x = TRUE, all.y = FALSE)
    }

    if ( nrow( summary_table ) > 0 ) {
        summary_table = tibble::remove_rownames( summary_table ) %>%
            arrange( method )
    }

    return(summary_table)
}
