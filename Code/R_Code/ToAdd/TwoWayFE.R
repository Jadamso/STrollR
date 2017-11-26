#http://darinchristensen.github.io/2014/11/21/twoway_fe/

#------------------------------------------------------------------
##################
#' Two Way clustered SE
################## 
#' @param y 
#' @param Xs
#' @param fe1
#' @param fe2
#' @param dt
#' @param cluster
#' @return data.frame of dij and (ei-ej)^2
#  @examples
#' @export

TwoWayFE <- function(y, Xs, fe1, fe2, dt, cluster=FALSE) {
    pkgs <- c("data.table", "RcppEigen")
    sapply(pkgs, require, character.only = TRUE)

    keep_cols <- c(y, Xs, fe1, fe2)
    model_dt <- dt[, keep_cols, with = FALSE]; rm(keep_cols)
    model_dt <- na.omit(model_dt)

    num_Xs <- length(Xs)
    new_names <- c("y", paste0("x", 1:num_Xs), "fe1", "fe2")
    setnames(model_dt, 1:ncol(model_dt), new_names)

    # Sample Means:
    cols <- new_names[!grepl("fe", new_names)]
    model_dt[, paste0("mean_", cols) :=
        lapply(.SD, mean, na.rm = TRUE), .SDcols = cols]

    # Means by FE1:
    setkey(model_dt, fe1)
    model_dt[,
        paste0("mean_", cols, "_fe1") :=
            lapply(.SD, mean, na.rm = TRUE), .SDcols = cols, by = fe1]
    M <- length(unique(model_dt$fe1))

    # Means by FE2:
    setkey(model_dt, fe2)
    model_dt[,
        paste0("mean_", cols, "_fe2") :=
        lapply(.SD, mean, na.rm = TRUE), .SDcols = cols, by = fe2]
    Y <- length(unique(model_dt$fe2))

    # Demeaning:
    model_dt[, "y_tilde" := y - mean_y_fe2 - mean_y_fe1 + mean_y]

    g <- function(i) {paste0("x",i,"_tilde")}
    LHS <- sapply(1:num_Xs, g)

    f <- function(i) {
        paste0("x",i," - mean_x",i,"_fe2 - mean_x",i,"_fe1 + mean_x", i)
    }
    RHS <- paste0("list(",paste(sapply(1:num_Xs, f), collapse = ", "), ")")

    model_dt[, eval(LHS) := eval(parse(text = RHS))]

    x_cols <- grep("x\\d{1}_tilde", names(model_dt), value = TRUE)
    model_dt <- model_dt[, c("y_tilde", eval(x_cols), "fe1", "fe2"),
        with = FALSE]

    y <- model_dt$y_tilde
    X <- model_dt[, x_cols, with = FALSE]
    cluster_vec <- model_dt$fe1
    rm(model_dt)

    m <- RcppEigen::fastLm(X, y)
    names(m$coefficients) <- Xs

    ##############################################
    DoF <- m$df.residual - (M - 1) - (Y - 1) - num_Xs + 1
    # No intercept in model.

    # SEs:
    if(cluster){
        N <- length(cluster_vec)
        K <- m$rank + (M - 1) + (Y - 1)  + 1

        dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
        est_fun <- residuals(m) * X

        dt <- data.table(est_fun, fe1 = cluster_vec)

        setkey(dt, fe1)
        dt <- dt[, lapply(.SD, sum), by = fe1][, 2:ncol(dt), with = FALSE]

        bread <- solve(crossprod(as.matrix(X))) * N
        meat <- crossprod(as.matrix(dt)) / N

        m$se <- as.numeric(sqrt(dfc * 1/N * diag(bread %*% meat %*% bread)))
        message("SEs Clustered on First FE")
    } else {
        m$se <- m$se * sqrt(m$df.residual / DoF)
        message("SEs Not Clustered")
    }

    # Correcting degrees of freedom:
    m$df.residual <- DoF

    return(m)
}
