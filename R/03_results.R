source("R/00_functions.R")

jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

runk <- load_output("LL_K")
runrw <- load_output("LL_RW")
runh <- load_output("LL_H")
runrwh <- load_output("LL_RWH")
llk <- -runk[[1]]; llh <- -runh[[1]]; llrw <- -runrw[[1]]; llrwh <- -runrwh[[1]]
park <- runk[[2]]; parh <- runh[[2]]; parrw <- runrw[[2]]; parrwh <- runrwh[[2]]

### AIC calculation 
### AIC = 2k - 2ln(L)
aic_k <- 2 * 6 - 2 * llk
aic_h <- 2 * 1 - 2 * llh
aic_rw <- 2 * 6 - 2 * llrw
aic_rwh <- 2 * 7 - 2 * llrwh

null_files <- paste0("null_", 1:njag, ".RDS")
ll0 <- -unlist(load_if_exists(null_files, "data/output/LL_null"))

dfk <- par_to_df(park)
dfrw <- par_to_df(parrw)
dfh <- par_to_df(parh)
dfrwh <- par_to_df(parrwh)

nmove <- sapply(1:njag, function(x) {
                  length(which(jag_move$ID == as.numeric(jag_id[x])))
                  })
ndays <- sapply(1:njag, function(x) {
    moves <- jag_move[ID == as.numeric(jag_id[x])]
    dates <- sort(as.Date(sapply(moves$timestamp, function(dt) {
            strsplit(as.character(dt), " ")[[1]][1]
        }), format = "%m/%d/%y"))
    return(as.numeric(difftime(dates[length(dates)], dates[1])))
})
meandist <- tapply(jag_move$dist, jag_move$ID, function(x) mean(x, na.rm = T))
totdist <- tapply(jag_move$dist, jag_move$ID, function(x) sum(x, na.rm = T))

res <- data.table(id = jag_id,
                  sex = as.factor(jagmeta_br$Sex),
                  age = as.numeric(jagmeta_br$Estimated.Age),
                  weight = as.numeric(jagmeta_br$Weight),
                  bio = as.factor(jagmeta_br$biome),
                  nmove = nmove,
                  ndays = ndays,
                  meandist = meandist,
                  distpd = totdist / ndays,
                  movespd = nmove / ndays,
                  llk = llk,
                  llrw = llrw,
                  llh = llh,
                  llrwh = llrwh,
                  ll0 = ll0,
                  dek = (1 - llk / ll0),
                  derw = (1 - llrw / ll0),
                  deh = (1 - llh / ll0),
                  derwh = (1 - llrwh / ll0),
                  aic_k = aic_k,
                  aic_rw = aic_rw,
                  aic_h = aic_h,
                  aic_rwh = aic_rwh)

res <- cbind(res, do.call(rbind, sapply(1:njag, function(i) {
    parrwh[[i]][[1]]
})))
names(res)[24:30] <- paste0("p", 1:7)

res <- res[-which(is.infinite(res$ll0)), ]
# res <- res[-which(res$nmove < 100), ]

msg("Results compiled")

if (FALSE) {
   hist(aic_k, 100, col = rgb(1, 0, 0, 0.4), border = NA)
    hist(aic_rw, 100, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA)
    hist(aic_rwh, 100, col = rgb(0, 1, 0, 0.4), add = TRUE, border = NA)

    plot(aic_k, pch = 19)
    points(aic_rw, pch = 19, col = "blue")
    points(aic_rwh, pch = 19, col = "green")
    points(aic_h, pch = 19, col = "red")

    barplot(aic_k, col = rgb(1, 0, 0, 0.4), border = NA)
    barplot(aic_rw, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA)
    barplot(aic_rwh, col = rgb(0, 1, 0, 0.4), add = TRUE, border = NA)
    barplot(aic_h, col = rgb(0.5, 0.5, 0, 0.4), add = TRUE, border = NA)

}