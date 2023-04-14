source("R/00_basics.R")

jagmeta_br <- jag_meta[jag_id[[1]], ]

## Load parameters and likelihood
load_output <- function(name) {
    dir <- paste0("data/output/", name)
    # ll_files <- list.files(dir)[grep("likelihood_", list.files(dir))]
    # par_files <- list.files(dir)[grep("par_out_", list.files(dir))]
    ll_files <- paste0("likelihood_", 1:njag, ".RDS")
    par_files <- paste0("par_out_", 1:njag, ".RDS")
    ll <- load_if_exists(ll_files, dir)
    par <- load_if_exists(par_files, dir)
    return(list(unlist(ll), par))
}

par_to_df <- function(par) {
    df <- do.call(rbind, lapply(par, function(x) {
        print(x[[1]])
    }))
}

runk <- load_output("LL_K")
runrw <- load_output("LL_RW")
runrwh <- load_output("LL_RWH")
llk <- -runk[[1]]; llrw <- -runrw[[1]]; llrwh <- -runrwh[[1]]
park <- runk[[2]]; parrw <- runrw[[2]]; parrwh <- runrwh[[2]]

null_files <- paste0("null_", 1:njag, ".RDS")
ll0 <- -unlist(load_if_exists(null_files, "data/output/LL_null"))

dfk <- par_to_df(park)
dfrw <- par_to_df(parrw)

nmove <- sapply(1:njag, function(x) length(which(jag_move$ID == 
                                                 as.numeric(jag_id[x]))))
ndays <- sapply(1:njag, function(x) {
    moves <- jag_move[ID == as.numeric(jag_id[x])]
    dates <- sort(as.Date(sapply(moves$timestamp, function(dt) {
            strsplit(as.character(dt), " ")[[1]][1]
        }), format = "%m/%d/%y"))
    return(as.numeric(difftime(dates[length(dates)], dates[1])))
})
meandist <- tapply(jag_move$dist, jag_move$ID, function(x) mean(x, na.rm = T))
totdist <- tapply(jag_move$dist, jag_move$ID, function(x) sum(x, na.rm = T))

res <- data.frame(id = jag_id,
                  sex = as.factor(jagmeta_br$Sex),
                  age = as.numeric(jagmeta_br$Estimated.Age),
                  weight = as.numeric(jagmeta_br$Weight),
                  nmove = nmove,
                  ndays = ndays,
                  meandist = meandist,
                  distpd = totdist / ndays,
                  movespd = nmove / ndays,
                  llk = llk,
                  llrw = llrw,
                  llrwh = llrwh,
                  ll0 = ll0,
                  dek = (1 - llk / ll0),
                  derw = (1 - llrw / ll0),
                  derwh = (1 - llrwh / ll0))

res <- res[-which(is.infinite(res$ll0)), ]
res <- res[-which(res$nmove < 100), ]
res <- res[order(res$derw - res$dek), ]

barplot(res$derw, col = rgb(0, 0, 1, 0.3), border = NA)
barplot(res$dek, col = rgb(1, 0, 0, 0.3), border = NA, add = T)
barplot(res$derwh, col = rgb(0, 1, 0, 0.3), border = NA, add = T)
## are sdms too broad of a scale? 