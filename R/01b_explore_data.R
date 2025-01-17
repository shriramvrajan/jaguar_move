source("R/00_functions.R")

meta <- jag_meta[, c(1, 2, 3, 4, 10, 11, 12)]
## Extra variables
meta$nmove <- sapply(1:njag, function(x) {
    length(which(jag_move$ID == as.numeric(jag_id[x])))
})
meta$ndays <- sapply(1:njag, function(x) {
    moves <- jag_move[ID == as.numeric(jag_id[x])]
    dates <- sort(as.Date(sapply(moves$timestamp, function(dt) {
            strsplit(as.character(dt), " ")[[1]][1]
        }), format = "%m/%d/%y"))
    return(as.numeric(difftime(dates[length(dates)], dates[1])))
})
meta$meandist <- tapply(jag_move$dist, jag_move$ID, function(x) mean(x, na.rm = TRUE))
meta$totdist <- tapply(jag_move$dist, jag_move$ID, function(x) sum(x, na.rm = TRUE))

meta$osize <- sapply(1:njag, function(x) {
    name <- paste0("sizeout_", x, ".rds")
    if (file.exists(paste0("data/output/", name))) {
        readRDS(paste0("data/output/", name)) %>% 
            gsub(" Mb", "", .) %>% as.numeric()
    } else {
        return(NA)
    }
})
