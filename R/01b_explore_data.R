source("R/00_functions.R")

## Mar 21, 2025 ----------------------------------------------------------------
# Diel and seasonal rhythms

tr <- make_full_track(jag_id[1] %>% as.numeric)


## Jan 29, 2025 ----------------------------------------------------------------
# Tackling the timestep and max_dist unevenness.

plot1 <- function(i, v = "dt") {
    moves <- make_full_track(jag_id[i] %>% as.numeric)
    # print(range(moves$timestamp))
    # hist(moves$sl, main = i)
    hist(moves$dt, main = i)
    thresh <- mean(moves$dt, na.rm = TRUE) + 4 * sd(moves$dt, na.rm = TRUE)
    abline(v = thresh, col = "red")
    # ggplot(data = moves) +
    #     geom_line(mapping = aes(x = timestamp, y = sl), col = "red") +
    #     geom_line(mapping = aes(x = timestamp, y = dt))
}

# i1 <- c(26, 41, 43, 44, 49, 64)
# for (i in i1) plot1(i)

par(mfrow = c(4, 4))
for (i in 1:82) plot1(i)



## Jan 2025 --------------------------------------------------------------------
# Adding extra fields to metadata table

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
