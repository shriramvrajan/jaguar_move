## Simulations
library(gstat)
library(raster)

### Switches ===================================================================



### Functions ==================================================================

gen_grid <- function(size = 100, 
                     beta = 1, psill = 0.03, range = 10, nugget = 0) {    
    xy <- expand.grid(1:size, 1:size); names(xy) <- c("x", "y")
    model <- gstat(formula = z ~ 1, locations = ~x + y, dummy = T, beta = beta, 
                model = vgm(psill = psill, range = range, model = "Exp",
                            nugget = nugget), nmax = 20)
    out <- predict(model, newdata = xy, nsim = 1)
    if (any(out < 0)) out[out < 0] <- 0
    gridded(out) <- ~x + y; out <- raster(out)
    plot(out)

    outdf <- as.data.frame(out)
    outdf <- cbind(outdf, rowColFromCell(out, seq_len(nrow(outdf))))
    return(list(raster = out, df = outdf))
}


### Generating two random fields with different correlation structures =========

env1 <- gen_grid(size = 100, psill = 0.02, range = 10)
env2 <- gen_grid(size = 100, psill = 0.003, range = 2)

par <- c(1, -1) # preferences for env1 and env2

jag_path <- function(x1, y1, nstep, par = c(1, 1), neighb = 2) {
    if (!(x1 %in% 1:100) | !(y1 %in% 1:100)) {
        print("Jaguar out of bounds")
        return(NULL)
    }
    path <- matrix(NA, nrow = nstep, ncol = 2)
    x <- x1; y <- y1; path[1, ] <- c(x, y)
    for (i in 2:nstep) {
        print(i)
        pos <- path[i - 1, ]
        nbhd <- make_nbhd(r = env1[[1]], rdf = env1[[2]], sz = neighb,
                          i = cellFromRowCol(env1[[1]], pos[2], pos[1]))
        attract <- env1[[1]][nbhd] * par[1] + env2[[1]][nbhd] * par[2]
        print(attract)
        if (any(is.na(attract))) attract[is.na(attract)] <- 0
        if (any(attract < 0)) attract[attract < 0] <- 0
        attract <- attract / sum(attract)
        step <- sample(seq_len(length(attract)), 1, prob = attract)
        path[i, ] <- rowColFromCell(env1[[1]], nbhd[step])
    }
    return(path)
}

path1 <- jag_path(20, 20, 100, par = c(0.1, 0.3), neighb = 3)
