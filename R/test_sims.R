source("R/00_functions.R")  # Existing functions
source("R/classes.R")       # New classes

for (i in 1:10) {
    set.seed(i + 2001)
    env_raster <- gen_landscape(
        size = 200,
        s = 10,
        r = 100,
    )

    x0 <- ceiling(200 / 2)
    y0 <- ceiling(200 / 2)
    
    path_i <- jag_path(x0, y0, 2000, 
                    par = c(-1.5, 1.5, -0.2), 
                    neighb = 1,
                    env_raster = env_raster)
                    
    plotpdf(nm = paste0("figs/land_", i, ".pdf"))
    plot_path(path_i, int = 5)
    dev.off()
}
