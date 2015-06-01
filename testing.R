# Create 
# Topography variables
x <- c(-0.11,  0.11, -0.29, -0.15,  0,  0.15,  0.29, -0.43, -0.26, 0.26,  0.43, -0.36, -0.18,
       0, 0.18, 0.36, -0.43, -0.26, 0.26, 0.43, -0.37, -0.29, -0.15, 0, 0.15, 0.29, 0.37, -0.11, 0, 
       0.11, 0
)
y <- c(0.35, 0.35, 0.21, 0.19, 0.18, 0.19, 0.21, 0.14, 0.10, 0.10, 0.14, 0, 0, 0, 0, 0, -0.14, 
       -0.10, -0.10, -0.14, -0.27, -0.21, -0.19, -0.18, -0.19, -0.21, -0.27, -0.35, -0.36, -0.35, 0.50
)
#radius <- 0.48
#x1 <- seq(-radius, radius, , 250)
#y1 <- seq(-radius, radius, , 250)
#grid1 <- data.frame(x1, y1)
#scale.col <- diverge_hcl(11, h = c(260, 0), c = 100, l = c(30, 100), power = 1, fixup = TRUE)

# Read 
# ERP data (from package csd toolbox)
path <- "/Users/diogoalmeida/Google Drive/Rcode/" 
erpdata <- as.matrix(read.delim(paste(path,"NR_C66_trr.dat", sep = ""), sep = "", header = F))


devtools::use_package("plyr")
devtools::use_package("akima")
devtools::use_package("rgcvpack")
devtools::use_package("tripack")
devtools::use_package("colorspace")
devtools::use_package("ggplot2", "Suggests")
devtools::use_package("fields", "Suggests")

test <- InterpolatePotentials(potentials = erpdata, x.elec = x, y.elec = y, grid.grain = 40, extrapolate = TRUE, method = "tps", electrode.dim = 1)
str(test)

with(test$'50', image(x = x, y = y, z = z, useRaster = TRUE))

with(test$'50', image(x = x, y = y, z = z, useRaster = TRUE, col = colorspace::diverge_hcl(11,
                                                                               h = c(260, 0), c = 100, l = c(30, 100),
                                                                               power = 1, fixup = TRUE
)))


test2 <- InterpolatePotentials(potentials = erpdata, x.elec = x, y.elec = y, grid.grain = 400, extrapolate = TRUE, method = "tps", electrode.dim = 1)
str(test)

with(test2$'50', image(x = x, y = y, z = z, useRaster = TRUE))

with(test2$'50', image(x = x, y = y, z = z, useRaster = TRUE, col = colorspace::diverge_hcl(11,
                                                                                           h = c(260, 0), c = 100, l = c(30, 100),
                                                                                           power = 1.5, fixup = TRUE
)))




test3 <- InterpolatePotentials(potentials = erpdata[21, , drop = FALSE], x.elec = x, y.elec = y, grid.grain = 400, extrapolate = TRUE, method = "tps", electrode.dim = 1)
str(test3)

with(test3$'1', image(x = x, y = y, z = z, useRaster = TRUE))

with(test3$'1', image(x = x, y = y, z = z, useRaster = TRUE, 
                      col = colorspace::diverge_hcl(11, h = c(260, 0), c = 100, 
                                                    l = c(30, 100), power = 1.5,
                                                    fixup = TRUE)))

library(fields)
with(test3$'1', 
     image.plot(x = x, y = y, z = z, useRaster = TRUE, nlevel = 11))
             
with(test3$'1', 
     image.plot(x = x, y = y, z = z, useRaster = TRUE, 
                col = colorspace::diverge_hcl(11, h = c(260, 0), c = 100,
                                              l = c(30, 100), power = 1.5,
                                              fixup = TRUE)))


library(ggplot2)
new.df <- SimplifyTopographyList(test3)
ggplot(new.df, aes(x=x, y=y, fill=z)) + geom_tile()
ggplot(new.df, aes(x=x, y=y, fill=z)) + geom_raster()
ggplot(new.df, aes(x=x, y=y, fill=z)) + geom_raster(interpolate = TRUE) + theme_bw()
ggplot(new.df[!is.na(new.df$z), ], aes(x=x, y=y, fill=z)) + geom_raster(interpolate = TRUE) + theme_bw()



nb.colors <- 11
my_pal <- function(x) {
  diverge_hcl(x, h = c(260, 0), c = 100, l = c(30, 100), power = 1, fixup = TRUE)
}
symm.breaks <- MakeSymmetricCuts(new.df$z[!is.na(new.df$z)], nb.colors)
ggplot(new.df[!is.na(new.df$z), ], 
       aes(x = x, y = y, fill = cut(z, breaks = symm.breaks, labels = FALSE))) + 
  geom_raster(interpolate = TRUE) + 
  scale_fill_gradient2(space = "Lab", midpoint = mean(1:nb.colors),
                       limits = range(1:nb.colors),
                       breaks = seq(1, nb.colors, length = 3),
                       labels = MakeSymmetricLabels(new.df$z, multiplication.factor = 1)
  ) +
  theme_bw() + 
  theme(line = element_blank(),
        title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        #legend.position = "bottom"
        legend.position = "right"
  ) + guides(fill = guide_colorbar(nbin = nb.colors))


nb.colors <- 11
symm.breaks <- MakeSymmetricCuts(new.df$z[!is.na(new.df$z)], nb.colors)
ggplot(new.df[!is.na(new.df$z), ], 
       aes(x = x, y = y, fill = cut(z, breaks = symm.breaks, labels = FALSE))) + 
  geom_raster(interpolate = TRUE) + 
  scale_fill_gradient2(midpoint = mean(1:nb.colors),
                       limits = range(1:nb.colors),
                       breaks = seq(1, nb.colors, length = 3),
                       labels = MakeSymmetricLabels(new.df$z, multiplication.factor = 1),
                       name = "goo"
  ) +
  theme_bw() + 
  theme(line = element_blank(),
        title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        #legend.position = "bottom",
        legend.position = "right",
        legend.title = element_text("goo")
  ) + guides(fill = guide_colorbar(nbin = nb.colors), name = "goo")


