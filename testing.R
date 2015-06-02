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
elecnames <- as.matrix(read.delim(paste(path,"E31.asc", sep = ""), sep = "", header = F))
sample.n  <- seq_len(nrow(erpdata))
time.t    <- seq(from = 0L, by = 5L, length.out = nrow(erpdata)) 
colnames(erpdata) <- elecnames
rownames(erpdata) <- time.t
erp.epoch <- reshape2::melt(erpdata, varnames = c("Time", "Electrode"), value.name = "mV")

x.loc <- rep(x, each = length(time.t))
y.loc <- rep(y, each = length(time.t))
erp.epoch <- data.frame(x = x.loc, y = y.loc, erp.epoch)



devtools::use_package("plyr")
devtools::use_package("akima")
devtools::use_package("rgcvpack")
devtools::use_package("tripack")
devtools::use_package("colorspace")
devtools::use_package("ggplot2", "Suggests")
devtools::use_package("fields", "Suggests")
devtools::use_package("reshape2", "Suggests")
devtools::use_data(erp.epoch, overwrite = TRUE)

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
     image.plot(x = x, y = y, z = z, useRaster = TRUE, nlevel = 64))
             
with(test3$'1', 
     image.plot(x = x, y = y, z = z, useRaster = TRUE, 
                col = colorspace::diverge_hcl(11, h = c(260, 0), c = 100,
                                              l = c(30, 100), power = 1.5,
                                              fixup = TRUE)))

with(test3$'1', 
     image.plot(x = x, y = y, z = z, useRaster = TRUE, bty = "n", axes = "F", ann = "F",
                col = colorspace::diverge_hcl(11, h = c(260, 0), c = 100,
                                              l = c(30, 100), power = 1.5,
                                              fixup = TRUE)))



library(ggplot2)
new.df <- SimplifyTopographyList(test3)
ggplot(new.df, aes(x=x, y=y, fill=z)) + geom_tile()
ggplot(new.df, aes(x=x, y=y, fill=z)) + geom_raster()
ggplot(new.df, aes(x=x, y=y, fill=z)) + geom_raster(interpolate = TRUE) + theme_bw()
ggplot(new.df[!is.na(new.df$z), ], aes(x=x, y=y, fill=z)) + geom_raster(interpolate = TRUE) + theme_bw()


new.df <- SimplifyTopographyList(test3)
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


new.df <- SimplifyTopographyList(test3)
nb.colors <- 11
symm.breaks <- MakeSymmetricCuts(new.df$z[!is.na(new.df$z)], nb.colors)
angle <- seq(-pi, pi, length = 50)
#df <- data.frame(
#  x = (sin(angle) - min(new.df[!is.na(new.df$z), ]$x))/(max(new.df[!is.na(new.df$z), ]$x) - min(new.df[!is.na(new.df$z), ]$x)), 
#  y = (cos(angle) - min(new.df[!is.na(new.df$z), ]$y))/(max(new.df[!is.na(new.df$z), ]$y) - min(new.df[!is.na(new.df$z), ]$y))
#)
#elec.sphere <- tripack::circumcircle(new.df[!is.na(new.df$z), ]$x, new.df[!is.na(new.df$z), ]$y, num.touch = 2, plot = FALSE, debug = FALSE)
df <- data.frame(
  x = sin(angle) * max(abs(new.df[!is.na(new.df$z), ]$x)) * 1.015, 
  y = cos(angle) * max(abs(new.df[!is.na(new.df$z), ]$y)) * 1.015
)
##
## for eeg plots, negative is generally blue
##
ggplot(new.df[!is.na(new.df$z), ], 
       aes(x = x, y = y, fill = cut(z, breaks = symm.breaks, labels = FALSE))) + 
  geom_raster(interpolate = TRUE) + 
  scale_fill_gradient2(space = "Lab", high = scales::muted("red"), low = scales::muted("blue"), 
                       midpoint = mean(1:nb.colors),
                       limits = range(1:nb.colors),
                       breaks = seq(1, nb.colors, length = 3),
                       labels = MakeSymmetricLabels(new.df$z, multiplication.factor = 1, unit = "mv"),
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
  ) + guides(fill = guide_colorbar(nbin = nb.colors)) +
  geom_path(aes(x, y), data = df , inherit.aes = FALSE) +
  coord_fixed() 


##
## SIMULATE THE PLOTS BY THE KAYSER & TENKE GROUP
## I think they use something very similar to the Spectral scale in ColorBrewer
##
# scales::show_col(scales::brewer_pal("qual", palette = "Spectral")(11))
# scales::show_col(scales::gradient_n_pal(my.pal)(seq(0, 1, length = 30)))
#
# 7 point scale
# magenta-ish, purple-ish, red-ish, yellowish, greenish, bluish, orchid.
# R colors
# plum1, 
# c("violet", "magenta", "red", "yellow", "green", "blue", "mediumpurple")
# orchid1, magenta3, red1, yellow1, green1, blue4, mediumpurple1
#170, 170, 252
#18, 24, 188
#43, 253, 165
#254, 241, 53
#252, 34, 28
#215, 32, 202
#246, 161, 243
#
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
nb.colors <- 64
symm.breaks <- MakeSymmetricCuts(new.df$z[!is.na(new.df$z)], nb.colors)
my.pal <- scales::brewer_pal("qual", palette = "Spectral")(11)
my.pal.2 <- scales::gradient_n_pal(my.pal)(seq(0, 1, length = nb.colors))
my.pal.3 <- scales::gradient_n_pal(c("violet", "magenta", "red", "yellow", "green", "blue", "mediumpurple"))(seq(0, 1, length = nb.colors))
my.pal.4 <- scales::gradient_n_pal(c("orchid1", "magenta3", "red1", "yellow1", "green1", "blue4", "mediumpurple1"))(seq(0, 1, length = nb.colors))
my.pal.5 <- rev(jet.colors(nb.colors))
#
ggplot(new.df[!is.na(new.df$z), ], 
       aes(x = x, y = y, fill = cut(z, breaks = symm.breaks, labels = FALSE))) + 
  geom_raster(interpolate = TRUE) + 
  scale_fill_gradientn(space = "Lab", colours = rev(my.pal.5),
                       limits = range(1:nb.colors),
                       breaks = seq(1, nb.colors, length = 3),
                       labels = MakeSymmetricLabels(new.df$z, multiplication.factor = 1, unit = "mv")
                       ) +
  theme_bw() + 
  theme(line = element_blank(),
        title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        #legend.position = "right",
        plot.margin = unit(c(-4, 0, -4, -1.5), "cm")
  ) + guides(fill = guide_colorbar(nbin = nb.colors)) +
  coord_fixed() + geom_path(aes(x, y), data = df , inherit.aes = FALSE)

#
#
test.plot <- ggplot(new.df[!is.na(new.df$z), ], 
       aes(x = x, y = y, fill = cut(z, breaks = symm.breaks, labels = FALSE))) + 
  geom_raster(interpolate = TRUE) + 
  scale_fill_gradientn(space = "Lab", colours = rev(my.pal.5),
                       limits = range(1:nb.colors),
                       breaks = seq(1, nb.colors, length = 3),
                       labels = MakeSymmetricLabels(new.df$z, multiplication.factor = 1, unit = "mv")
  ) +
  theme_bw() + 
  theme(line = element_blank(),
        title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        #legend.position = "bottom",
        legend.position = "right",
        #legend.key.height = unit(3, "cm"),
        #legend.key.height = unit(0.15, "npc"),        
        legend.justification = "left",
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + guides(fill = guide_colorbar(nbin = nb.colors)) +
  coord_fixed() + geom_path(aes(x, y), data = df , inherit.aes = FALSE)


#
# inspect the plot object
# http://stackoverflow.com/questions/29749661/ggplot2-how-to-match-legend-height-to-plot-area
# http://stackoverflow.com/questions/16255579/how-can-i-make-consistent-width-plots-in-ggplot-with-legends
# http://www.noamross.net/blog/2013/11/20/formatting-plots-for-pubs.html
gp <- ggplot_gtable(ggplot_build(test.plot)) 
leg <- gtable_filter(gp, "guide-box")
leg[[1]][[1]][[1]][[1]]$heights
leg[[1]][[1]]$heights[[2]]

#
# Independent color bar function
# http://colbyimaging.duckdns.org:8080/wiki/statistics/color-bars
#




###
### Unused lines from documentations
###
### FORMATPOSION    
#' #examples
#' ## Tests
#' ##Input x missing, should throw error
#' # yvec <- 1:10
#' # FormatPosition(y = yvec)
#'
#' # Inputs are vectors
#' # xvec <- 1:10
#' # yvec <- 1:10
#' ##FormatPosition(xvec, yvec)
#' ##xvec <- 1:10
#' ##yvec <- 1:11
#' ##FormatPosition(xvec, yvec)
#' ##xvec <- 1:10
#' ##FormatPosition(xvec)
#' ##xvec <- 1:10
#' ##yvec <- letters[1:10]
#' ##FormatPosition(xvec, yvec) # Error, y is a character string
#' #
#' ##Input is a matrix
#' ##xmat1 <- matrix(1:10, nrow = 2)
#' ##FormatPosition(xmat1)
#' ##xmat2 <- matrix(1:12, ncol = 3)
#' ##FormatPosition(xmat2)
#' ##Inputs are one matrix and one vector
#' #xmat1 <- matrix(1:12, ncol = 3)
#' #yvec1 <- 1:12
#' #yvec2 <- 1:4
#' #FormatPosition(xmat1, yvec)
#' #FormatPosition(yvec1, xmat1)
#' #FormatPosition(yvec2, xmat1)
#' ## Inputs are matrices
#' #xmat1 <- matrix(1:10, nrow = 2)
#' #xmat2 <- matrix(1:12, ncol = 3)
#' #FormatPosition(xmat1, xmat2)
#' ###Input x is a list
#' #xlist <- list(a = 1:10)
#' #FormatPosition(xlist)
#' #xlist <- list(a = 1:10, b = 1:5)
#' #FormatPosition(xlist)
#' #xlist <- list(a = 1:10)
#' #ylist <- list(c = 1:9)
#' #FormatPosition(xlist, ylist)
#' ## Inputs are two lists (second should be ignored)
#' #xlist <- list(a = 1:10)
#' #ylist <- list(c = 1:9)
#' #FormatPosition(xlist, ylist)
#' ## Input is a data.frame
#' #xdf <- data.frame(a = 1:10)
#' #FormatPosition(xdf)
#' #xdf <- data.frame(a = 1:10, b = 10:1)
#' #FormatPosition(xdf)
#' ## Input is a data.frames and a vector
#' #xvec <- 1:10
#' #xdf <- data.frame(a = 1:10, b = 10:1)
#' #FormatPosition(xdf, xvec)
#' #FormatPosition(xvec, xdf)
#' ## Input is a data.frame and a vector
#' #xvec <- 1:10
#' #xdf <- data.frame(a = 1:10, b = 10:1)
#' #FormatPosition(xdf, xvec)
#' #FormatPosition(xvec, xdf) # Error, y is a list/data.frame
# 
### FORMATPOSION    
