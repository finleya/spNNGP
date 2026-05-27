#!/usr/bin/env Rscript

set.seed(20260527)

dir.create("inst/figures", recursive = TRUE, showWarnings = FALSE)

width <- 1037
height <- 1200
scale <- 1
w <- width * scale
h <- height * scale

hex_border_size <- 1.2
hex_margin <- 1.02 + hex_border_size * 0.04
hex_half_width <- sqrt(3) / 2

x <- seq(-hex_half_width * hex_margin, hex_half_width * hex_margin, length.out = w)
y <- seq(hex_margin, -hex_margin, length.out = h)
x_mat <- matrix(rep(x, each = h), nrow = h, ncol = w)
y_mat <- matrix(rep(y, times = w), nrow = h, ncol = w)

inside_hex <- {
  ay <- abs(y_mat)
  ax <- abs(x_mat)
  (ay <= 0.5 & ax <= hex_half_width) |
    (ay > 0.5 & ay <= 1 & ax <= sqrt(3) * (1 - ay))
}

n_basis <- 70
freq <- matrix(rnorm(n_basis * 2, sd = 1.05), ncol = 2)
phase <- runif(n_basis, 0, 2 * pi)
amp <- exp(-0.34 * rowSums(freq^2)) * rnorm(n_basis)

z <- matrix(0, nrow = h, ncol = w)
for (i in seq_len(n_basis)) {
  z <- z + amp[i] * cos(2 * pi * (freq[i, 1] * x_mat + freq[i, 2] * y_mat) + phase[i])
}
z <- as.numeric(scale(z))
z[!inside_hex] <- NA_real_

palette <- c(
  "#07184f", # deep Parrish-like blue shadow
  "#0047ab", # cobalt / Parrish blue anchor
  "#2866df", # luminous cobalt blue
  "#4d3ca3", # blue violet
  "#6b2f8f", # purple glaze
  "#a7315f", # red violet
  "#c83e3e", # warm red
  "#ef6f2e", # orange-red highlight
  "#f6c45b"  # restrained warm glow
)
values <- c(-2.3, -1.45, -0.7, -0.15, 0.35, 0.85, 1.3, 1.8, 2.35)
color_fun <- grDevices::colorRamp(palette, space = "Lab")
z_scaled <- pmin(
  1,
  pmax(0, stats::approx(values, seq(0, 1, length.out = length(values)), z, rule = 2)$y)
)
z_scaled[is.na(z_scaled)] <- 0

rgb <- color_fun(z_scaled) / 255
img <- array(0, dim = c(h, w, 4))
img[, , 1] <- matrix(rgb[, 1], nrow = h, ncol = w)
img[, , 2] <- matrix(rgb[, 2], nrow = h, ncol = w)
img[, , 3] <- matrix(rgb[, 3], nrow = h, ncol = w)
img[, , 4] <- ifelse(inside_hex, 1, 0)

downsample <- function(a, factor) {
  dims <- dim(a)
  stopifnot(dims[1] %% factor == 0, dims[2] %% factor == 0)
  dim(a) <- c(factor, dims[1] / factor, factor, dims[2] / factor, dims[3])
  apply(a, c(2, 4, 5), mean)
}

if (scale > 1) {
  img <- downsample(img, scale)
}

outfile <- "inst/figures/spNNGP-hex.png"
grDevices::png(outfile, width = width, height = height, bg = "transparent", res = 300)
grid::grid.newpage()
grid::grid.raster(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)

hex_x <- c(0, -hex_half_width, -hex_half_width, 0, hex_half_width, hex_half_width)
hex_y <- c(1, 0.5, -0.5, -1, -0.5, 0.5)
vx <- 0.5 + hex_x / (2 * hex_half_width * hex_margin)
vy <- 0.5 + hex_y / (2 * hex_margin)

grid::grid.polygon(
  x = grid::unit(vx, "npc"),
  y = grid::unit(vy, "npc"),
  gp = grid::gpar(fill = NA, col = "#fff1c2", lwd = 7, linejoin = "round")
)
grid::grid.polygon(
  x = grid::unit(vx, "npc"),
  y = grid::unit(vy, "npc"),
  gp = grid::gpar(fill = NA, col = "#0b2459", lwd = 1.6, linejoin = "round")
)
grid::grid.text(
  "spNNGP",
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.48, "npc"),
  gp = grid::gpar(col = "#fff8e8", fontsize = 36, fontface = "bold", fontfamily = "sans")
)
grDevices::dev.off()

if (file.exists("Rplots.pdf")) {
  unlink("Rplots.pdf")
}
