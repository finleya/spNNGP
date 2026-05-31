rm(list = ls())

library(scico)
library(plot3D)

set.seed(8)

# ---- Logo controls ---------------------------------------------------------

out_file <- "spNNGP_gp_surface_plot3D.png"
img_width <- 2200
img_height <- 1800
gp_grid_n <- 96          # Number of GP sample locations; changing this changes the realization.
surface_render_n <- 150  # Render/interpolation resolution; increase this to smooth facets.

color_palette <- "roma" # Any scico palette name, e.g., "romaO", "batlow", "vikO"
draw_surface <- FALSE      # Set FALSE to hide the GP surface while keeping points/edges.
surface_alpha <- 0.01        # Surface opacity; background stays transparent.
surface_brightness <- 1 # Darkens the plot3D surface colors; 1 leaves palette unchanged.

# plot3D camera and surface controls.
# These affect spBayes_gp_surface_plot3D.png first; the hex sticker then embeds that PNG.
# The sticker stage uses a fresh temporary copy of this PNG so cached image paths
# do not make the hex sticker reuse an older surface.
plot3d_theta <- -38       # Rotate around the vertical axis. Try changing by 20-40 degrees.
plot3d_phi <- 45           # Tilt view: smaller tilts toward you; larger is more side-on.
plot3d_expand <- 0.78     # Additional plot3D z exaggeration; surface_height is usually more obvious.
plot3d_ltheta <- -35      # Light direction if plot3d_lighting is TRUE.
plot3d_lphi <- 55         # Light elevation if plot3d_lighting is TRUE.
plot3d_lighting <- FALSE  # FALSE gives a flatter/matte surface; TRUE uses controls below.
plot3d_ambient <- 0.32     # Base light level when plot3d_lighting is TRUE.
plot3d_diffuse <- 0.75     # Directional matte light strength.
plot3d_specular <- 0.05      # Highlight strength.
plot3d_exponent <- 0.18     # Highlight tightness; larger gives smaller/brighter highlights.
plot3d_sr <- 1            # Specular reflection mix; 0 is white shine, 1 tints shine by surface color.
plot3d_shade <- NA        # Use NA with lighting FALSE; set numeric values only if using lighting.

surface_height <- 1.55    # Main vertical GP scale before plotting; lower to flatten the surface.

mesh_every <- 8           # Draw every nth grid line over the surface.
mesh_color <- "#000000"
mesh_alpha <- 0.3
mesh_lwd <- 00

# NNGP point/edge layer on the surface.
point_n <- 100             # Exact number of jittered grid points.
point_grid_n_x <- NULL    # NULL chooses a near-square grid from point_n.
point_grid_n_y <- NULL    # NULL chooses a near-square grid from point_n.
point_seed <- 11          # Change this to get a different jittered point pattern.
point_jitter <- 0.28      # Fraction of grid spacing used for random x/y jitter; 0 gives a grid.
nngp_m <- 3               # Max previous-neighbor connections per point.
nngp_order_axis <- "x"    # "x", "y", "xy", or "yx"; determines predecessor ordering.
nngp_order_decreasing <- FALSE
nngp_max_edges <- Inf     # Set finite to randomly thin the plotted edges.
point_surface_z_offset <- 0.04 # Lift nodes/edges slightly above the surface so they remain visible.
edge_color <- "#EAEAEA"
edge_alpha <- 0.32
edge_lwd <- 1.0
edge_node_trim <- 0.025   # Fraction trimmed from each edge end; increase only if lines cross nodes.
edge_follow_surface <- TRUE # TRUE draws edges as curves following the GP surface.
edge_surface_steps <- 18  # Number of samples per edge when edge_follow_surface is TRUE.
point_style <- "sphere"   # "sphere" uses round highlighted nodes; "surface_sphere" is slower.
point_cex <- 1.15         # Used only when point_style is "point".
point_pch <- 21           # Used only when point_style is "point".
point_sphere_cex <- 2
point_sphere_highlight_cex <- 0.75
point_sphere_highlight_alpha <- 0.55
point_sphere_highlight_shift <- 0.018
point_sphere_radius <- 0.09
point_sphere_resolution <- 6
point_alpha <- 1
point_border_color <- "#000000"
point_border_alpha <- 0.8

sticker_file <- "spNNGP_hex_plot3D.png"
sticker_hex_fill <- "#000000"
sticker_hex_border <- "#000000"
sticker_hex_border_size <- 0.58

# Surface image placement inside the hex sticker.
# These do NOT change the 3D camera; they only place/scale the already-rendered PNG in the hex.
sticker_image_x <- 1.00
sticker_image_y <- 1.14
sticker_image_width <- 0.96  # Main control for surface size inside the hex.
sticker_image_height <- 0.74 # Less important for PNG subplots, but keep proportional.

sticker_text <- "spNNGP"
sticker_text_x <- 0.97
sticker_text_y <- 0.48
sticker_text_left <- "sp"
sticker_text_right <- "NNGP"
sticker_text_color_offset <- 0.28 # Distance from palette center; must be between 0 and 0.5.
sticker_text_join_x <- NULL        # NULL centers the combined two-color word.
sticker_text_width <- 1.0          # Approximate full word width in sticker coordinates.
sticker_text_family <- "Latin Modern Roman"
sticker_text_face <- "plain"
sticker_text_size <- 45
sticker_panel_pad <- 0.008    # Extra plotting-room so the hex border is not clipped.
sticker_padding_px <- 5       # Transparent padding added evenly around final sticker.

font_files <- c(
  "Latin Modern Roman" = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-regular.otf",
  "Liberation Mono" = "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf"
)

if (requireNamespace("sysfonts", quietly = TRUE) &&
    requireNamespace("showtext", quietly = TRUE) &&
    sticker_text_family %in% names(font_files) &&
    file.exists(font_files[[sticker_text_family]])) {
  sysfonts::font_add(sticker_text_family, regular = font_files[[sticker_text_family]])
  showtext::showtext_auto(TRUE)
}

gp_grid <- function(n = 96, sigma.sq = 1, range = 0.26, nugget = 1e-8) {
  x <- seq(-1, 1, length.out = n)
  y <- seq(-1, 1, length.out = n)

  # A squared-exponential GP on a regular grid can be sampled efficiently
  # with separable covariance matrices, avoiding a huge dense 2-D covariance.
  cov_1d <- function(u) {
    D <- as.matrix(dist(u))
    sigma.sq * exp(-0.5 * (D / range)^2)
  }

  Lx <- chol(cov_1d(x) + diag(nugget, n))
  Ly <- chol(cov_1d(y) + diag(nugget, n))
  z <- t(Lx) %*% matrix(rnorm(n * n), n, n) %*% Ly

  z <- z - min(z)
  z <- z / max(z)
  z <- z^1.18
  z <- z - 0.45
  z <- z * surface_height

  list(x = x, y = y, z = z)
}

resample_surface <- function(gp, n = surface_render_n) {
  if (is.null(n) || n <= length(gp$x)) {
    return(gp)
  }

  x_new <- seq(min(gp$x), max(gp$x), length.out = n)
  y_new <- seq(min(gp$y), max(gp$y), length.out = n)

  # Interpolate the already sampled GP. This increases render resolution
  # without drawing a new random surface.
  z_x <- apply(
    gp$z,
    2,
    function(z_col) {
      approx(gp$x, z_col, xout = x_new, ties = "ordered")$y
    }
  )
  z_new <- t(apply(
    z_x,
    1,
    function(z_row) {
      approx(gp$y, z_row, xout = y_new, ties = "ordered")$y
    }
  ))

  list(x = x_new, y = y_new, z = z_new)
}

surface_cols <- function(z, n_col = 256) {
  pal <- scico(n_col, palette = color_palette)
  z_scaled <- (z - min(z)) / diff(range(z))
  pal[pmax(1, pmin(n_col, floor(z_scaled * (n_col - 1)) + 1))]
}

surface_value_cols <- function(values, z_range, n_col = 256) {
  pal <- scico(n_col, palette = color_palette)
  scaled <- (values - z_range[1]) / diff(z_range)
  pal[pmax(1, pmin(n_col, floor(scaled * (n_col - 1)) + 1))]
}

surface_at_points <- function(gp, xout, yout) {
  ix <- findInterval(xout, gp$x, all.inside = TRUE)
  iy <- findInterval(yout, gp$y, all.inside = TRUE)
  ix <- pmin(ix, length(gp$x) - 1)
  iy <- pmin(iy, length(gp$y) - 1)

  x0 <- gp$x[ix]
  x1 <- gp$x[ix + 1]
  y0 <- gp$y[iy]
  y1 <- gp$y[iy + 1]
  tx <- (xout - x0) / (x1 - x0)
  ty <- (yout - y0) / (y1 - y0)

  z00 <- gp$z[cbind(ix, iy)]
  z10 <- gp$z[cbind(ix + 1, iy)]
  z01 <- gp$z[cbind(ix, iy + 1)]
  z11 <- gp$z[cbind(ix + 1, iy + 1)]

  (1 - tx) * (1 - ty) * z00 +
    tx * (1 - ty) * z10 +
    (1 - tx) * ty * z01 +
    tx * ty * z11
}

make_point_locations <- function(gp) {
  if (!is.null(point_seed)) {
    set.seed(point_seed)
  }

  n_points <- max(1, as.integer(point_n))
  n_x <- if (is.null(point_grid_n_x)) ceiling(sqrt(n_points)) else point_grid_n_x
  n_y <- if (is.null(point_grid_n_y)) ceiling(n_points / n_x) else point_grid_n_y
  n_x <- max(1, as.integer(n_x))
  n_y <- max(1, as.integer(n_y))

  x_grid <- seq(min(gp$x), max(gp$x), length.out = n_x)
  y_grid <- seq(min(gp$y), max(gp$y), length.out = n_y)
  pts <- expand.grid(x = x_grid, y = y_grid)
  if (nrow(pts) > n_points) {
    pts <- pts[seq_len(n_points), , drop = FALSE]
  }

  x_spacing <- if (length(x_grid) > 1) min(diff(x_grid)) else diff(range(gp$x))
  y_spacing <- if (length(y_grid) > 1) min(diff(y_grid)) else diff(range(gp$y))
  jitter <- max(0, point_jitter)

  pts$x <- pts$x + runif(nrow(pts), -jitter * x_spacing, jitter * x_spacing)
  pts$y <- pts$y + runif(nrow(pts), -jitter * y_spacing, jitter * y_spacing)
  pts$x <- pmin(max(gp$x), pmax(min(gp$x), pts$x))
  pts$y <- pmin(max(gp$y), pmax(min(gp$y), pts$y))

  pts$surface_z <- surface_at_points(gp, pts$x, pts$y)
  pts$z <- pts$surface_z + point_surface_z_offset
  pts$col <- surface_value_cols(pts$surface_z, range(gp$z))
  pts
}

order_points_for_nngp <- function(pts) {
  if (nngp_order_axis == "x") {
    ord <- order(pts$x, pts$y, decreasing = nngp_order_decreasing)
  } else if (nngp_order_axis == "y") {
    ord <- order(pts$y, pts$x, decreasing = nngp_order_decreasing)
  } else if (nngp_order_axis == "xy") {
    ord <- order(pts$x + pts$y, pts$x, decreasing = nngp_order_decreasing)
  } else if (nngp_order_axis == "yx") {
    ord <- order(pts$y + pts$x, pts$y, decreasing = nngp_order_decreasing)
  } else {
    stop("nngp_order_axis must be one of 'x', 'y', 'xy', or 'yx'.")
  }

  pts[ord, , drop = FALSE]
}

make_nngp_edges <- function(pts, m = nngp_m) {
  pts_ordered <- order_points_for_nngp(pts)

  if (requireNamespace("spNNGP", quietly = TRUE)) {
    nngp_ns <- asNamespace("spNNGP")
    nn <- get("mkNNIndx", nngp_ns)(
      coords = as.matrix(pts_ordered[, c("x", "y")]),
      m = m,
      n.omp.threads = 1
    )
    nn_list <- get("mk.n.indx.list", nngp_ns)(nn$nnIndx, nrow(pts_ordered), m)
  } else {
    nn_list <- vector("list", nrow(pts_ordered))
    nn_list[[1]] <- NA_integer_
    for (i in seq_len(nrow(pts_ordered))) {
      if (i == 1) {
        next
      }

      previous <- seq_len(i - 1)
      d2 <- (pts_ordered$x[previous] - pts_ordered$x[i])^2 +
        (pts_ordered$y[previous] - pts_ordered$y[i])^2
      nn_list[[i]] <- previous[order(d2)][seq_len(min(m, length(previous)))]
    }
  }

  edges <- vector("list", max(0, nrow(pts_ordered) - 1))
  for (i in seq_len(nrow(pts_ordered))) {
    if (i == 1) {
      next
    }

    keep <- nn_list[[i]]
    keep <- keep[!is.na(keep)]

    edges[[i - 1]] <- data.frame(
      x0 = pts_ordered$x[i],
      y0 = pts_ordered$y[i],
      z0 = pts_ordered$z[i],
      x1 = pts_ordered$x[keep],
      y1 = pts_ordered$y[keep],
      z1 = pts_ordered$z[keep]
    )
  }

  edges <- do.call(rbind, edges)
  if (is.finite(nngp_max_edges) && nrow(edges) > nngp_max_edges) {
    edges <- edges[sort(sample(seq_len(nrow(edges)), nngp_max_edges)), , drop = FALSE]
  }

  trim <- max(0, min(0.45, edge_node_trim))
  if (trim > 0 && nrow(edges) > 0) {
    dx <- edges$x1 - edges$x0
    dy <- edges$y1 - edges$y0
    dz <- edges$z1 - edges$z0
    edges$x0 <- edges$x0 + trim * dx
    edges$y0 <- edges$y0 + trim * dy
    edges$z0 <- edges$z0 + trim * dz
    edges$x1 <- edges$x1 - trim * dx
    edges$y1 <- edges$y1 - trim * dy
    edges$z1 <- edges$z1 - trim * dz
  }

  edges
}

draw_nngp_edges <- function(gp, edges) {
  edge_col <- grDevices::adjustcolor(edge_color, alpha.f = edge_alpha)

  if (isTRUE(edge_follow_surface)) {
    n_steps <- max(2, as.integer(edge_surface_steps))
    for (i in seq_len(nrow(edges))) {
      x <- seq(edges$x0[i], edges$x1[i], length.out = n_steps)
      y <- seq(edges$y0[i], edges$y1[i], length.out = n_steps)
      z <- surface_at_points(gp, x, y) + point_surface_z_offset
      plot3D::lines3D(
        x = x,
        y = y,
        z = z,
        add = TRUE,
        colvar = NULL,
        col = edge_col,
        lwd = edge_lwd,
        colkey = FALSE
      )
    }
    return(invisible(TRUE))
  }

  plot3D::segments3D(
    x0 = edges$x0,
    y0 = edges$y0,
    z0 = edges$z0,
    x1 = edges$x1,
    y1 = edges$y1,
    z1 = edges$z1,
    add = TRUE,
    colvar = NULL,
    col = edge_col,
    lwd = edge_lwd,
    colkey = FALSE
  )

  invisible(TRUE)
}

draw_sphere_points <- function(pts) {
  if (identical(point_style, "sphere")) {
    plot3D::points3D(
      x = pts$x,
      y = pts$y,
      z = pts$z,
      add = TRUE,
      colvar = NULL,
      pch = 21,
      cex = point_sphere_cex,
      col = grDevices::adjustcolor(point_border_color, alpha.f = point_border_alpha),
      bg = grDevices::adjustcolor(pts$col, alpha.f = point_alpha),
      colkey = FALSE
    )
    plot3D::points3D(
      x = pts$x - point_sphere_highlight_shift,
      y = pts$y + point_sphere_highlight_shift,
      z = pts$z + point_sphere_highlight_shift,
      add = TRUE,
      colvar = NULL,
      pch = 16,
      cex = point_sphere_highlight_cex,
      col = grDevices::adjustcolor("#FFFFFF", alpha.f = point_sphere_highlight_alpha),
      colkey = FALSE
    )
    return(invisible(TRUE))
  }

  if (point_sphere_radius <= 0 || point_alpha <= 0) {
    return(invisible(FALSE))
  }

  u <- seq(0, 2 * pi, length.out = point_sphere_resolution)
  v <- seq(0, pi, length.out = point_sphere_resolution)
  uv <- expand.grid(u = u, v = v)
  sx <- matrix(cos(uv$u) * sin(uv$v), nrow = length(u), ncol = length(v))
  sy <- matrix(sin(uv$u) * sin(uv$v), nrow = length(u), ncol = length(v))
  sz <- matrix(cos(uv$v), nrow = length(u), ncol = length(v))
  colvar <- matrix(1, nrow = length(u), ncol = length(v))

  for (i in seq_len(nrow(pts))) {
    plot3D::surf3D(
      x = pts$x[i] + point_sphere_radius * sx,
      y = pts$y[i] + point_sphere_radius * sy,
      z = pts$z[i] + point_sphere_radius * sz,
      colvar = colvar,
      col = grDevices::adjustcolor(pts$col[i], alpha.f = point_alpha),
      border = grDevices::adjustcolor(point_border_color, alpha.f = point_border_alpha),
      facets = TRUE,
      add = TRUE,
      colkey = FALSE,
      lighting = FALSE,
      shade = NA
    )
  }

  invisible(TRUE)
}

darken_cols <- function(cols, brightness = surface_brightness) {
  rgb <- matrix(grDevices::col2rgb(cols), nrow = 3) / 255
  rgb <- matrix(pmax(0, pmin(1, rgb * brightness)), nrow = 3)
  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ])
}

darkest_right_palette_color <- function(n_col = 256) {
  pal <- scico(n_col, palette = color_palette)
  right <- pal[seq(floor(n_col / 2), n_col)]
  rgb <- grDevices::col2rgb(right)
  luminance <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  right[which.min(luminance)]
}

plot3d_lighting_control <- function() {
  if (!isTRUE(plot3d_lighting)) {
    return(plot3d_lighting)
  }

  list(
    type = "light",
    ambient = plot3d_ambient,
    diffuse = plot3d_diffuse,
    specular = plot3d_specular,
    exponent = plot3d_exponent,
    sr = plot3d_sr,
    alpha = surface_alpha
  )
}

sticker_text_colors <- function(n_col = 256) {
  offset <- max(0, min(0.5, sticker_text_color_offset))
  pal <- scico(n_col, palette = color_palette)
  positions <- c(0.5 - offset, 0.5 + offset)
  idx <- pmax(1, pmin(n_col, round(positions * (n_col - 1)) + 1))
  pal[idx]
}

draw_with_plot3D <- function(gp, file) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(file, width = img_width, height = img_height,
                  units = "px", background = "transparent", res = 300)
  } else {
    png(file, width = img_width, height = img_height,
        bg = "transparent", res = 300)
  }
  on.exit(dev.off(), add = TRUE)

  par(mar = rep(0, 4), bg = NA)

  z <- gp$z
  pal <- darken_cols(scico(512, palette = color_palette))
  surface_col <- if (isTRUE(draw_surface)) pal else "transparent"
  surface_plot_alpha <- if (isTRUE(draw_surface)) surface_alpha else 0
  plot3D::persp3D(
    x = gp$x,
    y = gp$y,
    z = z,
    colvar = z,
    theta = plot3d_theta,
    phi = plot3d_phi,
    scale = FALSE,
    expand = plot3d_expand,
    col = surface_col,
    border = NA,
    shade = plot3d_shade,
    lighting = plot3d_lighting_control(),
    alpha = surface_plot_alpha,
    ltheta = plot3d_ltheta,
    lphi = plot3d_lphi,
    bty = "n",
    axes = FALSE,
    ticktype = "simple",
    colkey = FALSE,
    xlab = NA,
    ylab = NA,
    zlab = NA
  )

  if (isTRUE(draw_surface) && mesh_every > 0 && mesh_lwd > 0 && mesh_alpha > 0) {
    mesh_col <- grDevices::adjustcolor(mesh_color, alpha.f = mesh_alpha)
    mesh_idx_x <- seq(1, length(gp$x), by = mesh_every)
    mesh_idx_y <- seq(1, length(gp$y), by = mesh_every)

    for (j in mesh_idx_y) {
      plot3D::lines3D(
        x = gp$x,
        y = rep(gp$y[j], length(gp$x)),
        z = z[, j],
        add = TRUE,
        colvar = NULL,
        col = mesh_col,
        lwd = mesh_lwd,
        colkey = FALSE
      )
    }
    for (i in mesh_idx_x) {
      plot3D::lines3D(
        x = rep(gp$x[i], length(gp$y)),
        y = gp$y,
        z = z[i, ],
        add = TRUE,
        colvar = NULL,
        col = mesh_col,
        lwd = mesh_lwd,
        colkey = FALSE
      )
    }
  }

  pts <- make_point_locations(gp)
  edges <- make_nngp_edges(pts)

  # Draw edges before nodes so the sphere-like points cover line endpoints.
  if (!is.null(edges) && nrow(edges) > 0 && edge_lwd > 0 && edge_alpha > 0) {
    draw_nngp_edges(gp, edges)
  }

  # Draw nodes last so they sit visually on top of the NNGP edges.
  if (nrow(pts) > 0 && point_style %in% c("sphere", "surface_sphere")) {
    draw_sphere_points(pts)
  } else if (nrow(pts) > 0 && point_cex > 0 && point_alpha > 0) {
    plot3D::points3D(
      x = pts$x,
      y = pts$y,
      z = pts$z,
      add = TRUE,
      colvar = NULL,
      pch = point_pch,
      cex = point_cex,
      col = grDevices::adjustcolor(point_border_color, alpha.f = point_border_alpha),
      bg = grDevices::adjustcolor(pts$col, alpha.f = point_alpha),
      colkey = FALSE
    )
  }

  invisible(TRUE)
}

save_sticker_plot <- function(sticker_plot, output_file, panel_pad) {
  center <- 1
  radius <- 1
  half_width <- sqrt(3) / 2 * radius
  border_room <- sticker_hex_border_size * 0.04
  x_pad <- half_width * panel_pad + border_room
  y_pad <- radius * panel_pad + border_room

  sticker_plot <- sticker_plot +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      limits = c(center - half_width - x_pad, center + half_width + x_pad)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      limits = c(center - radius - y_pad, center + radius + y_pad)
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "lines"))

  ggplot2::ggsave(
    filename = output_file,
    plot = sticker_plot,
    width = 43.9,
    height = 50.8,
    units = "mm",
    bg = "transparent",
    dpi = 600
  )
}

add_sticker_padding <- function(input_file, output_file, padding_px) {
  if (!requireNamespace("magick", quietly = TRUE) || padding_px <= 0) {
    file.copy(input_file, output_file, overwrite = TRUE)
    return(invisible(output_file))
  }

  img <- magick::image_read(input_file)
  info <- magick::image_info(img)
  geometry <- sprintf("%dx%d", info$width + 2 * padding_px, info$height + 2 * padding_px)
  img <- magick::image_extent(img, geometry = geometry, gravity = "center", color = "none")
  magick::image_write(img, path = output_file, format = "png")

  invisible(output_file)
}

add_sticker_text <- function(sticker_plot) {
  cols <- sticker_text_colors()
  split <- if (is.null(sticker_text_join_x)) {
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      widths <- systemfonts::string_width(
        c(sticker_text_left, sticker_text_right),
        family = sticker_text_family,
        size = sticker_text_size
      )
      widths[1] / sum(widths)
    } else {
      nchar(sticker_text_left) / nchar(paste0(sticker_text_left, sticker_text_right))
    }
  } else {
    0.5
  }
  join_x <- if (is.null(sticker_text_join_x)) {
    sticker_text_x + (split - 0.5) * sticker_text_width
  } else {
    sticker_text_join_x
  }
  sticker_plot +
    ggplot2::annotate(
      "text",
      x = join_x,
      y = sticker_text_y,
      label = sticker_text_left,
      hjust = 1,
      color = cols[1],
      family = sticker_text_family,
      fontface = sticker_text_face,
      size = sticker_text_size
    ) +
    ggplot2::annotate(
      "text",
      x = join_x,
      y = sticker_text_y,
      label = sticker_text_right,
      hjust = 0,
      color = cols[2],
      family = sticker_text_family,
      fontface = sticker_text_face,
      size = sticker_text_size
    )
}

gp <- resample_surface(gp_grid(gp_grid_n), surface_render_n)

if (is.null(sticker_hex_border)) {
  sticker_hex_border <- darkest_right_palette_color()
}

draw_with_plot3D(gp, out_file)

message("Wrote surface PNG: ", normalizePath(out_file, mustWork = FALSE))

sticker_surface_tmp <- tempfile(fileext = ".png")
sticker_tmp <- tempfile(fileext = ".png")
sticker_panel_tmp <- tempfile(fileext = ".png")

if (!file.copy(out_file, sticker_surface_tmp, overwrite = TRUE)) {
  stop("Could not copy the rendered surface PNG for the sticker.")
}

sticker_plot <- hexSticker::sticker(
  subplot = sticker_surface_tmp,
  s_x = sticker_image_x,
  s_y = sticker_image_y,
  s_width = sticker_image_width,
  s_height = sticker_image_height,
  package = "",
  p_x = sticker_text_x,
  p_y = sticker_text_y,
  p_color = "transparent",
  p_family = sticker_text_family,
  p_fontface = sticker_text_face,
  p_size = 0,
  h_fill = sticker_hex_fill,
  h_color = sticker_hex_border,
  h_size = sticker_hex_border_size,
  white_around_sticker = FALSE,
  filename = sticker_tmp,
  dpi = 600
)

sticker_plot <- add_sticker_text(sticker_plot)
save_sticker_plot(sticker_plot, sticker_panel_tmp, sticker_panel_pad)
add_sticker_padding(sticker_panel_tmp, sticker_file, sticker_padding_px)

message("Wrote hex sticker PNG: ", normalizePath(sticker_file, mustWork = FALSE))

system(paste0("rm ", out_file))
