# =============================================================================
# Full corrected script for bivariate mapping (Thailand)
# Based on: https://www.datawim.com/post/creating-professional-bivariate-maps-in-r/
# Modifications: fixed caption margin, removed invalid params, adjusted ggsave
# =============================================================================

# Install and load packages ---------------------------------------------------
#pak::pkg_install("mikejohnson51/climateR")

pacman::p_load(
  tidyverse, sf, terra, geobounds, climateR, geodata, 
  biscale, ggplot2, dplyr, cowplot, countrycode, httr, pkgconfig, elevatr
)

# Set cache directory for geoBoundaries ---------------------------------------
newdir <- file.path(tempdir(), "/geoboundvignette")
gb_set_cache_dir(newdir)  # optionally with install = TRUE for persistence

# Check current cache
current <- gb_detect_cache_dir()
print(current)
print(current == gb_detect_cache_dir())

# Load Thailand boundary data ------------------------------------------------
# ADM0 (country)

countrycode(
  c("Thailand", "Japan", "Pakistan", "Srilanka", "Myanmar", "Lao", "Vietnam", "Nepal", "China"),
  origin = "country.name",
  destination = "iso3c"
)

# Select country boundary data
iso3c <- "NPL"

roi_simp0 <- gb_get_adm0(country = iso3c, simplified = TRUE) |> 
  mutate(res = "Simplified")

print(object.size(roi_simp0), units = "Mb")

# ADM1 (provinces/regions)
roi_simp1 <- gb_get_adm1(country = iso3c, simplified = TRUE) |> 
  mutate(res = "Simplified")

print(object.size(roi_simp1), units = "Mb")

# Quick plot to inspect geometry
plot(st_geometry(roi_simp1), border = "gray")
plot(st_geometry(roi_simp0), add = TRUE, lwd = 2)   # base plot uses lwd, not size

roi_simp0 <- st_transform(roi_simp0, 4326)

# After loading roi_simp1 (Japan ADM1)
bbox <- st_bbox(roi_simp0)

# Temperature data ------------------------------------------------------------
# Japan's wide bounding box causes "DAP DATADDS packet is apparently too short"
# even with 5-year chunks.  Fix: yearly chunks (12 layers each) + retry with
# exponential back-off so transient server hiccups don't abort the whole run.
httr::set_config(httr::timeout(600))
options(terra.remote = "https")

# Retry wrapper: tries up to max_tries times with increasing sleep between attempts
safe_terraclim <- function(AOI, varname, startDate, endDate, max_tries = 5) {
  for (attempt in seq_len(max_tries)) {
    result <- tryCatch(
      getTerraClim(AOI = AOI, varname = varname, startDate = startDate, endDate = endDate),
      error = function(e) {
        message(sprintf("  [%s] attempt %d failed: %s", varname, attempt, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(result)) return(result)
    Sys.sleep(10 * attempt)   # back-off: 10 s, 20 s, 30 s, 40 s
  }
  stop("All ", max_tries, " attempts failed for ", varname, " ", startDate, " to ", endDate)
}

years     <- 1995:2025
tmax_list <- vector("list", length(years))
tmin_list <- vector("list", length(years))

for (i in seq_along(years)) {
  s <- paste0(years[i], "-01-01")
  e <- paste0(years[i], "-12-31")
  message("Downloading tmax ", years[i])
  tmax_list[[i]] <- safe_terraclim(roi_simp0, "tmax", s, e)$tmax
  message("Downloading tmin ", years[i])
  tmin_list[[i]] <- safe_terraclim(roi_simp0, "tmin", s, e)$tmin
}

tmax_all <- do.call(c, tmax_list)
tmin_all <- do.call(c, tmin_list)

# Monthly mean temperature
mean_temp_monthly <- (tmax_all + tmin_all) / 2
annual_mean_temp  <- mean(mean_temp_monthly, na.rm = TRUE)

# Project and resample to ~1 km resolution
crs1 <- crs(roi_simp1)
annual_mean_temp <- terra::project(annual_mean_temp, "EPSG:4326")
new_res <- 0.025 # degree
new_raster <- rast(ext(annual_mean_temp), resolution = new_res, crs = crs(annual_mean_temp))
annual_mean_temp <- resample(annual_mean_temp, new_raster, method = "bilinear")
roi_temp <- terra::crop(annual_mean_temp, roi_simp1, mask = TRUE)

# Precipitation data ----------------------------------------------------------
ppt_list <- vector("list", length(years))

for (i in seq_along(years)) {
  s <- paste0(years[i], "-01-01")
  e <- paste0(years[i], "-12-31")
  message("Downloading ppt ", years[i])
  ppt_list[[i]] <- safe_terraclim(roi_simp0, "ppt", s, e)$ppt
}

ppt_all <- do.call(c, ppt_list)
annual_mean_ppt <- mean(ppt_all, na.rm = TRUE)
annual_mean_ppt <- terra::project(annual_mean_ppt, "EPSG:4326")
new_raster <- rast(ext(annual_mean_ppt), resolution = new_res, crs = crs(annual_mean_ppt))
annual_mean_ppt <- resample(annual_mean_ppt, new_raster, method = "bilinear")
roi_ppt <- terra::crop(annual_mean_ppt, roi_simp1, mask = TRUE)

# Stack rasters ---------------------------------------------------------------
# Align roi_ppt to roi_temp grid to ensure identical extents before stacking
roi_ppt <- terra::resample(roi_ppt, roi_temp, method = "bilinear")
temp_ppt <- c(roi_temp, roi_ppt)
names(temp_ppt) <- c("temp", "ppt")

# Convert to data frame -------------------------------------------------------

temp_ppt_df <- temp_ppt |> 
  project(roi_simp1) |> 
  as.data.frame(xy = TRUE)
head(temp_ppt_df)

temp_ppt_df <- temp_ppt_df |> drop_na(temp, ppt)

# Bivariate classification ----------------------------------------------------
data <- bi_class(temp_ppt_df,
                 x = temp, 
                 y = ppt, 
                 style = "quantile", dim = 4)

# Check class distribution
data |> 
  count(bi_class) |> 
  ggplot(aes(x = bi_class, y = n)) +
  geom_col() +
  labs(title = "Distribution of Bivariate Classes", 
       x = "Bivariate Class", y = "Frequency")

# =============================================================================
# MAPPING – CORRECTED SECTION
# =============================================================================

pallet <- "BlueOr"

# Main map
map <- ggplot() +
  theme_void(base_size = 14) +
  xlim(bbox$xmin, bbox$xmax) +
  ylim(bbox$ymin, bbox$ymax) +
  geom_raster(data = data, 
              mapping = aes(x = x, y = y, fill = bi_class), 
              show.legend = FALSE) +
  bi_scale_fill(pal = pallet, dim = 4, flip_axes = FALSE, rotate_pal = FALSE) +
  geom_sf(data = roi_simp1, fill = NA, colour = "black", linewidth = 0.20) +
  geom_sf(data = roi_simp0, fill = NA, colour = "black", linewidth = 0.40) +
  labs(
    title = "Nepal: Temperature and Precipitation Patterns", 
    subtitle = "Mean temperature and precipitation patterns based on 30 years of data.",
    caption = "Source: Terra Climate Data (1995-2025)     Author: Waruth POJSILAPACHAI"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 10)),
    plot.caption = element_text(size = 10, face = "bold", hjust = 0.5, 
                                margin = margin(t = 15)),
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10),
    plot.background = element_rect(fill = "white", color = NA)  # Ensure white background
  )

legend <- bi_legend(
  pal = pallet,   
  flip_axes = FALSE,
  rotate_pal = FALSE,
  dim = 4,
  xlab = "Temperature (\u00B0C)",
  ylab = "Precipitation (mm)",
  size = 8          # larger text (try 10–12)
)

finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) + # (0, 0) is the bottom-left and (1, 1) is the top-right corner
  #draw_plot(legend, 0.6, 0.2, 0.2, 0.2)   # smaller legend box
  draw_plot(legend, 0.7, 0.6, 0.2, 0.2)   # For Nepal

# 3D render: embed legend via annotation_custom instead of cowplot.
# annotation_custom places content inside the plot panel's data-coordinate
# space and does NOT change panel dimensions, so dem_map alignment is preserved.
# Adjust xmin/xmax/ymin/ymax (in degrees) to reposition the legend if needed.
map_with_legend <- map +
  annotation_custom(
    grob = ggplotGrob(legend),
    xmin = bbox$xmin,
    xmax = bbox$xmin + (bbox$xmax - bbox$xmin) * 0.28,
    ymin = bbox$ymax - (bbox$ymax - bbox$ymin) * 0.32,
    ymax = bbox$ymax
  )

finalPlot

out.dir <- "./output_png"

dir.create(out.dir)

ggsave(paste0(out.dir, "/", "NPL_Temp_PPT.png"), finalPlot, 
       width = 12, height = 11, dpi = 300, bg = "white")



# CREATE TERRAIN LAYER  ----------------------------------------------------------------------------

# Obtain AWS tiles DEM data from elevatr
# convert to terra SpatRaster and crop
dem <- elevatr::get_elev_raster(
  locations = roi_simp0, z = 8,
  clip = "locations"
) |> terra::rast() |>
  terra::crop(roi_simp0, mask = TRUE)

# project and convert DEM to dataframe — must use EPSG:4326 to match finalPlot
dem_df <- dem |>
  terra::project("EPSG:4326") |>
  as.data.frame(xy = TRUE, na.rm = TRUE)

# rename the third column to "dem"
names(dem_df)[3] <- "dem"

# Height map for 3D render — must be a structural mirror of `map` so that
# rayshader's pixel-by-pixel overlay aligns correctly.
# Key requirements:
#   1. Same theme_void(base_size=14) so title/subtitle/caption take identical space
#   2. Same xlim/ylim and labs/theme margins as `map`
#   3. A transparent geom_sf to trigger the same implicit coord_sf (CRS + expansion)
#      that `map` gets from its own geom_sf calls
dem_map <- ggplot() +
  theme_void(base_size = 14) +
  xlim(bbox$xmin, bbox$xmax) +
  ylim(bbox$ymin, bbox$ymax) +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = dem), show.legend = FALSE) +
  scale_fill_gradientn(colors = gray.colors(256, start = 0, end = 1)) +
  geom_sf(data = roi_simp0, fill = NA, colour = NA, linewidth = 0) +   # invisible; triggers coord_sf
  labs(
    title = "Nepal: Temperature and Precipitation Patterns",
    subtitle = "Mean temperature and precipitation patterns based on 30 years of data.",
    caption = "Source: Terra Climate Data (1995-2025)     Author: Waruth POJSILAPACHAI"
  ) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 10)),
    plot.caption  = element_text(size = 10, face = "bold", hjust = 0.5,
                                 margin = margin(t = 15)),
    plot.margin   = margin(t = 10, r = 10, b = 20, l = 10),
    plot.background = element_rect(fill = "white", color = NA)
  )


# Render 3D Scene -----------------------------------------------------------------------------
# Set white background
rgl::bg3d("white")

rayshader::plot_gg(
  ggobj = map_with_legend,          # use map, not finalPlot — cowplot overlay shifts the panel position
  ggobj_height = dem_map,
  width = 7,
  height = 7,
  windowsize = c(600, 600),
  scale = 100,
  shadow = TRUE,
  shadow_intensity = 1,
  phi = 87, theta = 0, zoom = .56,
  multicore = FALSE  # forking unreliable on macOS Apple Silicon
)

# zoom out
rayshader::render_camera(zoom = .6)


# Render 3D Objects ---------------------------------------------------------------------------
# Using explicit directional lights instead of an HDRI environment_light:
# - environment_light overrides backgroundhigh/backgroundlow, producing a dark background
# - preview = FALSE avoids opening an rgl window on macOS that can stall the render

rayshader::render_highquality(
  filename = "./output_png/NPL_Temp_PPT-3d.png",
  preview = FALSE,
  light = TRUE,
  lightdirection  = c(315, 45),   # key light NW, fill light NE
  lightaltitude   = c(60, 30),    # key light high, fill light low
  lightcolor      = c("white", "white"),
  lightintensity  = c(600, 150),
  parallel = FALSE,               # forking unreliable on macOS Apple Silicon
  width = 1500, height = 1500,
  interactive = FALSE,
  ground_material = rayrender::diffuse(color = "white"),
  backgroundhigh  = "white",
  backgroundlow   = "white"
)
