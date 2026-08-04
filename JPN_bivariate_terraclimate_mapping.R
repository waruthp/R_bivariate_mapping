# =============================================================================
# Bivariate mapping (Japan) with hillshade background (2D)
# =============================================================================

# Install and load packages ---------------------------------------------------
pacman::p_load(
  tidyverse, sf, terra, geobounds, climateR, geodata, 
  biscale, ggplot2, dplyr, cowplot, countrycode, elevatr, ggnewscale
)

# Set cache directory for geoBoundaries ---------------------------------------
newdir <- file.path(tempdir(), "/geoboundvignette")
gb_set_cache_dir(newdir)

countrycode(
  c("Thailand", "Japan", "Pakistan", "Srilanka", "Myanmar", "Lao", "Vietnam"),
  origin = "country.name",
  destination = "iso3c"
)

# Select Japan boundary data ----------------------------------------------------
iso3c <- "JPN"

roi_simp0 <- gb_get_adm0(country = iso3c, simplified = TRUE) |> 
  mutate(res = "Simplified")
roi_simp1 <- gb_get_adm1(country = iso3c, simplified = TRUE) |> 
  mutate(res = "Simplified")

# Project to Japan Albers Equal Area (metres)
crs_proj <- "EPSG:6690"   # JGD2011 / Japan Albers Equal Area
roi_simp0_proj <- st_transform(roi_simp0, crs_proj)
roi_simp1_proj <- st_transform(roi_simp1, crs_proj)

bbox <- st_bbox(roi_simp1_proj)

# Temperature data ------------------------------------------------------------
temp_terra <- getTerraClim(AOI = roi_simp0, 
                           varname = c("tmax", "tmin"),
                           startDate = "1995-01-01",
                           endDate   = "2025-12-31")

mean_temp_monthly <- (temp_terra$tmax + temp_terra$tmin) / 2
annual_mean_temp <- mean(mean_temp_monthly, na.rm = TRUE)

# Project and resample
annual_mean_temp <- terra::project(annual_mean_temp, crs_proj, method = "bilinear")
new_res <- 1000   # metres (adjust for desired resolution)
new_raster <- rast(ext(annual_mean_temp), resolution = new_res, crs = crs_proj)
annual_mean_temp <- resample(annual_mean_temp, new_raster, method = "bilinear")
roi_temp <- terra::crop(annual_mean_temp, roi_simp0_proj, mask = TRUE)

# Precipitation data ----------------------------------------------------------
ppt_terra <- getTerraClim(AOI = roi_simp0, 
                          varname = "ppt",
                          startDate = "1995-01-01",
                          endDate   = "2025-12-31")

annual_mean_ppt <- mean(ppt_terra$ppt, na.rm = TRUE)
annual_mean_ppt <- terra::project(annual_mean_ppt, crs_proj, method = "bilinear")
annual_mean_ppt <- resample(annual_mean_ppt, new_raster, method = "bilinear")
roi_ppt <- terra::crop(annual_mean_ppt, roi_simp0_proj, mask = TRUE)

# Stack rasters ---------------------------------------------------------------
temp_ppt <- c(roi_temp, roi_ppt)
names(temp_ppt) <- c("temp", "ppt")

# Convert to data frame for bivariate classification
temp_ppt_df <- temp_ppt |> 
  as.data.frame(xy = TRUE, na.rm = TRUE)

# Bivariate classification ----------------------------------------------------
data <- bi_class(temp_ppt_df,
                 x = temp, 
                 y = ppt, 
                 style = "quantile", dim = 4)

# =============================================================================
# HILLSHADE from DEM
# =============================================================================
message("Downloading and processing DEM for hillshade...")
dem <- elevatr::get_elev_raster(
  locations = roi_simp0, z = 9,        # resolution ~30 m
  clip = "locations"
) |> terra::rast() |>
  terra::crop(roi_simp0, mask = TRUE)

# Project and align with the bivariate grid
dem <- terra::project(dem, crs_proj, method = "bilinear")
dem <- terra::resample(dem, temp_ppt, method = "bilinear")   # same grid

# Compute slope and aspect, then hillshade
slope <- terra::terrain(dem, "slope", unit = "radians")
aspect <- terra::terrain(dem, "aspect", unit = "radians")
hill <- terra::shade(slope, aspect, angle = 45, direction = 315)   # sun angle and direction

# Convert hillshade to data frame for plotting
hill_df <- hill |> as.data.frame(xy = TRUE, na.rm = TRUE)
names(hill_df)[3] <- "hill"

# =============================================================================
# 2D MAPPING with hillshade background
# =============================================================================
pallet <- "BlueOr"

map <- ggplot() +
  # Hillshade layer (grayscale)
  geom_tile(data = hill_df, aes(x = x, y = y, fill = hill), show.legend = FALSE) +
  scale_fill_gradient(low = "black", high = "white", na.value = NA) +
  new_scale_fill() +   # to allow a second fill scale for bivariate tiles
  # Bivariate tiles on top with transparency
  geom_tile(data = data, 
            aes(x = x, y = y, fill = bi_class),
            alpha = 0.6, show.legend = FALSE) +
  bi_scale_fill(pal = pallet, dim = 4, flip_axes = FALSE, rotate_pal = FALSE) +
  # Administrative boundaries
  geom_sf(data = roi_simp1_proj, fill = NA, colour = "black", linewidth = 0.20) +
  geom_sf(data = roi_simp0_proj, fill = NA, colour = "black", linewidth = 0.40) +
  # Limits and theme
  xlim(bbox$xmin, bbox$xmax) +
  ylim(bbox$ymin, bbox$ymax) +
  theme_void(base_size = 14) +
  labs(
    title = "Japan: Temperature and Precipitation Patterns", 
    subtitle = "Mean temperature and precipitation patterns based on 30 years of data.\nHillshade from 30m DEM",
    caption = "Source: Terra Climate Data (1995-2025) | DEM: elevatr | Author: Waruth POJSILAPACHAI"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 10)),
    plot.caption = element_text(size = 10, face = "bold", hjust = 0.5, 
                                margin = margin(t = 15)),
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10)
  )

# Create separate legend using bi_legend
legend <- bi_legend(
  pal = pallet,   
  flip_axes = FALSE,
  rotate_pal = FALSE,
  dim = 4,
  xlab = "Temperature (\u00B0C)",
  ylab = "Precipitation (mm)",
  size = 8
)

# Combine map and legend
finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.75, 0.2, 0.2, 0.2)

# Display the plot
print(finalPlot)

# Save output ----------------------------------------------------------------
out.dir <- "./output_png"
dir.create(out.dir, showWarnings = FALSE)

ggsave(paste0(out.dir, "/", "JPN_Temp_PPT_Hillshade.png"), finalPlot, 
       width = 12, height = 11, dpi = 300, bg = "white")

message("Hillshade map saved to ./output_png/JPN_Temp_PPT_Hillshade.png")
