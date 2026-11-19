# ============================================================
# Georgia Sampling Locations + Ecoregions Map
# ============================================================

library(ggplot2)
library(sf)
library(maps)
library(ggrepel)

# ------------------------------------------------------------
# Load sampling coordinates
# ------------------------------------------------------------

maize_map <- read.csv(
  "C:/Aduragbemi/Manuscript/Review 5/map/coordinates.txt",
  header = TRUE,
  sep = ""
)

# ------------------------------------------------------------
# Load Georgia state boundary
# ------------------------------------------------------------

ga <- map_data("state", region = "georgia")

# ------------------------------------------------------------
# Load EPA ecoregions shapefile
# ------------------------------------------------------------

ecoreg <- st_read(
  "C:/Aduragbemi/Manuscript/Review 5/map/eco_files/ga_eco_l3.shp",
  quiet = TRUE
)

# Convert CRS to longitude/latitude
ecoreg <- st_transform(ecoreg, 4326)

# ------------------------------------------------------------
# Create label positions for ecoregions
# ------------------------------------------------------------

eco_labels <- st_centroid(ecoreg)

coords <- st_coordinates(eco_labels)

eco_labels$X <- coords[,1]
eco_labels$Y <- coords[,2]

# ------------------------------------------------------------
# Final map
# ------------------------------------------------------------

Final_map <- ggplot() +
  
  # Ecoregions
  geom_sf(
    data = ecoreg,
    aes(fill = US_L3NAME),
    color = "grey30",
    linewidth = 0.4,
    alpha = 0.35
  ) +
  
  # Georgia outline
  geom_polygon(
    data = ga,
    aes(long, lat, group = group),
    fill = NA,
    color = "black",
    linewidth = 1.2
  ) +
  
  # Sampling points
  geom_point(
    data = maize_map,
    aes(
      GPSlongitude,
      GPSlatitude
    ),
    shape = 21,
    fill = "dodgerblue3",
    color = "white",
    stroke = 0.6,
    size = 6
  ) +
  
  # Sampling location names
  geom_text_repel(
    data = maize_map,
    aes(
      GPSlongitude,
      GPSlatitude,
      label = Location
    ),
    size = 3.5,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  
  # Ecoregion labels
  geom_text_repel(
    data = eco_labels,
    aes(
      X,
      Y,
      label = US_L3NAME
    ),
    color = "black",
    fontface = "bold",
    size = 5.5,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  
  # Ecoregion colors
  scale_fill_manual(
    values = c(
      "Blue Ridge" = "#4daf4a",
      "Piedmont" = "#f4a582",
      "Southeastern Plains" = "#92c5de",
      "Southern Coastal Plain" = "#d9d9a3",
      "Ridge and Valley" = "#b39ddb",
      "Southwestern Appalachians" = "#80cdc1"
    )
  ) +
  
  coord_sf(
    xlim = c(-86, -80.5),
    ylim = c(30.5, 35.2),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Georgia Maize Sampling Locations Across Ecoregions"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    legend.position = "none",
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.2
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    )
  )

# ------------------------------------------------------------
# Display map
# ------------------------------------------------------------

Final_map

# ------------------------------------------------------------
# Save image
# ------------------------------------------------------------

ggsave(
  "Final_map.png",
  plot = Final_map,
  path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5",
  dpi = 700,
  width = 8,
  height = 7,
  units = "in"
)

ggsave(
  "Final_map.pdf",
  plot = Final_map,
  path = "C:/Aduragbemi/Manuscript/Review 5/PDF5",
  dpi = 700,
  width = 8,
  height = 7,
  units = "in"
)
