#Load Libraries
library(ggplot2)
library(ggmap)
library(maps)
library(mapproj)
library(mapdata)
library(maptools)
library(sp)
library(raster)
library(rgdal)
library(dismo)

#register the API key
register_stadiamaps(key = "5c3883dc-5df8-40a6-bd7a-8e052bfe68bc") 

#input the coordinates of each location
maize_map = read.csv("C:/Aduragbemi/Manuscript/Review 5/map/coordinates.txt", h=T, sep="")

#calculate the range of the latitude and longitude
range(maize_map$GPSlatitude)
range(maize_map$GPSlongitude)

#use the calculated range to set for the coordinate range size for the map
maize_base = get_map(location=c(-87,30,-80,36), source="stadia", zoom=7, maptype="stamen_terrain")
maize_map1 = ggmap(maize_base)
maize_map1


#make final plot using ggplot
Final_map <- maize_map1 +
  geom_point(
    data = maize_map,
    aes(x = GPSlongitude, y = GPSlatitude, fill = Location),
    shape = 21,
    color = "white",
    size = 7
  ) +
  scale_fill_manual(
    values = c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
               "#e6ab02","#a6761d","#666666","#1f78b4","#b2df8a",
               "#fb9a99","#cab2d6","#fdbf6f","#3c3c3c","#6a3d9a","#ff7f00"),
    labels = c("Arlington","Athens","Blairsville","Cave Spring","Dawson",
               "Fitzgerald","Fort Valley","Hawkinsville","Midville","Pavo",
               "Plains","Rome","Tennile","Tifton","Wadley","Watkinsville"),
    name = NULL
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Sampling sites"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = rel(0.75)),
    legend.key = element_rect(colour = "white"),
    axis.text.x = element_text(angle = 45, vjust = 0.5)
  )

Final_map

#save as image
ggsave("Final_map.png", plot = Final_map , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")


#save as pdf
ggsave("Final_map.pdf", plot = Final_map, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


