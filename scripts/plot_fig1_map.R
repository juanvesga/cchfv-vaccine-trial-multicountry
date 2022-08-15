
# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script creates maps for the three countries of study, marking specific provinces

rm(list = ls())
library(maptools)
library(raster)
library(here)
library(rgdal)
library(ggplot2)
library(gridExtra)

shp <- readOGR(
  dsn = here("spatial", "ne_10m_admin_0_countries"),
  layer = "ne_10m_admin_0_countries"
)



# Colors by country
afcol <- "gold"
sacol <- "olivedrab3"
tkcol <- "tomato2"

# AFG map
shp.afg <- readOGR(dsn = here("spatial", "AFG"), layer = "gadm40_AFG_1")
shp.afg2 <- shp.afg[shp.afg$NAME_1 == "Hirat", ]

dfw <- fortify(shp)
df <- fortify(shp.afg)
df2 <- fortify(shp.afg2)

AFG <- ggplot() +
  geom_polygon(
    data = dfw,
    aes(x = long, y = lat, group = group),
    fill = "grey98",
    color = "grey88"
  ) +
  geom_polygon(
    data = df,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "black"
  ) +
  geom_polygon(
    data = df2,
    aes(x = long, y = lat, group = group),
    fill = afcol,
    color = "black"
  ) +
  theme_bw() +
  labs(title = "Afghanistan") +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  coord_cartesian(xlim = c(55, 80), ylim = c(25, 45))

# RSA map
shp.sa <- readOGR(
  dsn = here("spatial", "SA"),
  layer = "gadm40_ZAF_1"
)
shp.sa2 <- shp.sa[
  shp.sa$NAME_1 == "Northern Cape" |
    shp.sa$NAME_1 == "Free State" |
    shp.sa$NAME_1 == "North West",
]


df <- fortify(shp.sa)
df2 <- fortify(shp.sa2)

RSA <- ggplot() +
  geom_polygon(
    data = dfw,
    aes(x = long, y = lat, group = group),
    fill = "grey98",
    color = "grey88"
  ) +
  geom_polygon(
    data = df,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "black"
  ) +
  geom_polygon(data = df2, aes(x = long, y = lat, group = group), fill = sacol, color = "black") +
  theme_bw() +
  labs(title = "South Africa") +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  coord_cartesian(xlim = c(10, 40), ylim = c(-40, -15))


# Turkey map

shp.tk <- readOGR(dsn = here("spatial", "TKY"), layer = "gadm40_TUR_1")
shp.tk2 <- shp.tk[shp.tk$NAME_1 == "Tokat" |
  shp.tk$NAME_1 == "Sivas" |
  shp.tk$NAME_1 == "Erzincan" |
  shp.tk$NAME_1 == "Erzurum" |
  shp.tk$NAME_1 == "GÃ¼mÃ¼shane", ]

df <- fortify(shp.tk)
df2 <- fortify(shp.tk2)

png(here("output", "file.png"), width = 800, height = 800)

TKY <- ggplot() +
  geom_polygon(
    data = dfw,
    aes(x = long, y = lat, group = group),
    fill = "grey98",
    color = "grey88"
  ) +
  geom_polygon(
    data = df,
    aes(x = long, y = lat, group = group),
    fill = "grey96",
    color = "black"
  ) +
  geom_polygon(
    data = df2,
    aes(x = long, y = lat, group = group),
    fill = tkcol,
    color = "black"
  ) +
  theme_bw() +
  labs(title = "Turkey") +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  coord_cartesian(xlim = c(21, 49), ylim = c(31, 49))

TKY
dev.off()

# World map with country outlines
df <- fortify(shp)
df.af <- fortify(shp.afg)
df.sa <- fortify(shp.sa)
df.tk <- fortify(shp.tk)

wrld <- ggplot() +
  geom_polygon(
    data = df,
    aes(x = long, y = lat, group = group),
    fill = "grey85",
    color = "grey85"
  ) +
  geom_polygon(
    data = df.af,
    aes(x = long, y = lat, group = group),
    fill = afcol,
    color = afcol
  ) +
  geom_polygon(
    data = df.sa,
    aes(x = long, y = lat, group = group),
    fill = sacol,
    color = sacol
  ) +
  geom_polygon(
    data = df.tk,
    aes(x = long, y = lat, group = group),
    fill = tkcol, color = tkcol
  ) +
  theme_bw() +
  labs(title = "Study locations") +
  theme(
    # Hide panel borders and remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  coord_cartesian(ylim = c(-100, 100))

## Get everything together in one fig
windows()
g<-gridExtra::grid.arrange(wrld, AFG, RSA, TKY, ncol = 2)

# Save figure 1
ggsave(path = here("output"), filename = "Fig1.tiff", device = "tiff",plot=g)

#######
## End of code
#######
