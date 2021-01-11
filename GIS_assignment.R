library(here)
library(spatstat)
library(here)
library(sp)
library(rgeos)
library(maptools)
library(GISTools)
library(tmap)
library(sf)
library(geojson)
library(geojsonio)
library(tmaptools)
library(rgdal)

#read the Zhejiang data in
ZhejiangTB <- st_read(here::here("Zhejiang", "Zhejiang.shp"))

#View the ZhejiangTB, add scale bar and compass   
tmap_mode("view") 
tm_shape(ZhejiangTB) +
  tm_polygons(col = NA, alpha = 0.5) +
  tm_scale_bar(position=c("left", "bottom"),text.size = 0.4) +
  tm_compass(type = "4star", position=c("left", "top"))

  
#Read the dataset of the TB incidence rate in Zhejiang
library(janitor)
TBData <- read.csv("ZJTB.csv")%>%
  clean_names()

#Merge the TBData and ZhejiangTB
ZhejiangTBMerged <- ZhejiangTB%>%
  merge(.,
        TBData,
        by.x="PAC",
        by.y="area_code")
tmap_mode("plot")

#View the distribution of TB incidences in 2010
ZhejiangTBMerged %>%
  qtm(.,fill = "x2010_incidence_rate")+
  tm_scale_bar(position=c("left", "bottom"),text.size = 0.4) + 
  tm_compass(type = "4star", position=c("left", "top"))

#View the distribution of TB incidences in 2019
ZhejiangTBMerged %>%
  qtm(.,fill = "x2019_incidence_rate")+
  tm_scale_bar(position=c("left", "bottom"),text.size = 0.4) + 
  tm_compass(type = "4star", position=c("left", "top"))

#First calculate the centroids of all Wards in Zhejiang
library(spdep)
coordsW <- ZhejiangTBMerged%>%
  st_centroid()%>%
  st_geometry()

plot(coordsW,axes=TRUE)

#create a neighbours list

Lward_nb <- ZhejiangTBMerged %>%
  poly2nb(., queen=T)

#plot them
plot(Lward_nb, st_geometry(coordsW), col="red")
#add a map underneath
plot(ZhejiangTBMerged$geometry, add=T)

#create a spatial weights object from these weights
Lward.lw <- Lward_nb %>%
  nb2listw(., style="C",zero.policy=TRUE)

head(Lward.lw$neighbours)

#############in 2010
#Global Moran's I statistics
I_LWard_Global_Density <- ZhejiangTBMerged %>%
  pull('x2010_incidence_rate') %>%
  as.vector()%>%
  moran.test(., Lward.lw,zero.policy=TRUE)
I_LWard_Global_Density

# Geary's C
C_LWard_Global_Density <- ZhejiangTBMerged %>%
  pull('x2010_incidence_rate') %>%
  as.vector()%>%
  geary.test(., Lward.lw,zero.policy=TRUE)
C_LWard_Global_Density

#General G
G_LWard_Global_Density <- ZhejiangTBMerged %>%
  pull('x2010_incidence_rate') %>%
  as.vector()%>%
  globalG.test(., Lward.lw,zero.policy=TRUE)
G_LWard_Global_Density

#############in 2019
#Global statistics in 2019
I_LWard_Global_Density <- ZhejiangTBMerged %>%
  pull('x2019_incidence_rate') %>%
  as.vector()%>%
  moran.test(., Lward.lw,zero.policy=TRUE)
I_LWard_Global_Density

C_LWard_Global_Density <- ZhejiangTBMerged %>%
  pull('x2019_incidence_rate') %>%
  as.vector()%>%
  geary.test(., Lward.lw,zero.policy=TRUE)
C_LWard_Global_Density

G_LWard_Global_Density <- ZhejiangWards %>%
  pull('x2019_incidence_rate') %>%
  as.vector()%>%
  globalG.test(., Lward.lw,zero.policy=TRUE)
G_LWard_Global_Density

###2010
#use the localmoran function to generate I for each ward in the city

I_LWard_Local_Density <- ZhejiangTBMerged %>%
  pull('x2010_incidence_rate') %>%
  as.vector()%>%
  localmoran(., Lward.lw)%>%
  as_tibble()

#what does the output (the localMoran object) look like?
slice_head(I_LWard_Local_Density, n=90)

I_LWard_Local_Density


ZhejiangTBMerged <- ZhejiangTBMerged %>%
  mutate(density_I =as.numeric(I_LWard_Local_Density$Ii))%>%
  mutate(density_Iz =as.numeric(I_LWard_Local_Density$Z.Ii))

#Mapping outputs
breaks1<-c(-1000,-2.58,-1.96,-1.65,1.65,1.96,2.58,1000)
MoranColours<- rev(brewer.pal(8, "RdGy"))

tm_shape(ZhejiangTBMerged) +
  tm_polygons("density_Iz",
              style="fixed",
              breaks=breaks1,
              palette=MoranColours,
              midpoint=NA,
              title="2010 TB indicate rate
(Local Moran's I)")+
  tm_scale_bar(position=c("left", "bottom"),text.size = 0.4) +
  tm_compass(type = "4star", position=c("left", "top"))

#Getis Ord  
Gi_LWard_Local_Density <- ZhejiangTBMerged %>%
  pull('x2010_incidence_rate') %>%
  as.vector()%>%
  localG(., Lward.lw,zero.policy=TRUE)

head(Gi_LWard_Local_Density)

Gi_LWard_Local_Density

ZhejiangTBMerged <- ZhejiangTBMerged %>%
  mutate(density_G = as.numeric(Gi_LWard_Local_Density))

GIColours<- rev(brewer.pal(8, "RdBu"))

#now plot on an interactive map
tm_shape(ZhejiangTBMerged) +
  tm_polygons("density_G",
              style="fixed",
              breaks=breaks1,
              palette=GIColours,
              midpoint=NA,
              title="2010 TB indicate rate 
              (Gi*)")+
  tm_scale_bar(position=c("left", "bottom"),text.size = 0.4) +
  tm_compass(type = "4star", position=c("left", "top"))

###2019
#use the localmoran function to generate I for each ward in the city
I_LWard_Local_Density <- ZhejiangTBMerged %>%
  pull('x2019_incidence_rate') %>%
  as.vector()%>%
  localmoran(., Lward.lw)%>%
  as_tibble()

#what does the output (the localMoran object) look like?
slice_head(I_LWard_Local_Density, n=90)

I_LWard_Local_Density


ZhejiangTBMerged <- ZhejiangTBMerged %>%
  mutate(density_I =as.numeric(I_LWard_Local_Density$Ii))%>%
  mutate(density_Iz =as.numeric(I_LWard_Local_Density$Z.Ii))

#Mapping outputs
breaks1<-c(-1000,-2.58,-1.96,-1.65,1.65,1.96,2.58,1000)
MoranColours<- rev(brewer.pal(8, "RdGy"))

tm_shape(ZhejiangTBMerged) +
  tm_polygons("density_Iz",
              style="fixed",
              breaks=breaks1,
              palette=MoranColours,
              midpoint=NA,
              title="2019 TB indicate rate
(Local Moran's I)")+
  tm_scale_bar(position=c("left", "bottom"),text.size = 0.4) +
  tm_compass(type = "4star", position=c("left", "top"))

#Getis Ord  
Gi_LWard_Local_Density <- ZhejiangTBMerged %>%
  pull('x2019_incidence_rate') %>%
  as.vector()%>%
  localG(., Lward.lw,zero.policy=TRUE)

head(Gi_LWard_Local_Density)

Gi_LWard_Local_Density

ZhejiangTBMerged <- ZhejiangTBMerged %>%
  mutate(density_G = as.numeric(Gi_LWard_Local_Density))

GIColours<- rev(brewer.pal(8, "RdBu"))

#now plot on an interactive map
tm_shape(ZhejiangTBMerged) +
  tm_polygons("density_G",
              style="fixed",
              breaks=breaks1,
              palette=GIColours,
              midpoint=NA,
              title="2019 TB indicate rate 
              (Gi*)")+
  tm_scale_bar(position=c("left", "bottom"),text.size = 0.4) +
  tm_compass(type = "4star", position=c("left", "top"))


