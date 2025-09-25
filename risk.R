################################################################################
####      ECOLOGICAL RISK-BASED APPROACH TO FACILITATE LICENSING OFFSHORE   ####
####                          WIND DEVELOPMENT                              ####
####                                                                        ####
####                          PART 3 (of 3): RISK                           ####
################################################################################

# This script relates to work in Bolam, S.G., Cooper, K.M., and 
# Downie, A-L. Developing an Ecological Risk-Based Approach to Facilitate
# Licensing Offshore Wind Development. Ecosphere.

# In this study we create a Combined Risk layer based on 3 element layers for: 
# Biodiversity, Sensitivity and Assemblage Rarity.

# This script 'PART 2: RISK is used to develop the final Combined Risk layer,
# together with an associated confidence layer. 

# Input data are raster files for Sensitivity (see PART 1: SENSITIVITY;
# sensitivity.R), Biodiversity (see https://github.com/keithmcooper/OneBenthicBiodiversity),
# and Assemblage Rarity (see PART 2: ASSEMBLAGES; assemblages.R).

# Note than individual maps are formatted for use in paper figures (2, 4, 5, 6). 
# As such, the ggplot code made need tweaking if other maps are plotted 
# individually or in combination with others.
#_______________________________________________________________________________
#### RISK ELEMENT MAPS ####

## Load packages
library(sf)
library(terra)
library(tidyterra)
library(rasterVis)
library(ggplot2)
library(colorRamps)
library(raster)
library(dplyr)
library(RColorBrewer)

## Load countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))
#_______________________________________________________________________________
#### SENSITIVITY: MODEL MAP ####

## Load sensitivity raster (spatraster)
sensitivity_model_rast <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Sensitivity/Model/sensitivity_Mean_SQ.tif')
#plot(sensitivity_model )

## Reduce size of raster
#sensitivity_model_agg <- aggregate(sensitivity_model_rast, fact=7,fun = modal)
#sensitivity_model_rast <- sensitivity_model_agg

## Sensitivity map
sensitivity_model <- ggplot() +
  geom_spatraster(data = sensitivity_model_rast) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 20)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 20))+#make legend background transparent and text white
  theme(legend.position =c(0.9,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.key.size = unit(0.6, "cm"))+#t, r, b, l
  theme(
    panel.background = element_blank(),
    axis.title.y = element_text(colour = "white"),
    axis.title.x = element_blank(),
    axis.text.x=element_text(colour = "white"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

sensitivity_model

## Save sensitivity map.
ggsave(plot = sensitivity_model,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\sensitivity_model.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
#### SENSITIVITY: CONFIDENCE MAP ####

## Load sensitivity raster (spatraster)
sensitivity_conf_rast <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Sensitivity/Confidence/sensitivity_CV_SQ.tif')
#plot(sensitivity_model )

## Reduce size of raster
#sensitivity_conf_rast_agg <- aggregate(sensitivity_conf_rast, fact=7,fun = modal)
#sensitivity_conf_rast <- sensitivity_conf_rast_agg

## Sensitivity confidence map
sensitivity_conf <- ggplot() +
  geom_spatraster(data = sensitivity_conf_rast) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name="CV")+#OrangesGreys
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw(base_size = 18)+
  annotate("text",x=c(-7.8),y=c(60.05),label=c("Sensitivity"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12),legend.title = element_text(color = "white"))+#make legend background transparent and text white #,legend.title = element_blank()
  theme(plot.margin = unit(c(0,0.2,0,1), "cm"),legend.key.size = unit(0.45, "cm"))+#t, r, b, l
  theme(legend.position=c(0.87,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(
    axis.title.x = element_text(colour = "white"),
    axis.title.y = element_text(colour = "white"),
    axis.text.x = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "white", size = 12))

sensitivity_conf

## Save sensitivity map
ggsave(plot = sensitivity_conf,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\sensitivity_confidence.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,

#_______________________________________________________________________________
### SENSITIVITY: MODEL AND CONFIDENCE MAPS ####

## Load package
library(ggpubr)

## Stitch plots together
sensitivity_model_conf <- ggpubr::ggarrange(sensitivity_model,sensitivity_conf,nrow=1,labels = c("a)", "b)"),font.label=list(color="black",size=21,face='plain'),align="v")#font.label=list(color="black",size=6,face='plain'),align="v") ,labels = c("D0", "D1", "D2","N")

## Add x and y labels
sensitivity_model_conf2 <- annotate_figure(
  sensitivity_model_conf, 
  bottom = text_grob("          Longitude", 
                     color = "black", face = "plain", size = 21),#,
  left = text_grob("Latitude", 
                   color = "black", face = "plain", size = 21,rot = 90)
)

## Save combined plot
ggsave(plot = sensitivity_model_conf2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\sensitivity_model_conf.png"),
       height = 240, width =480, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285

#_______________________________________________________________________________
#### BIODIVERSITY: MODEL MAP ####

## Load biodiversity spatraster
biodiv_model_rast <- rast("C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Bio (Combined Biodiversity)/Model/BiodiversityClusterMaxClass_Aug25.tif")

## Make cluster a factor
values(biodiv_model_rast) <- as.factor(values(biodiv_model_rast))
unique(biodiv_model_rast)

## Assemblage map
biodiv_model <- ggplot() +
  geom_spatraster(data = biodiv_model_rast) +
  scale_fill_manual(breaks = c('2','8','5','1','4','6','7','3'), values = c(
    "#37C331", "#6FD326" ,"#A8E21B", "#E0F210", "#F3F223","#E2E256",  "#D1D189",  "#C0C1BC", "white"),
    labels = c('Bio-A','Bio-B','Bio-C','Bio-D','Bio-E','Bio-F','Bio-G','Bio-H'),
    na.translate = F)+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 20)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 20))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.key.size = unit(0.6, "cm"))+#t, r, b, l
  theme(
    panel.background = element_blank(),
    axis.title.y = element_text(colour = "white"),
    axis.title.x = element_blank(),
    axis.text.x=element_text(colour = "white")
  )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

biodiv_model

## Save biodiversity map
ggsave(plot = biodiv_model,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\biodiv_model.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
#### BIODIVERSITY: CONFIDENCE MAP ####

## Load biodiversity spatraster
biodiv_conf_rast <-  rast("C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Bio (Combined Biodiversity)/Confidence/BiodiversityClusterConfidence_Aug25.tif")

## Assemblage map
biodiv_conf <- ggplot() +
  geom_spatraster(data = biodiv_conf_rast) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name="Confidence")+#Oranges
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw(base_size = 18)+
  annotate("text",x=c(-7.5),y=c(60.05),label=c("Biodiversity"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12),legend.title = element_text(color = "white"))+#make legend background transparent and text white #,legend.title = element_blank()
  theme(plot.margin = unit(c(0,0.2,0,1), "cm"),legend.key.size = unit(0.45, "cm"))+#t, r, b, l
  theme(legend.position=c(0.87,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(
    axis.title.x = element_text(colour = "white"),
    axis.title.y = element_text(colour = "white"),
    axis.text.x = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "white", size = 12))

biodiv_conf

## Save biodiversity map
ggsave(plot = biodiv_conf,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\biodiv_conf.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
### BIODIVERSITY: MODEL AND CONFIDENCE MAPS ####

## Load package
library(ggpubr)

## Stitch plots together
biodiv_model_conf <- ggpubr::ggarrange(biodiv_model,biodiv_conf,nrow=1,labels = c("a)", "b)"),font.label=list(color="black",size=21,face='plain'),align="v")#font.label=list(color="black",size=6,face='plain'),align="v") ,labels = c("D0", "D1", "D2","N")

## Add x and y labels
biodiv_model_conf2 <- annotate_figure(
  biodiv_model_conf, 
  bottom = text_grob("          Longitude", 
                     color = "black", face = "plain", size = 21),#,
  left = text_grob("Latitude", 
                   color = "black", face = "plain", size = 21,rot = 90))

biodiv_model_conf2

## Save combined plot
ggsave(plot = biodiv_model_conf2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\biodiv_model_conf.png"),
       height = 240, width =480, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285
#_______________________________________________________________________________
#### ASSEMBLAGES: MODEL MAP ####

## Load assemblage raster
assemblage_model_rast<- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Assemblage/Model/AssemblageMaxClass_Nov24.tif')

## Add attibute cluster group labels
clus <- data.frame(id=1:11, cluster=c('A1','A2a','A2b','B1b','C1a','C1b','D1','D2a','D2b','D2c','D2d'))
levels(assemblage_model_rast) <- clus
is.factor(assemblage_model_rast)
assemblage_model_rast

## Assemblage map
assemblage_model <- ggplot() +
  geom_spatraster(data = assemblage_model_rast) +
  scale_fill_manual(values = c(  "#0000ee","#00ffff","#05aac1","#9a32cd","#00cd00","#9aff9a","#b40202","#ff0000","#ff8c00","#ffff00","#b4b404",'white'),name="cluster", labels=c("A1","A2a","A2b","B1b","C1a","C1b","D1","D2a","D2b","D2c","D2d"),na.value = "transparent", na.translate = FALSE)+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 20)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 20))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.22))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.key.size = unit(0.6, "cm"))+#t, r, b, l
  theme(
    axis.title.y = element_text(colour = "white"),
    axis.title.x = element_blank())+
  theme(legend.title = element_text( size=2), legend.text=element_text(size=18))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

assemblage_model

## Save assemblage map
ggsave(plot = assemblage_model,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\assemblage_model.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
#### ASSEMBLAGES: CONFIDENCE MAP ####

## Load assemblage raster
assemblage_conf_rast <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Assemblage/Confidence/AssemblageConfidence_Nov24.tif')

## Assemblage map
assemblage_conf <- ggplot() +
  geom_spatraster(data = assemblage_conf_rast) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name="Confidence")+#Oranges
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw(base_size = 18)+
  annotate("text",x=c(-9),y=c(60.05),label=c("Rarity"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12),legend.title = element_text(color = "white"))+#make legend background transparent and text white #legend.title = element_blank()
  theme(plot.margin = unit(c(0,0.2,0,1), "cm"),legend.key.size = unit(0.45, "cm"))+#t, r, b, l
  theme(legend.position=c(0.87,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(colour = "white"),
        axis.title.x = element_text(colour = "white"),
        legend.title = element_text(color = "white", size = 12))


assemblage_conf

## Save assemblage map
ggsave(plot = assemblage_conf,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\assemblage_conf.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
### ASSEMBLAGES: MODEL AND CONFIDENCE MAPS ####

## Load package
library(ggpubr)

## Stitch plots together
assemblage_model_conf <- ggpubr::ggarrange(assemblage_model,assemblage_conf,nrow=1,labels = c("a)", "b)"),font.label=list(color="black",size=21,face='plain'),align="v")#font.label=list(color="black",size=6,face='plain'),align="v") ,labels = c("D0", "D1", "D2","N")

## Add x and y labels
assemblage_model_conf2 <- annotate_figure(
  assemblage_model_conf, 
  bottom = text_grob("          Longitude", 
                     color = "black", face = "plain", size = 21),#,
  left = text_grob("Latitude", 
                   color = "black", face = "plain", size = 21,rot = 90))

assemblage_model_conf2

## Save combined plot
ggsave(plot = assemblage_model_conf2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\assemblage_model_conf.png"),
       height = 240, width =480, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285

#_______________________________________________________________________________
#### SENSITIVITY, BIODIVERSITY, ASSEMBLAGES: MODEL AND CONFIDENCE PLOTS (3 x 2 PLOT) ####
library(egg)
library(ggpubr)

## Elements:
s1 <- sensitivity_model
s2 <- sensitivity_conf

b1 <- biodiv_model
b2 <- biodiv_conf

a1 <- assemblage_model#a1
a2 <- assemblage_conf#a2

# Stitch plots together
s_stitch <- egg::ggarrange(s1,s2, labels = c("", ""),nrow=1)
b_stitch <- egg::ggarrange(b1,b2, labels = c("", ""),nrow=1)
a_stitch <- egg::ggarrange(a1,a2, labels = c("", ""),nrow=1)

sens_biodiv_assem_models_conf <- ggpubr::ggarrange(s_stitch,b_stitch,a_stitch, labels = c("a)", "b)", "c)"),nrow=3,font.label=list(color="black",size=24,face='plain'),align="h",widths = c(0.81,0.81,1))

sens_biodiv_assem_models_conf2 <- annotate_figure(
  sens_biodiv_assem_models_conf, 
  bottom = text_grob("      Longitude", color = "black", face = "plain", size = 24),
  left = text_grob("Latitude", color = "black", face = "plain", size = 24, rot = 90))

sens_biodiv_assem_models_conf2

## Save
ggsave(plot = sens_biodiv_assem_models_conf2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\sens_biodiv_assem_models_conf.png"),
       height = 600, width =450, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285

#_______________________________________________________________________________
#### MANIPULATION OF RASTER ELEMENT LAYERS ####
#_______________________________________________________________________________
#### RESCALE FUNCTION ####

#https://gis.stackexchange.com/questions/437520/normalize-raster-in-r
# Function to be applied on the raster values; return: SpatRaster object
rescale01 <- function(x) {
  val <- values(x)
  values(x) <- (val - min(val, na.rm = TRUE)) / (max(val, na.rm = TRUE) - min(val, na.rm = TRUE))
  x
}
#_______________________________________________________________________________
#### LOAD RISK ELEMENT MODELS ####

## Load package
library(terra)

## Load rasters
sensitivity_model_rast <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Sensitivity/Model/sensitivity_Mean_SQ.tif')
biodiv_model_rast <- rast("C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Bio (Combined Biodiversity)/Model/BiodiversityClusterMaxClass_Aug25.tif")
assemblage_model_rast<- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Assemblage/Model/AssemblageMaxClass_Nov24.tif')
#_______________________________________________________________________________
#### SENSITIVITY: RESCALE  ####

## Rescale values of SpatRaster object with nlyr = 1
sensitivity_model_scaled <- rescale01(sensitivity_model_rast)
plot(sensitivity_model_scaled)

# Save as .tiff
writeRaster(sensitivity_model_scaled,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\sensitivity_Mean_SQ_scaled.tif',overwrite=TRUE)#,format = "GTiff"
#_______________________________________________________________________________
#### ASSEMBLAGES: CALCULATE AREA (BY CLASS, TABLE 3) ####

## Load packages
library(DT)
library(dplyr)
library(knitr)
library(kableExtra)
library(scales)
library(kableExtra)
library(knitr)
library(kableExtra)
library(magick)
library(webshot2)
library(magick)

## Identify unique values
unique_values <- unique(values(assemblage_model_rast, na.rm=TRUE))#
unique_values

## Initialize a data frame to store results
assem_area_df <- data.frame(value = unique_values, area = NA)

## Area by assemblage class
for (i in 1:length(unique_values)) {
  value <- unique_values[i]
  mask <- assemblage_model_rast == value
  assem_area_df$area[i] <- sum(cellSize(mask, unit="km")[values(mask, na.rm=TRUE)])
}

# Print the results
print(assem_area_df)

## Calculate the 'total area' and 'percentage'
assem_area_df %>%
  select_if(is.numeric) %>%
  summarize(across(everything(), sum))

## Add col for total area
assem_area_df$tot_area <- 903466.3
assem_area_df

## Create and format
assemblage_area_df2 <-assem_area_df %>%
  
  mutate(
    AssemblageMaxClass_Nov24 = case_when(
      AssemblageMaxClass_Nov24 == 1 ~ 'A1', 
      AssemblageMaxClass_Nov24 == 2 ~ 'A2a',
      AssemblageMaxClass_Nov24 == 3 ~ 'A2b',
      AssemblageMaxClass_Nov24 == 4 ~ 'B1b',
      AssemblageMaxClass_Nov24 == 5 ~ 'C1a',
      AssemblageMaxClass_Nov24 == 6 ~ 'C1b',
      AssemblageMaxClass_Nov24 == 7 ~ 'D1',
      AssemblageMaxClass_Nov24 == 8 ~ 'D2a',
      AssemblageMaxClass_Nov24 == 9 ~ 'D2b',
      AssemblageMaxClass_Nov24 == 10 ~ 'D2c',
      AssemblageMaxClass_Nov24 == 11 ~ 'D2d'),
    score = tot_area - area,
    `Class % of Total` = (area / tot_area) * 100,
    `Rarity %` = (score / tot_area) * 100) %>%
  arrange(desc(score)) %>% # order by rarity (area not occupied)
  mutate(area = round(area,0)) %>%
  mutate(tot_area = round(tot_area,0)) %>%
  mutate(score = round(score,0)) %>%
  mutate(`Class % of Total` = round(`Class % of Total`,2)) %>%
  mutate(`Rarity %` = round(`Rarity %`,2),
         area = comma(area),
         tot_area = comma(tot_area),
         score = comma(score))%>%
  rename(
    'Group' = 'AssemblageMaxClass_Nov24',
    'Area (km²)' = 'area',
    'Total Area (km²)' = 'tot_area',
    'Area Not Occupied (km²)' = 'score')

assemblage_area_df2 

## Update column order
assemblage_area_df3 <- assemblage_area_df2[,c(1,2,3,5,4,6)]
assemblage_area_df3

## Add colouring to Group column in kable output (optional)
assemblage_area_df3$Group = cell_spec(
  assemblage_area_df3$Group, #color = "black", align = "c", angle = 45, 
  background = factor(assemblage_area_df2$Group, c("A1","A2a","A2b","B1b","C1a","C1b","D1","D2a","D2b","D2c","D2d"), c(
    "#0000ee","#00ffff","#05aac1","#9a32cd","#00cd00","#9aff9a","#b40202","#ff0000","#ff8c00","#ffff00","#b4b404")))

# Define the path to save the HTML file
html_file <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_faunal_assemblage_area_scores.html"

## Save table
kable(assemblage_area_df3, escape=FALSE,align = c("l", "c","c"),format = "html",caption = "") %>%#Table 3. Spatial extent scores (Inverse Area km2) for faunal assemblage clusters.
  column_spec (1,color = "white",width = "8em", bold = F)%>%
  column_spec(2, width = "15em")%>%  # Adjust the width as needed 
  column_spec(3, width = "15em")%>%  # Adjust the width as needed 
  column_spec(4, width = "15em")%>%  # Adjust the width as needed 
  column_spec(5, width = "20em")%>%  # Adjust the width as needed 
  column_spec(6, width = "15em")%>%  # Adjust the width as needed 
  kable_styling(full_width = FALSE)%>%#"striped","hover",,bootstrap_options = c( "condensed")
  #add_header_above(c("Table 3. Spatial extent scores for faunal assemblage clusters." =3))%>%
  # row_spec(nrow(assemblage_area_df2), color = "black",bold = TRUE)%>%
  save_kable(file = html_file)

# Capture the screenshot of the HTML file
png_file <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_faunal_assemblage_area_scores.png"
webshot(html_file, file = png_file)

# Read the PNG file
img <- image_read("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_faunal_assemblage_area_scores.png")

# Trim the white space
img_trimmed <- image_trim(img)

# Save the trimmed image
image_write(img_trimmed, "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_3.png")
#_______________________________________________________________________________
#### ASSEMBLAGES: SUBSTITUTE (RASTER CLUSTERS WITH AREA NOT OCCUPIED BY CLASS) ####

## Substitute cluster classes for Area Score (see Table 3)
assemblage_model_rast1 <- subst(assemblage_model_rast, from = 1, to = 901893)# A1
assemblage_model_rast2 <- subst(assemblage_model_rast1, from = 2, to = 902574)# A2a
assemblage_model_rast3 <- subst(assemblage_model_rast2, from = 3, to = 892190)# A2b
assemblage_model_rast4 <- subst(assemblage_model_rast3, from = 4, to = 886010)# B1b
assemblage_model_rast5 <- subst(assemblage_model_rast4, from = 5, to = 876874)# C1a
assemblage_model_rast6 <- subst(assemblage_model_rast5, from = 6, to = 893828)# C1b
assemblage_model_rast7 <- subst(assemblage_model_rast6, from = 7, to = 816637)# D1
assemblage_model_rast8 <- subst(assemblage_model_rast7, from = 8, to = 741526)# D2a
assemblage_model_rast9 <- subst(assemblage_model_rast8, from = 9, to = 579850)# D2b
assemblage_model_rast10 <- subst(assemblage_model_rast9, from = 10, to = 754748)# D2c
assemblage_model_rast11 <- subst(assemblage_model_rast10, from = 11, to = 788533)# D2d
plot(assemblage_model_rast11)

# Save as .tiff
writeRaster(assemblage_model_rast11,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\assemblage_rarity.tif',overwrite=TRUE)#,format =
#_______________________________________________________________________________
#### ASSEMBLAGE RARITY: MAP (AREA NOT OCCUPIED BY CLASS) ####

## Load package
library(viridis)
library(ggplot2)
library(terra)
library(tidyterra)

## Plot of assemblagesbased on 'Area Not Occupied (km²)'
rarity_plot <- ggplot() +
  geom_spatraster(data = assemblage_model_rast11) +
  scale_fill_viridis(na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 22)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 22))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.20))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_text(colour = "white"))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.key.size = unit(0.8, "cm"))+#t, r, b, l
  theme(
    axis.title.y = element_text(colour = "white"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(fill = expression(km^2))

rarity_plot

## Save assemblage map
ggsave(plot = rarity_plot,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\rarity_plot.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
#### ASSEMBLAGE RARITY: RESCALE (AREA NOT OCCUPIED BY CLASS) ####

## Rescale values of SpatRaster object 
assemblage_model_rast_rarity_scaled <- rescale01(assemblage_model_rast11)
plot(assemblage_model_rast_rarity_scaled)

## Save as .tiff
writeRaster(assemblage_model_rast_rarity_scaled,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\assemblage_rarity_scaled.tif',overwrite=TRUE)#,format = "GTiff"
#_______________________________________________________________________________
#### BIODIVERSITY: SUBSTITUTE (CLASS TO CLUSTER CENTRE SUM) ####

## Identify unique values
unique_values <- unique(values(biodiv_model_rast, na.rm=TRUE))
unique_values

## Dataframe
biodiv_centre_sums <- data.frame(unique_values)
biodiv_centre_sums

## Change nameof 1st column
colnames(biodiv_centre_sums)[1] <- 'Cluster'

## Make Cluster a factor
biodiv_centre_sums$Cluster <- factor(biodiv_centre_sums$Cluster,levels=c(2,8,5,1,4,6,7,3))

## Order by cluster (highest to lowest biodiversity)
biodiv_centre_sums <- biodiv_centre_sums%>% arrange(Cluster)

## Update Cluster labels
levels(biodiv_centre_sums$Cluster)[levels(biodiv_centre_sums$Cluster) == "2"] <- "Bio-A"
levels(biodiv_centre_sums$Cluster)[levels(biodiv_centre_sums$Cluster) == "8"] <- "Bio-B"
levels(biodiv_centre_sums$Cluster)[levels(biodiv_centre_sums$Cluster) == "5"] <- "Bio-B"
levels(biodiv_centre_sums$Cluster)[levels(biodiv_centre_sums$Cluster) == "1"] <- "Bio-D"
levels(biodiv_centre_sums$Cluster)[levels(biodiv_centre_sums$Cluster) == "4"] <- "Bio-E"
levels(biodiv_centre_sums$Cluster)[levels(biodiv_centre_sums$Cluster) == "6"] <- "Bio-F"
levels(biodiv_centre_sums$Cluster)[levels(biodiv_centre_sums$Cluster) == "7"] <- "Bio-G"
levels(biodiv_centre_sums$Cluster)[levels(biodiv_centre_sums$Cluster) == "3"] <- "Bio-H"

## Add in biodiversity values (From hotspots paper - see Table 4)
biodiv_centre_sums$Biodiversity <- c(4.83,4.34,4.27,4.06,3.85,3.28,3.14,2.39)
head(biodiv_centre_sums)

## Add colouring to Cluster column in kable output (optional)
biodiv_centre_sums$Cluster = cell_spec(
  biodiv_centre_sums$Cluster, #color = "black", align = "c", angle = 45, 
  background = factor(biodiv_centre_sums$Cluster, c("Bio-A","Bio-B","Bio-C","Bio-D","Bio-E","Bio-F","Bio-G","Bio-H"), c(
    "#37C331", "#6FD326" ,"#A8E21B", "#E0F210", "#F3F223","#E2E256",  "#D1D189",  "#C0C1BC")
  ))

## Define the path to save the HTML file
html_file <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_biodiversity_values.html"

## Save table
kable(biodiv_centre_sums, escape=FALSE,format = "html",align = c("l", "c"),caption = "") %>%#Table 2. Biodiversity values of each macrofaunal biodiversity cluster group (Cooper et al., submitted), derived from the ranked sum of cluster centres.
  column_spec (1,color = "black", bold = F,width = "12em")%>%
  column_spec(2, width = "12em")%>%  # Adjust the width as needed 
  kable_styling(full_width = FALSE)%>% #"striped", "condensed ,bootstrap_options = c("hover")
  #add_header_above(c( "Table 2. Biodiversity by cluster group. Values based on summed cluster centres (from Cooper et al. (2025)." = 2))%>%
  
  save_kable(file = html_file)

# Capture the screenshot of the HTML file
png_file <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_biodiversity_values.png"
webshot(html_file, file = png_file)

# Read the PNG file
img <- image_read("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_biodiversity_values.png")

# Trim the white space
img_trimmed <- image_trim(img)

# Save the trimmed image
image_write(img_trimmed, "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_2.png")
#_______________________________________________________________________________
#### BIODIVERSITY: SUBSTITUTE (RASTER CLUSTERS WITH CLUSTER CENTRE SUM) ####

## Initialize a data frame to store results
biodiv_area_df <- data.frame(value = unique_values, area = NA)

## Update values
biodiv_model_rast1 <- subst(biodiv_model_rast, from = 1, to = 4.06)
biodiv_model_rast2 <- subst(biodiv_model_rast1, from = 2, to = 4.83)
biodiv_model_rast3 <- subst(biodiv_model_rast2, from = 3, to = 2.39)
biodiv_model_rast4 <- subst(biodiv_model_rast3, from = 4, to = 3.85)
biodiv_model_rast5 <- subst(biodiv_model_rast4, from = 5, to = 4.27)
biodiv_model_rast6 <- subst(biodiv_model_rast5, from = 6, to = 3.28)
biodiv_model_rast7 <- subst(biodiv_model_rast6, from = 7, to = 3.14)
biodiv_model_rast8 <- subst(biodiv_model_rast7, from = 8, to = 4.34)

#plot(biodiv_model_rast)
#_______________________________________________________________________________
#### BIODIVERSITY: MAP (CLUSTER CENTRE SUM) ####

## Load package
library(viridis)

## Create plot
biodiversity_numeric_plot <- ggplot() +
  geom_spatraster(data = biodiv_model_rast8) +
  scale_fill_viridis(na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 22)+
  #annotate("text",x=c(-9.52),y=c(60.4),label=c("Combined Risk"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 22))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.20))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.key.size = unit(0.8, "cm"))+#t, r, b, l
  theme(
    axis.title.y = element_text(colour = "white"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x=element_text(colour = "white")
  )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

biodiversity_numeric_plot

## Save numeric biodiversity map
ggsave(plot = biodiversity_numeric_plot,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\biodiversity_numeric_plot.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
#### BIODIVERSITY: RESCALE (CLUSTER CENTRE SUM) ####

## Rescale values of SpatRaster object
biodiv_model_rast_numeric_scaled <- rescale01(biodiv_model_rast8)
plot(biodiv_model_rast_numeric_scaled)

# Save as .tiff
writeRaster(biodiv_model_rast_numeric_scaled,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\biodiv_numeric_scaled.tif',overwrite=TRUE)#,format = "GTiff"
#_______________________________________________________________________________
#### RISK LAYER ####
#_______________________________________________________________________________
#### RISK: LOAD SCALED RISK ELEMENT MODELS ####

## Load rasters
sensitivity_model_scaled <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Sensitivity/Model/sensitivity_Mean_SQ_scaled.tif')
assemblage_model_rast_scaled <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Assemblage rarity/Model/assemblage_rarity_scaled.tif')
biodiv_model_rast_numeric_scaled <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Biodiversity/Model/biodiv_numeric_scaled.tif')
#_______________________________________________________________________________
#### RISK: CREATE LAYER (SUMMATION OF NORMALISED ELEMENT LAYERS) ####

## Sum the three layers
r_cont_risk <- sensitivity_model_scaled + assemblage_model_rast_scaled + biodiv_model_rast_numeric_scaled
plot(r_cont_risk)# Note zero and NAs are same

#_______________________________________________________________________________
#### RISK: RESCALE (SUM OF NORMALISED ELEMENTS) ####
## Rescale values of SpatRaster object
r_cont_risk_scaled <- rescale01(r_cont_risk)

# Check the current name of the layer
print(names(r_cont_risk_scaled))

# Change the name of the layer
names(r_cont_risk_scaled) <- "combined_risk"

# Verify the new name
r_cont_risk_scaled

## Save as .tiff
writeRaster(r_cont_risk_scaled,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\combined_risk_scaled.tif',overwrite=TRUE)#,format = "GTiff"
#_______________________________________________________________________________
#### RISK: MAP (NORMALISED) ####

## Load package
library(viridis)

## Load raster
r_cont_risk_scaled <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Risk/Model/combined_risk_scaled.tif')

## Risk map
risk_plot <- ggplot() +
  geom_spatraster(data = r_cont_risk_scaled) +
  scale_fill_viridis(na.value="transparent")+
   geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 18)+
  annotate("text",x=c(-9.52),y=c(60.4),label=c("Combined Risk"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 18))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.15))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.key.size = unit(1.0, "cm"))+#t, r, b, l
  theme(
    axis.title.y = element_text(colour = "white"),
    axis.title.x = element_text(colour = "white")
  )+
  #theme(legend.title = element_text( size=2), legend.text=element_text(size=20))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

risk_plot

## Save assemblage map
ggsave(plot = risk_plot,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\combined_risk_plot.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,

#_______________________________________________________________________________
#### BIODIVERSITY: MAP (NORMALISED)  ####

# Load biodiv norm raster
biodiv_norm <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Biodiversity/Model/biodiv_numeric_scaled.tif')

## Create plot
biodiv_norm_plot <- ggplot() +
  geom_spatraster(data = biodiv_norm) +
  scale_fill_viridis(na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw(base_size = 18)+
  annotate("text",x=c(-7.5),y=c(60.05),label=c("Biodiversity"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10),legend.title = element_blank())+#make legend background transparent and text white
  theme(plot.margin = unit(c(0,0.2,0,1), "cm"),legend.key.size = unit(0.45, "cm"))+#t, r, b, l
  theme(legend.position=c(0.87,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(
    axis.title.x = element_text(colour = "white"),
    axis.title.y = element_text(colour = "white"),
    axis.text.x = element_blank()
  )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### SENSITIVITY: MAP (NORMALISED) ####

# Load biodiv norm raster
sens_norm <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Sensitivity/Model/sensitivity_Mean_SQ_scaled.tif')

## Create plot
sens_norm_plot <- ggplot() +
  geom_spatraster(data = sens_norm) +
  scale_fill_viridis(na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw(base_size = 18)+
  annotate("text",x=c(-7.8),y=c(60.05),label=c("Sensitivity"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10),legend.title = element_blank())+#make legend background transparent and text white
  theme(plot.margin = unit(c(0,0.2,0,1), "cm"),legend.key.size = unit(0.45, "cm"))+#t, r, b, l
  theme(legend.position=c(0.87,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(
    axis.title.x = element_text(colour = "white"),
    axis.title.y = element_text(colour = "white"),
    axis.text.x = element_blank()
  )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### ASSEMBLAGEs: MAP (NORMALISED) ####

## Load raster
assem_norm <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Assemblage rarity/Model/assemblage_rarity_scaled.tif')

## Create plot
assem_norm_plot <- ggplot() +
  geom_spatraster(data = assem_norm) +
  scale_fill_viridis(na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw(base_size = 18)+
  annotate("text",x=c(-9),y=c(60.05),label=c("Rarity"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10),legend.title = element_blank())+#make legend background transparent and text white
  theme(plot.margin = unit(c(0,0.2,0,1), "cm"),legend.key.size = unit(0.45, "cm"))+#t, r, b, l
  theme(legend.position=c(0.87,0.21))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(colour = "white"),
        axis.title.x = element_text(colour = "white"))
#_______________________________________________________________________________
### NORMALISED RISK AND ELEMENTS PLOT (FIGURE 5) ####

## Load package
library(ggpubr)

## Combined Risk Plot
risk_plot2 <- ggpubr::ggarrange(risk_plot,   
                                nrow=1, font.label=list(color="black",size=16,face='plain'),align="v")#labels = c("d)"),

## Risk elements and combined risk plots
figure4 <-ggarrange(
  ggpubr::ggarrange(sens_norm_plot,biodiv_norm_plot,assem_norm_plot ,nrow=3,font.label=list(color="black",size=18,face='plain'),align="v",heights = c(0.94,0.94, 1)),
  risk_plot2,
  nrow = 1, 
  widths = c(0.39,1),
  #labels = c("a)","b)")
  labels = c("","")
)

## Add Lat and Lon labels
figure4_2 <- annotate_figure(
  figure4, 
  bottom = text_grob("          Longitude", 
                     color = "black", face = "plain", size = 18),#,
  left = text_grob("Latitude", 
                   color = "black", face = "plain", size = 18,rot = 90)
)

## Save as png
ggsave(plot = figure4_2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\Figure_5.png"),
       height = 310, width =450, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285

#_______________________________________________________________________________
#### RISK: CONFIDENCE ####

## Load libraries
library(terra)

## Load rasters
#sensitivity_mean <- rast("C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Sensitivity/Model/sensitivity_Mean_SQ.tif")   # predicted sensitivity (mean across runs)
biodiv_conf <-  rast("C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Bio (Combined Biodiversity)/Confidence/BiodiversityClusterConfidence_Aug25.tif")# coefficient of variation (CV)
rarity_conf <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Assemblage/Confidence/AssemblageConfidence_Nov24.tif')# biodiversity confidence (0–1)
sensitivity_cv <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Sensitivity/Confidence/sensitivity_CV_SQ.tif')# rarity confidence (0–1)

# --- Step 1: Convert sensitivity CV -> Confidence ---
# CV on a scale from 0 to 1 where 0 is high
# This is the opposite for categorical models, so need to subtract CV values from 1 to ensure comparability
sensitivity_conf <- 1 - sensitivity_cv

# --- Step 2: Combine confidence layers to get an avaerage
risk_conf <- 
  (sensitivity_conf +
    biodiv_conf +
    rarity_conf)/3

# --- Step 3: Optional classification into 5 confidence categories ---
risk_conf_class <- classify(
  risk_conf,
  rcl = matrix(c(
    0.0, 0.2, 1,  # Low 
    0.2, 0.4, 2,  # Low - Moderate
    0.4, 0.6, 3,  # Moderate
    0.6, 0.8, 4,  # High - Moderate
    0.8, 1.0, 5   # High
  ), ncol = 3, byrow = TRUE),
  include.lowest = TRUE
)

levels(risk_conf_class) <- data.frame(
  ID = 1:5,
  Class = c("Low","Low - Moderate","Moderate", "High - Moderate", "High")
)

plot(risk_conf_class)
# --- Step 7: Save outputs ---
writeRaster(risk_conf, "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicPOSEIDON_Risk/R/www/risk_confidence_continuous.tif", overwrite = TRUE)
writeRaster(risk_conf_class, "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicPOSEIDON_Risk/R/www/risk_confidence_5class.tif", overwrite = TRUE)

## Plot rasters
plot(risk_conf)
plot(risk_conf_class)
#_______________________________________________________________________________
#### RISK: MAP (CONFIDENCE) ####

## Load risk confidence raster
risk_conf_rast <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Risk/Confidence/risk_confidence_continuous.tif')

## Combined Risk confidence map
risk_conf <- ggplot() +
  geom_spatraster(data = risk_conf_rast) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name="Confidence")+#Oranges
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 18)+
  annotate("text",x=c(-9.52),y=c(60.4),label=c("Combined Risk"),color="#36454F", size=5)+#3
  theme(legend.background=element_blank(),legend.text = element_text(color="white"),legend.title = element_text(color = "white"),size= 20)+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.21))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.key.size = unit(1.0, "cm"))+#t, r, b, l
  theme(
    axis.title.y = element_text(colour = "white"),
    axis.title.x = element_text(colour = "white"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(color = "white", size = 18))

risk_conf

## Save
ggsave(plot = risk_conf_plot,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\risk_conf_plot.png"),
       height = 200, width =200,units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#height = 65, width =180,
#_______________________________________________________________________________
### CONFIDENCE: MAP (SENSITIVITY, BIODIVERSITY, RARITY, RISK - 2X2 FORMAT) ####

## Load package
library(ggpubr)

## Check plots are ggplot
class(biodiv_conf)
class(sensitivity_conf)
class(assemblage_conf)
class(risk_conf)

## Create combined plot
biodiv_sens_assem_risk_conf<- ggpubr::ggarrange(biodiv_conf,sensitivity_conf,assemblage_conf,risk_conf,ncol=2,nrow=2,labels = c("a)", "b)","c)","d)"),font.label = list(size = 24), label.x = 0.02,label.y = 1.05 )

top_row <- ggarrange(
  biodiv_conf, sensitivity_conf,
  ncol = 2,
  labels = c("a)", "b)")#,
  #font.label = list(size = 24),
  #label.y = 1.02  # slightly above for top row
)

bottom_row <- ggarrange(
  assemblage_conf, risk_conf,
  ncol = 2,
  labels = c("c)", "d)")#,
  # font.label = list(size = 24),
  #label.y = 1.05  # even higher to avoid y-axis overlap
)

biodiv_sens_assem_risk_conf <- ggarrange(
  top_row, bottom_row,
  ncol = 1, nrow = 2
)

risk_conf_plots <- annotate_figure(
  biodiv_sens_assem_risk_conf, 
  bottom = text_grob("      Longitude", color = "black", face = "plain", size = 24),
  left = text_grob("Latitude", color = "black", face = "plain", size = 24, rot = 90))

## Save
ggsave(plot =risk_conf_plots,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\sens_rarity_biodiv_risk_conf_plot_4by4.png"),
       height = 420, width =420, units = "mm", dpi = 500,
       
       device = "png",limitsize = FALSE,bg="white")#width =285
#_______________________________________________________________________________
### CONFIDENCE: MAP (SENSITIVITY, BIODIVERSITY, RARITY, RISK - 3 & 1 FORMAT; FIGURE 6) ####

## Load package
library(ggpubr)

## Combined Risk Plot
risk_conf_plot2 <- ggpubr::ggarrange(risk_conf,   
                                     nrow=1, font.label=list(color="black",size=16,face='plain'),align="v")#labels = c("d)"),

## Risk elements and combined risk plots
figure4_conf <-ggarrange(
  ggpubr::ggarrange(sensitivity_conf,biodiv_conf,assemblage_conf ,nrow=3,font.label=list(color="black",size=18,face='plain'),align="v",heights = c(0.94,0.94, 1)),
  risk_conf_plot2,
  nrow = 1, 
  widths = c(0.39,1),
  #labels = c("a)","b)")
  labels = c("",""))

## Add Lat and Lon labels
figure4_conf2 <- annotate_figure(
  figure4_conf, 
  bottom = text_grob("          Longitude", 
                     color = "black", face = "plain", size = 18),#,
  left = text_grob("Latitude", 
                   color = "black", face = "plain", size = 18,rot = 90))

## Save as png
ggsave(plot = figure4_conf2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\Figure_6.png"),
       height = 310, width =450, units = "mm", dpi = 500,
       device = "png",limitsize = FALSE,bg="white")#width =285
#_______________________________________________________________________________
### RISK ELEMENTS: MAP (MODELS AND NUMERIC DERIVATIVES; FIGURE 4) ####

## Load packages
library(ggpubr)
library(egg)

## Blank plot
blank_plot <- ggplot() + 
  geom_blank(data=countries) + 
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme_void()+
  theme(
    axis.title.y = element_text(colour = "white"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "white")#,
    #axis.text.x=element_text(colour = "white")
  )

## Add an arrow between the plots
arrow_plot <- ggplot() +
  geom_segment(aes(x=1, y=1, xend=3, yend=1), 
               arrow = arrow(length = unit(0.5, "cm")), color = "black",lineend = "round", # See available arrow types in example above
               linejoin = "round",
               linewidth = 1.5) +
  xlim(0, 3) + 
  theme_void()+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme_void()+
  theme(
    axis.title.y = element_text(colour = "white"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "white")#,
    #axis.text.x=element_text(colour = "white")
  )
arrow_plot 

white_arrow_plot <- ggplot() +
  geom_segment(aes(x=1, y=1, xend=3, yend=1), 
               arrow = arrow(length = unit(0.5, "cm")), color = "white") +
  xlim(0, 3) + 
  theme_void()+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme_void()+
  theme(
    axis.title.y = element_text(colour = "white"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "white")#,
    #axis.text.x=element_text(colour = "white")
  )
class(assemblage_model)

## Create each row
col1 <- ggpubr::ggarrange(sensitivity_model,biodiv_model, assemblage_model, labels = c("a)", "b)","c)"),nrow=3,font.label = list(size = 24,face = "plain"))#ggpubr
col2 <-ggpubr::ggarrange(white_arrow_plot ,arrow_plot,arrow_plot, labels = c("", "",""),nrow=3,font.label = list(size = 24,face = "plain"))#ggpubr
col3 <- ggpubr::ggarrange(blank_plot,biodiversity_numeric_plot,rarity_plot, labels = c("", "d)","e)"),nrow=3,font.label = list(size = 24,face = "plain"))#ggpubr blank_plot
## Stitch rows together

figure2 <- ggpubr::ggarrange(col1, col2, col3, ncol=3,font.label=list(color="black",size=22,face='plain'),align="h",widths = c(1,0.1,0.94), heights = c(1,0.94,1))#font.label=list(color="black",size=6,face='plain'),align="v") ,labels = c("D0", "D1", "D2","N")

## Add x and y labels
fig2 <- annotate_figure(
  figure2, 
  bottom = text_grob("          Longitude", 
                     color = "black", face = "plain", size = 22),#,
  left = text_grob("Latitude", 
                   color = "black", face = "plain", size = 22,rot = 90)
)

## Save combined plot
ggsave(plot = fig2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\Figure_4.png"),
       height = 600, width =480, units = "mm", dpi = 500,
       device = "png",limitsize = TRUE,bg="white")#width =285
#_______________________________________________________________________________
#### GRAPHICAL ABSTRACT: PLOTS (FIGURE 2) ####
library(terra)
library(viridis)
library(easypackages)
library(RColorBrewer)
easypackages::packages("sf",
                       "raster",
                       "stars",
                       "r5r",
                       "geobr",
                       "aopdata",
                       "gtfs2gps",
                       "ggplot2",
                       "osmdata",
                       "h3jsr",
                       "viridisLite",
                       "ggnewscale",
                       "dplyr",
                       "magrittr",
                       prompt = FALSE
)

dem <- rast("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/assemblage_model_rast_rarity_scaled.tif")
dem2 <- aggregate(dem, fact=7,fun = modal)
dem3 <- st_as_stars(dem2)
dem4 <- st_as_sf(dem3)
plot(dem4)

sens <- rast("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/sensitivity_model_scaled.tif")
sens2 <- aggregate(sens, fact=7,fun = modal)
sens3 <- st_as_stars(sens2)
sens4 <- st_as_sf(sens3)
plot(sens4)

## Biodiversity (normalised)
biodiv <- rast("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/biodiv_model_rast_numeric_scaled.tif")
biodiv2 <- aggregate(biodiv, fact=7,fun = modal)
biodiv3 <- st_as_stars(biodiv2)
biodiv4 <- st_as_sf(biodiv3)
plot(biodiv4)

## Assemblages (raw)
assem_raw <- rast("C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Risk paper/Assemblage/Model/AssemblageMaxClass_Nov24.tif")
assem_raw2 <- aggregate(assem_raw, fact=7,fun = modal)
values(assem_raw2) <- as.factor(values(assem_raw2))## Make cluster a factor
assem_raw3 <- st_as_stars(assem_raw2)
assem_raw4 <- st_as_sf(assem_raw3)
plot(assem_raw4)

## Biodiversity (raw)
abiodiv_raw <- rast("C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Biodiversity paper/Bio (Combined Biodiversity)/Model/BiodiversityClusterMaxClass_Aug25.tif")
abiodiv_raw2 <- aggregate(abiodiv_raw, fact=7,fun = modal)
values(abiodiv_raw2) <- as.factor(values(abiodiv_raw2))## Make cluster a factor
abiodiv_raw3 <- st_as_stars(abiodiv_raw2)
abiodiv_raw4 <- st_as_sf(abiodiv_raw3)
plot(abiodiv_raw4)
#_______________________________________________________________________________
#### GRAPHICAL ABSTRACT: FUNCTIONS ####
rotate_data <- function(data, x_add = 0, y_add = 0) {
  
  shear_matrix <- function(){ matrix(c(2, 1.2, 0, 1), 2, 2) }
  
  rotate_matrix <- function(x){ 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geometry = .$geometry * shear_matrix() * rotate_matrix(pi/20) + c(x_add, y_add)
    )
}

rotate_data_geom <- function(data, x_add = 0, y_add = 0) {
  shear_matrix <- function(){ matrix(c(2, 1.2, 0, 1), 2, 2) }
  
  rotate_matrix <- function(x) { 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geom = .$geom * shear_matrix() * rotate_matrix(pi/80) + c(x_add, y_add)
    )
}
#_______________________________________________________________________________
#### GRAPHICAL ABSTRACT: PLOT (RISK ELEMENTS) ####

# annotate parameters
x = 94
color = 'gray40'
p1 <-ggplot() +
  
  # terrain
  geom_sf(data = dem4 %>% rotate_data(), aes(fill=AssemblageMaxClass_Nov24), color=NA, show.legend = FALSE) +
  #scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  scale_fill_viridis(na.value="transparent")+
  annotate("text", label='Rarity', x=x, y= 49, hjust = 0, color=color,size=8) #+
#labs(caption = "image by @UrbanDemog")

p2 <- p1 +
  new_scale_fill() + 
  new_scale_color() +
  #ggplot() +
  geom_sf(data = sens4 %>% rotate_data(y_add = 10),aes(fill=sensitivity_Mean_SQ), color=NA, show.legend = FALSE) +
  scale_fill_viridis(na.value="transparent")+
  #scale_fill_distiller(palette = "RdYlGn", direction = 1)+ 
  annotate("text", label='Sensitivity', x=x, y= 59, hjust = 0, color=color,size=8)+ 
  
  new_scale_fill() + 
  new_scale_color() +
  #ggplot() +
  geom_sf(data = biodiv4 %>% rotate_data(y_add = 20),aes(fill=BiodiversityClusterMaxClass_Aug25), color=NA, show.legend = FALSE) +
  scale_fill_viridis(na.value="transparent")+
  #scale_fill_distiller(palette = "PuBu", direction = 1)+
  annotate("text", label='Biodiversity', x=x, y= 69, hjust = 0, color=color,size=8)+
  theme_void(
    #plot.background=element_rect(fill = "black"),
    #panel.background = element_rect(fill = 'black'),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.border = element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank()
  )+
  scale_x_continuous(limits = c(43, 108))
p2

## Save plot
ggsave(plot = p2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\raster_stack_norm.png"),
       width = 32,height = 18,units = "cm", pointsize = 48,
       device = "png",limitsize = FALSE,bg="white")

#_______________________________________________________________________________
#### GRAPHICAL ABSTRACT: MAP (OVERALL RISK) ####

## Load raster
risk <- rast("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/cont_risk.tif")

## Reduce size of raster
risk2 <- aggregate(risk, fact=7,fun = modal)

## Convert to starts object
risk3 <- st_as_stars(risk2)

## Convert to sf object
risk4 <- st_as_sf(risk3)
plot(risk4)

## Combined Risk plot
p3 <-ggplot() +
  geom_sf(data = risk4, aes(fill=combined_risk), color=NA, show.legend = FALSE) +
  scale_fill_viridis(na.value="transparent")+
  theme_void()
p3

## Save plot
ggsave(plot = p3,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\combined_risk.png"),
       width = 41,height = 40,units = "cm", pointsize = 48,
       device = "png",limitsize = FALSE,bg="white")
#_______________________________________________________________________________
#### GRAPHICAL ABSTRCT: PLOT (MODELS) ####

# Annotate parameters
x = 94
color = 'gray40'

## Assemblages
p1 <-ggplot() +
  geom_sf(data = assem_raw4 %>% rotate_data(),aes(fill=AssemblageMaxClass_Nov24), color=NA, show.legend = FALSE) +
  scale_fill_manual(values = c(  "#0000ee","#00ffff","#05aac1","#9a32cd","#00cd00","#9aff9a","#b40202","#ff0000","#ff8c00","#ffff00","#b4b404",'white'),name="cluster", 
                    labels=c("A1","A2a","A2b","B1b","C1a","C1b","D1","D2a","D2b","D2c","D2d"),na.value = "transparent", na.translate = FALSE)+
  annotate("text", label='Assemblages', x=x, y= 49, hjust = 0, color=color,size=8)

p2 <- p1 +
  
  ## Sensitivity
  new_scale_fill() + 
  new_scale_color() +
  geom_sf(data = sens4 %>% rotate_data(y_add = 10),aes(fill=sensitivity_Mean_SQ), color=NA, show.legend = FALSE) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  annotate("text", label='Sensitivity', x=x, y= 59, hjust = 0, color=color,size=8)+ 
  
  ## Biodiversity
  new_scale_fill() + 
  new_scale_color() +
  geom_sf(data = abiodiv_raw4 %>% rotate_data(y_add = 20), aes(fill=BiodiversityClusterMaxClass_Aug25), color=NA, show.legend = FALSE) +
  scale_fill_manual(breaks = c('2','8','5','1','4','6','7','3'), values = c(
    "#37C331", "#6FD326" ,"#A8E21B", "#E0F210", "#F3F223","#E2E256",  "#D1D189",  "#C0C1BC", "white"),
    na.translate = F)+
  annotate("text", label='Biodiversity', x=x, y= 69, hjust = 0, color=color,size=8)+
  theme_void()+
  scale_x_continuous(limits = c(43, 108))
p2

## Save plot
ggsave(plot = p2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\raster_stack_models.png"),
       width = 32,height = 18,units = "cm", pointsize = 48,
       device = "png",limitsize = FALSE,bg="white")
#_______________________________________________________________________________
#### GRAPHICAL ABSTRCT: MAP (SAMPLES) ####

PSam3=ggplot()+ 
  geom_raster(data = hill_spdf, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_raster(data = dem_spdf, aes(x = x, y = y, fill = value), alpha=0.7)+
  scale_fill_gradientn(
    colours = c("#051650","#02367b","#006ca5","#0496c7", "#04bade","#55e2e9"),
    limits  = c(-3800,0),
    values  = scales::rescale(c(-3800,-1000,-100,-50,-25,0), from = c(-3800,0)))+
  geom_point(data = points2,aes(x = X, y = Y), fill="yellow",alpha = 0.7,colour="yellow",size=0.7)+
  geom_sf(data=countries, fill ="white",col ="white")+ 
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_void()+
  guides(fill = "none")

## Save plot
ggsave(plot = PSam3,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\sample_locations.png"),
       width = 41,height = 40,units = "cm", pointsize = 48,

       device = "png",limitsize = FALSE,bg="white")

