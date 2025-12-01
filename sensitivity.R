################################################################################
####      ECOLOGICAL RISK-BASED APPROACH TO FACILITATE LICENSING OFFSHORE   ####
####                        WIND DEVELOPMENT                                ####
####                                                                        ####
####                PART 1 (of 3): SENSITIVITY (VERSION 1.0, 2025)          ####
################################################################################

# This script relates to work in Bolam, S.G., Cooper, K.M., and 
# Downie, A-L. Developing an Ecological Risk-Based Approach to Facilitate
# Licensing Offshore Wind Development. Ecosphere.

# In this study we create a 'Combined Risk' layer based on 3 element layers for: 
# 'Biodiversity', 'Sensitivity' and 'Assemblage Rarity'.

# This script 'PART 1: SENSITIVITY (VERSION 1.0. 2025; see sensitivity.R)' is
# used to develop the Sensitivity layer. Whilst the code includes lines for 
# generating a random forest Sensitivity model, this is intended only as a quick
# look see. The final Sensitivity model, and associated confidence layer should 
# be created using code from Risk_ContinuousVariablesModel_2025.R, using input
# data (i.e. point sample sensitivity scores) generated in this file.

# Benthic data used in the script is sourced from the OneBenthic 
# (https://rconnect.cefas.co.uk/onebenthic_portal/) database using sql 
# queries. For users without direct access to this database, data
# can be sourced using either OneBenthic APIs:
# Faunal data: https://rconnect.cefas.co.uk/onebenthic_api_1/__docs__/)
# or using the OneBenthic Data Extraction tool: Grab/Core
# (https://rconnect.cefas.co.uk/onebenthic_dataextractiongrabcore/)
#_______________________________________________________________________________
#### SENSITIVITY: LOAD REQUIRED DATA ####

## Set working directory
setwd("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R")

## Load packages
library(pool)
library(DBI)
library (RPostgres)
library(dplyr)

## Create connection to OneBenthic
Sys.setenv(R_CONFIG_ACTIVE = "one_benthic")

dw <- config::get()

pool <- dbPool(drv = dbDriver(dw$driver),
               dbname = dw$database,
               host = dw$server,
               port =  dw$port,
               user = dw$uid,
               password = dw$pwd)

## Get sensitivity scores by aphia id ( from View 'taxon_sensitivity_score_phy_generic')
traits = dbGetQuery(pool,"SELECT 
                    worrms_aphiaid as aphiaid,
                    taxon_sensitvity_score_phy_generic as sensitivity
                    FROM
                    faunal_data.taxon_sensitivity_score_phy_generic;")

## Inspect sensitivity scores
head(traits)

## Find min and max sensitivity scores
min(traits$sensitivity)#20
max(traits$sensitivity)#59

## Get taxon abundance data by sample
taxa_by_st = dbGetQuery(pool,"
SELECT
s.samplecode,
s.samplelat,
s.samplelong,
ts.worrms_aphiaid as aphiaid,
ts.abund

FROM 
associations.survey as su
INNER JOIN associations.surveysample as ss ON ss.survey_surveyname = su.surveyname 
INNER JOIN samples.sample as s ON ss.sample_samplecode = s.samplecode
INNER JOIN faunal_data.taxasample as ts ON s.samplecode= ts.sample_samplecode 
INNER JOIN faunal_data.worrms as w ON w.aphiaid = ts.worrms_aphiaid 

WHERE (s.gear_gearcode = 'MHN' OR
s.gear_gearcode = 'DG' OR
s.gear_gearcode = 'VV' OR
s.gear_gearcode = 'SM' OR
s.gear_gearcode = 'NIOZ' OR
s.gear_gearcode = 'BC_0.1' OR
s.gear_gearcode = 'C/VV'OR
s.gear_gearcode = 'BC' OR
s.gear_gearcode = 'DVV' OR
s.gear_gearcode = 'DG/VV')
AND (ts.taxaqual_qualifier NOT IN ('J','E','EP', 'L', 'MEGA', 'PR', 'PU', 'TAIL', 'Z', 'SP', 'FRAG','FUR', 'DEAD') OR ts.taxaqual_qualifier IS NULL)
AND (s.treatment = 'R' or s.treatment IS NULL)
AND s.macrosieve= 1
AND w.include = TRUE
AND s.samplelat > 47.92938
AND ts.abund IS NOT NULL
AND s.id <= 52951
AND su.surveyname NOT IN ('Long-Term Monitoring Program Muddy-sandy intertidal flats in Kandalaksha Bay (White Sea)')
--AND su.datapubliclyavailable = TRUE
ORDER by s.samplecode,ts.abund desc;")

## Inspect data
head(taxa_by_st)
dim(taxa_by_st)# 1194355       5

## Examine both dataframes
names(traits)
names(taxa_by_st)

## Merge dataframes
df_merge <- merge(taxa_by_st, traits, by = 'aphiaid', all = F)
head(df_merge)

## Transform abundances (sqrt)
df_merge$abund_trans <- sqrt(df_merge$abund)

## Take cols of interest
df_merge2 <- df_merge[,c(2,3,4,1,6,7)]# samplecode, samplelat, samplelong, aphiaid, sensitvity, abund_trans
head(df_merge2)

## Order by samplecode
df_merge3 <- df_merge2[order(df_merge2$samplecode),]
head(df_merge3)

## Add column 'sensitvity_abund_trans' for 'sensitivity' multiplied by 'abund_trans'
df_merge3$sensitvity_abund_trans <- df_merge3$sensitivity * df_merge3$abund_trans
head(df_merge3)
View(df_merge3)

## Calculate mean individual sensitivity score (sum of 'sensitvity_abund_trans' across all taxa divided by total number of individuals in the sample)
library(dplyr)
df_merge4 <-  df_merge3 %>% group_by(samplecode,samplelong,samplelat)%>% summarise(sum_indiv = sum(abund_trans), score =sum(sensitvity_abund_trans),sample_sensitivity = score/sum_indiv)
head(df_merge4)

## Remove unwanted cols (i.e.'sum_indiv' and 'score')
df_merge5 <- df_merge4[,c(1:3,6)]
head(df_merge5)
dim(df_merge5)#37925     4
#_______________________________________________________________________________
#### SENSITIVITY: SAMPLE LOCATIONS (FIGURE 1)  ####

## Load libraries
library(sf)
library(rasterVis)# use raster in ggplot
library(raster)
library(ggnewscale)
library(scales)
library(ggplot2)
library(dplyr)
library(shadowtext)

## Load countries polygon and Norther Ireland border
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))
ni_border <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\ni_border.shp"))

## Load DEM
dem <- raster("C:\\Users\\kmc00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\GEBCO_15_Nov_2023_be52dc9c9b2d\\gebco_2023_n61.0_s48.0_w-11.0_e11.0.tif")

## Reduce size of DEM (optional)
#dem2 <- aggregate(dem, fact=5)
#dem <- dem2

## Remove elevations above sea level.
dem[dem>0] <- 0

## Create slope and hillshade
slope = terrain(dem, opt='slope')
aspect = terrain(dem, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)

dem_spdf <- as(dem, "SpatialPixelsDataFrame")
dem_spdf <- as.data.frame(dem_spdf)
colnames(dem_spdf) <- c("value", "x", "y")

hill_spdf <- as(hill, "SpatialPixelsDataFrame")
hill_spdf <- as.data.frame(hill_spdf)
colnames(hill_spdf) <- c("value", "x", "y")

## Get unique position coordinates
unique_pos <- df_merge5 %>% dplyr::select(samplecode,samplelong, samplelat)%>%unique

## Convert latitude and longitude into geometries using st_as_sf(). 
points <- unique_pos %>%
  st_as_sf(coords = c("samplelong", "samplelat"), crs = 4326)

## st_coordinates() extracts the lon/lat values as a data frame with X and Y columns so you can use geom_point (ability to resize points etc)
points2 <- st_coordinates(points)
points2 <- as.data.frame(st_coordinates(points))
dim(points2)
## Create plot (map of sample locations, place names and background bathymetry)
# For bathy colours and breakpoints: https://stackoverflow.com/questions/70739780/is-there-a-scale-function-with-which-i-can-use-4-breaks-points

PSam2=ggplot()+ 
  geom_raster(data = hill_spdf, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_raster(data = dem_spdf, aes(x = x, y = y, fill = value), alpha=0.7)+
  scale_fill_gradientn(
    colours = c("#051650","#02367b","#006ca5","#0496c7", "#04bade","#55e2e9"),
    limits  = c(-3800,0),
    values  = scales::rescale(c(-3800,-1000,-100,-50,-25,0), from = c(-3800,0)))+
  geom_point(data = points2,aes(x = X, y = Y), fill="yellow",alpha = 0.4,colour="yellow",size=0.1)+
  geom_sf(data=countries, fill ="black",col ="black")+ 
  geom_sf(data=ni_border, ,col ="grey",linewidth = 0.2)+ 
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw()+#text = element_text(size=40)
  xlab("Longitude") +
  ylab("Latitude")+
  annotate("text",x=c(-1.3),y=c(52.5),label=c("UNITED \nKINGDOM"),color="white", size=5)+#3
  annotate("text",x=c(-6.7),y=c(51.2),label=c("Celtic \nSea"),color="white", size=7)+
  annotate("text",x=c(-1.1),y=c(50.15),label=c("English Channel"),color="white", size=7)+
  annotate("text",x=c(-4.7),y=c(53.8),label=c("Irish Sea"),color="white", size=6)+
  annotate("text",x=c(3),y=c(56.5),label=c("North Sea"),color="white", size=7)+
  annotate("text",x=c(6.3),y=c(54.5),label=c("German Bight"),color="white", size=5)+
  annotate("text",x=c(4.3),y=c(54.3),label=c("Oyster Ground"),color="white", size=5)+
  annotate("text",x=c(2.9),y=c(51.9),label=c("Southern \nBight"),color="white", size=5)+
  annotate("text",x=c(0),y=c(58.5),label=c("Fladden Ground"),color="white", size=5)+
  annotate("text",x=c(2.33),y=c(54.9),label=c("Dogger Bank"),color="white", size=5)+
  annotate("text",x=c(5.6),y=c(58.3),label=c("Norwegian Trench"),color="white", size=5, angle = -40)+
  annotate("text",x=c(-1),y=c(56.5),label=c("Scalp \nBank"),color="white", size=5)+
  annotate("text",x=c(3),y=c(58),label=c("Ling Bank"),color="white", size=5)+
  annotate("text",x=c(2.5),y=c(59),label=c("Utsira \nHigh"),color="white", size=5)+
  annotate("text",x=c(5.1),y=c(57),label=c("Fisher \nBanks"),color="white", size=5)+
  annotate("text",x=c(7.5),y=c(57),label=c("Jutland \nBank"),color="white", size=5)+
  annotate("text",x=c(-5.8),y=c(52),label=c("St George's \nChannel"),color="white", size=4.5)+
  annotate("text",x=c(-7.6),y=c(53.2),label=c("IRELAND"),color="white", size=5)+
  annotate("text",x=c(1),y=c(49),label=c("FRANCE"),color="white", size=5)+
  annotate("text",x=c(6.5),y=c(52.5),label=c("NETHERLANDS"),color="white", size=5)+
  annotate("text",x=c(9.05),y=c(56),label=c("DENMARK"),color="white", size=5)+
  annotate("text",x=c(3.9),y=c(51),label=c("BELGIUM"),color="white", size=5)+
  annotate("text",x=c(7.6),y=c(59),label=c("NORWAY"),color="white", size=5)+
  annotate("text",x=c(9),y=c(53),label=c("GERMANY"),color="white", size=5)+
  annotate("text",x=c(-4.5),y=c(54.3),label=c("Isle \nof \nMan"),color="white", size=4)+
  annotate("text",x=c(-5.5),y=c(54.97),label=c("North\nChannel"),color="white", size=3)+
  annotate("text",x=c(-6.7),y=c(56.9),label=c("Hebrides"),color="white", size=5)+
  annotate("text",x=c(0.32),y=c(56.5),label=c("Devil's \nHole"),color="white", size=5)+
  annotate("text",x=c(0.1666),y=c(55.883),label=c("Swallow \nHole"),color="white", size=5)+
  annotate("text",x=c(-2.4),y=c(56.12),label=c("Firth of Forth"),color="white", size=4)+
  annotate("text",x=c(-9),y=c(50),label=c("SW \nApproaches"),color="white", size=6)+
  shadowtext::geom_shadowtext(aes(label = 'Inner \nSilver \nPit'),x=1,y=53.5, size=5,bg.color="#0496c7", bg.r=0.1, color = "white")+#hjust=1,
  shadowtext::geom_shadowtext(aes(label = 'Bristol Channel'),x=-4.8,y=51.4, size=5,bg.color="#0496c7", bg.r=0.1, color = "white")+#hjust=1,
  shadowtext::geom_shadowtext(aes(label = 'Outher \nThames'),x=1.35,y=51.57,angle = 0, size=4,bg.color="#0496c7", bg.r=0.1, color = "white")+#hjust=1,
  shadowtext::geom_shadowtext(aes(label = 'Strait of Dover'),x=1.5,y=51,angle = 45, size=5,bg.color="#0496c7", bg.r=0.1, color = "white")+#hjust=1,
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 20),legend.title = element_text(color = "white", size = 20))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.16))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.key.size = unit(1.5, "cm"))+
  theme(text = element_text(size = 22))+
  labs(fill = "Bathymetry (m)\n")

## Save Figure 1
ggsave(plot = PSam2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\Figure_1.tif"),
       width = 41,height = 40,units = "cm", pointsize = 48,
       device = "tiff",limitsize = FALSE,bg="white")
#_______________________________________________________________________________
#### SENSITIVITY: REMOVE 'REPLICATES' TO ADDRESS SPATIAL AUTOCORRELATION ####

## Load packages
library(sp)

## Set coordinates
coordinates(df_merge5) <- c("samplelong", "samplelat")

## Work out 50m distance in decimal degrees. 1 degree of latitude =111,000m
50/111000# degrees for 50m #0.0004504505

## Set distance within which to remove replicates (zero)
zd <- zerodist(df_merge5,zero = 0.0004504505)

## Drop replicates
df_merge5_norep <- df_merge5[-zd[,2], ]
dim(df_merge5_norep)# 22790     2

## Change class to df
df_merge5_norep_2=data.frame(df_merge5_norep)
class(df_merge5_norep_2)
names(df_merge5_norep_2)

## Drop col 'optional'
df_merge5_norep_3=df_merge5_norep_2[,1:(ncol(df_merge5_norep_2)-1)]
head(df_merge5_norep_3)
#View(df_merge5_norep_3)
#_______________________________________________________________________________
#### SENSITIVITY: SAVE DATA FOR MODELLING (DATA FOR ANNA) ####
df_merge5_norep_3_2 <- df_merge5_norep_3

## Add new cols
df_merge5_norep_3_2$type <- 'numeric'
df_merge5_norep_3_2$paper <- 'vulnerability'
df_merge5_norep_3_2$metric <- 'sensitivity'
head(df_merge5_norep_3_2)

## Update column names
colnames(df_merge5_norep_3_2) <- c('sample','x','y','value','type','paper','metric')
df_merge5_norep_3_2 <- df_merge5_norep_3_2[,c(6,1:3,7,5,4)]
head(df_merge5_norep_3_2)

## Save sensitivity data for use with RF modelling script (see xx.R)
write.csv(df_merge5_norep_3_2, "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/vulnerability_metric_sensitivity_4_modelling.csv", row.names=FALSE)
#_______________________________________________________________________________
#### SENSITIVITY: SPATIAL MODELLING (QUICK LOOK) ####
#_______________________________________________________________________________
#### SENSITIVITY: ENVIRONMENTAL PREDICTORS ####

## Load libraries
library(raster)

## Load rasters
bathy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/bathy3.tif")
cur <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/cur3.tif")
gravel <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Predicted_Gravel_Fraction.tif")
light <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/light3.tif")
mud <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Predicted_Mud_Fraction.tif")
oxy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/oxy3.tif")
phyto <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/phyto3.tif")
sal <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/sal3.tif")
sil <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/sil3.tif")
spm <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/SPM_MEAN.tif")
temp <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/temp3.tif")
wov <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Wave_veloc.tif")

## Update crs
crs(bathy) <- "+proj=longlat +datum=WGS84 +no_defs" 

## Create raster stack
predictors <- stack(bathy,cur,gravel,light,mud,oxy,phyto,sal,sil,spm,temp,wov)
#predictors

# Simple names for predictor variables
names(predictors)=c("Bathymetry","Current","Gravel","Light","Mud","Oxygen","Phyto","Salinity","Silicate","SPM","Temperature","WOV")#"Sand",
#names(predictors)

## Plot raster stack
#plot(predictors)
#plot(bathy)

#_______________________________________________________________________________
#### SENSITIVITY: PREPARE DATA ####

## Data for modelling
metric <- df_merge5_norep_3

## Data required for modelling sensitivity
head(metric)

## Unload extract function from tidyr package (otherwise it won't work)
#.rs.unloadPackage("tidyr")

## Extract predictor variables from raster stack
sdata <- raster::extract(predictors, metric[,2:3])

## Change from matrix to df
class(sdata)
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(metric$samplecode,metric$sample_sensitivity,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "value"
head(sdata2)
#View(sdata2)
dim(sdata2)# 22790    14

## Remove rows with NA (outside spatial extent of rasters used for modelling)
sdata3 <- na.omit(sdata2)
dim(sdata3)# 22294    14
#____________________________________________________________________________________________________________________
#### SENSITIVITY: MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Call library
#install.packages("caTools")
library(caTools)

## Vector for response variable
Y <- sdata3[,2]

## Take a random 90% of samples in proportion across the 11 groups
set.seed(2)
msk= sample.split( Y, SplitRatio = 9/10, group = NULL )

## Check it was 90% 10% split

## The training set
train = sdata3[ msk,] 
dim(train)# 20064    14
head(train)

## Remove sample labels for train
train2 =train[,2:14]#was 11
head(train2)

## The test set
test  = sdata3[ !msk,]
dim(test) # 2230   14
head(test)

## Remove station labels for test
test2 =test[,2:14]#was 11
head(test2)
str(test2)

## Check number of observations for train (TRUE) and test (FALSE) sets
print(table(Y, msk)) 

## Check number of samples in train and test sets is equal to total
dim(sdata3) # 22294    14
dim(train)+dim(test)# 22294    28
#____________________________________________________________________________________________________________________
#### SENSITIVITY: DO MODELLING ####

## Call library
library(randomForest)

## Model
model <- value ~Bathymetry+Current+Gravel+Light+Mud+Oxygen+Phyto+Salinity+Silicate+SPM+Temperature+WOV

## Prepare training data
train3 <- train2[complete.cases(train2), ]
train3$value <- as.numeric(train3$value)
str(train3)
head(train3)

## Run model
rf2 <- randomForest(model, data=train3)
#_______________________________________________________________________________
#### SENSITIVITY: EXAMINE HOW VARIABLES ARE INFLUENCING THE MODEL ####
library(gridExtra)
library(pdp)
## Produce plots 
varImpPlot(rf2)

#png('OUTPUTS/STRUCTURE/variables_affecting_model_structure.png') # height=nrow(pr), width=ncol(pr) EFFECTS TRAITS
varImpPlot(rf2)
#dev.off()

## Partial Dependence Plots
preds <- names(rf2$forest$xlevels)

for (i in 1:length(preds)){
  partialPlot(rf2, train2, preds[i], which.class ='1')
  next
}
#____________________________________________________________________________________________________________________
#### SENSITIVITY: PRODUCE FULL COVERAGE RASTER ####

## Use model to predict cluster group for each raster cell
pr <- predict(predictors, rf2)
plot(pr)

## Save as .tiff
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicCertain\\OUTPUTS\\sensitivity__quick_look.tif',overwrite=TRUE,format = "GTiff")





