################################################################################
####      ECOLOGICAL RISK-BASED APPROACH TO FACILITATE LICENSING OFFSHORE   ####
####                        WIND DEVELOPMENT                                ####
####                                                                        ####
####                PART 2: ASSEMBLAGES (VERSION 3.0, 2025)                 ####
################################################################################

# This script relates to work in Bolam, S.G., Cooper, K.M., and 
# Downie, A-L. Developing an Ecological Risk-Based Approach to Facilitate
# Licensing Offshore Wind Development. Ecosphere.

# In this study we create a Combined Risk layer based on 3 element layers for: 
# Biodiversity, Sensitivity and Assemblage Rarity.

# This script 'PART 2: ASSEMBLAGES (VERSION 3.0. 2025)' is used to develop the  
# Assemblages layer. Whilst the code includes lines for # generating a random 
# forest Assemblages model, this is intended only as a quick look-see. The
# final Assemblages model, and associated confidence layer should be created 
# using code from Risk_ClusterModel_2025.R, with input data (i.e. point sample
# assemblages) generated in this file.

# Benthic data used in the script is sourced from the OneBenthic 
# (https://rconnect.cefas.co.uk/onebenthic_portal/) database using sql 
# queries. For users without direct access to this database, data
# can be sourced using either OneBenthic APIs:
# Faunal data: https://rconnect.cefas.co.uk/onebenthic_api_1/__docs__/)
# or using the OneBenthic Data Extraction tool: Grab/Core
# (https://rconnect.cefas.co.uk/onebenthic_dataextractiongrabcore/)
#_______________________________________________________________________________
#### #### LOAD REQUIRED DATA ####

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R')

## Load packages
library(pool)
library(DBI)
library (RPostgres)
library(dplyr)
#_______________________________________________________________________________
#### LOAD SPATIAL DATA FOR USE WITH GGPLOT ####

## Load packages
library(sf)

## Bring in countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))
#_______________________________________________________________________________
#### CREATE A CONNECTION TO OneBenthic LIVE ####
Sys.setenv(R_CONFIG_ACTIVE = "one_benthic")

dw <- config::get()

pool <- dbPool(drv = dbDriver(dw$driver),
               dbname = dw$database,
               host = dw$server,
               port =  dw$port,
               user = dw$uid,
               password = dw$pwd)
#_______________________________________________________________________________
#### RETRIEVE DATA FROM ONEBENTHIC DB #### 

## SQL select query
abunddata = dbGetQuery(pool,"SELECT
su.surveyname,
s.samplecode,
s.samplelat,
s.samplelong,
w.validname,
w.family,
ts.abund,
w.phylum

FROM 
associations.survey as su
INNER JOIN associations.surveysample as ss ON ss.survey_surveyname = su.surveyname 
INNER JOIN samples.sample as s ON ss.sample_samplecode = s.samplecode
INNER JOIN faunal_data.taxasample as ts ON s.samplecode= ts.sample_samplecode 
--INNER JOIN faunal_data.taxaqual as tq ON ts.taxaqual_qualifier = tq.qualifier 
INNER JOIN faunal_data.worrms as w ON w.aphiaid = ts.worrms_aphiaid 
INNER JOIN associations.sampleowner as so ON so.sample_samplecode = s.samplecode
INNER JOIN associations.owner as o ON so.owner_ownername = o.ownername

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
AND s.id <= 52951
AND su.surveyname NOT IN ('Long-Term Monitoring Program Muddy-sandy intertidal flats in Kandalaksha Bay (White Sea)')
--AND su.datapubliclyavailable = TRUE
ORDER by su.surveyname, s.samplecode,ts.abund desc;")

## Check data
head(abunddata)

## Check the number of samples
length(unique(abunddata$samplecode))# 37933
#_______________________________________________________________________________
#### CONVERT RAW FAUNAL DATA TO FAMILY LEVEL (USE ALL FAMILIES) ####

## Load packages
library(dplyr)
library(tidyr)

## Examine data
head(abunddata)

## Make a copy
data=abunddata

# Note missing values in col 'family'. These missing values indicate taxon name is at higher level. Infill with appropriate name from col 'scientificname'
data$family <- ifelse(data$family == "", data$validname, data$family)
data$family <- ifelse( is.na(data$family), data$validname, data$family)
data$family <- ifelse( data$family == "NA", data$validname, data$family)

## Check no rows wit NA
data %>% filter(data$family == "NA")

## Collapse data by family name
data2 <- data %>% 
  dplyr::group_by(family,samplecode,samplelat,samplelong) %>% 
  dplyr::summarise(abund = sum(abund))%>%
  dplyr::ungroup()

## Check dimensions
dim(data2)# 949376      5

## Change data2 from long to wide format: data, key (headers), values
#https://www.bing.com/videos/search?q=tidyr+long+to+wide&&view=detail&mid=9397360FBC77F8AE08099397360FBC77F8AE0809&&FORM=VRDGAR
data3 <-tidyr::spread(data2,family,abund)
dim(data3) #37933   822
head(data3)

## Update column names
colnames(data3)
colnames(data3)[1]="Sample"
colnames(data3)[2]="Latitude_WGS84"
colnames(data3)[3]="Longitude_WGS84"

## Convert from tibble to a df
class(data3)
data3 <- as.data.frame(data3)

## Finally change all NAs to zero
data3[is.na(data3)] <- 0

## Change FALSE to zero
data6 <- data3 %>% mutate_if(is.logical,as.numeric)
names(data6)
#_______________________________________________________________________________
#### REMOVE 'REPLICATES' TO ADDRESS SPATIAL AUTOCORRELATION ####

## Load packages
library(sp)

## Set coordinates
coordinates(data6) <- c("Longitude_WGS84", "Latitude_WGS84")

## Work out 50m distance in decimal degrees. 1 degree of latitude =111,000m
50/111000# degrees for 50m #0.0004504505

## Set distance within which to remove replicates (zero)
zd <- zerodist(data6,zero = 0.0004504505)

## Drop replicates
data6norep <- data6[-zd[,2], ]
dim(data6norep)# 22814   820

## Change class to df
data6.2=data.frame(data6norep)
class(data6.2)
names(data6.2)

## Drop col 'optional'
data6.2=data6.2[,1:(ncol(data6.2)-2)]
#_______________________________________________________________________________
##### PREPARE DATA FOR FAUNAL ANALYSIS #####

## Faunal subset (ie remove Sample,Latitude_WGS84, Longitude_WGS84, month and year)
data7=data6.2[,4:ncol(data6.2)]

## Check dimensions of df 'data7'
dim(data7) # 22814   818

## Check df 'data7' is just the faunal data 
names(data7)# it is 

## Change class of df data7 to a matrix 
data8=data.matrix(data7) 

## Create a df 'pos' for Sample, Latitude_WGS84 and Longitude_WGS84
pos=data6.2[,1:3] 

## Check names
names(pos)

## Make df
pos <- as.data.frame(pos)
class(pos)
#_______________________________________________________________________________
#### TRANSFORM FAUNAL DATA PRIOR TO CLUSTERING #### 

## Transform the data (fourth-root transformation) 
datat=data8^(0.25)
#View(datat)
#_______________________________________________________________________________
#### ELBOW PLOT  ####

## Load packages
library(factoextra)

## Generate plot
elbow <- fviz_nbclust(datat,kmeans, method = "wss",k.max = 30,linecolor = "black")+
  geom_vline(xintercept = c(11), linetype = 2)+theme_classic(base_size = 14)
plot(elbow)

#_______________________________________________________________________________
#### KMEANS CLUSTERING #####

## Perform Kmeans clusterinig of data. Results (cluster group) to the object 'results' 
set.seed(1234) 
results=kmeans(datat,11,algorithm="MacQueen",iter.max=100,nstart=25) 

## Number of samples belonging to each cluster group 
results$size # 1004 2755 1328 1154 6590 1458 2332  619  338 2208 3028
#_______________________________________________________________________________
#### DENDROGRAM (RELATIONSHIP BETWEEN CLUSTERS) ####

## Load packages
library(ggplot2)
library(ggdendro)
library(ggplot2)
library(dplyr)
library(dendextend)

## Function to calculate absolute differences between cluster centres over all variables.
nclusters = 11
absdiff = matrix(0, nrow=nclusters, ncol=nclusters)
centers = results$centers
for (j in 1:nclusters) {
  for (k in 1:nclusters) {
    absdiff[j,k] = sum(abs(centers[j,] - centers[k,]))
  }
}
d=round(absdiff, 1)

## Find distance matrix
d1 <- dist(as.matrix(d))

## Produce dendrogram
d2 <- d1%>% hclust %>% as.dendrogram %>%set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 9) %>%  # node point size
  set("leaves_col", c('#0000ee',
                      '#05aac1',
                      '#00ffff',
                      '#b40202',
                      '#ff8c00',
                      '#ff0000',
                      '#ffff00',
                      '#b4b404',
                      '#9a32cd',
                      '#9aff9a',
                      '#00cd00'
                      )) %>%
  
  set("labels_cex",0.8)%>%
  set("labels", c( "      A1","      A2b","      A2a","      D1","      D2b","      D2a","      D2c","      D2d","      B1b","      C1b","      C1a"))%>% 
  set("branches_lwd", 0.7)

## Change dendrogram into a ggplot
ggd1 <- as.ggdend(d2)
dendrogram <-  ggplot(ggd1, horiz = T)+theme_classic(base_size = 14)+theme(axis.title.y=element_blank(),
                                                                           axis.text.y=element_blank(),
                                                                           axis.ticks.y=element_blank(),
                                                                           axis.line.y=element_blank())+labs(y='Height')
dendrogram
#_______________________________________________________________________________
#### 2.3.5 FIGURE 3 (ELBOW & DENDROGRAM) ####

## Combined elbow plot and dendrogram
png("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/elbow_dendro.png", width = 33, height = 20, units = "cm", res = 800,pointsize = 12)
ggpubr::ggarrange(elbow,NULL,dendrogram,labels = c("a)","", "b)"),nrow=1,widths = c(1, 0.05, 1))
dev.off()
#_______________________________________________________________________________
#### DF FOR PRODUCING FAUNAL CLUSTER MAPS ####

## Add cluster group from kmeans results file to df 'pos' which includes 'Sample', 
# 'Latitude_WGS84' and 'Longitude_WGS84' 
faunal.cluster=cbind(pos,results$cluster) 

## Change name of col 'results$cluster' to 'ClusterNum' 
names(faunal.cluster)[4]<-paste("ClusterNum") 

## Add a new empty col 'FaunalCluster' to df 'faunal.cluster 
faunal.cluster["FaunalCluster"]=NA 

## Populate FaunalCluster col with new names (see dendrogram from Step 21) 
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 1]<-"A2b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 2] <-"D2a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 3] <-"D1"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 4]<-"B1b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 5] <-"D2c"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 6] <-"C1b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 7] <-"D2b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 8] <-"A2a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 9]<-"A1"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 10] <-"C1a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 11] <-"D2d"

## Save file and bring in to save running above code
#write.csv(faunal.cluster,'C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/DATA/faunal.cluster.csv', row.names=F)

## Bring in cluster results file
#faunal.cluster <- read.csv2(file='C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/DATA/faunal.cluster.csv',row.names = NULL,header=T,sep= ",")
head(faunal.cluster)
#View(faunal.cluster)

## Change column type
faunal.cluster$Latitude_WGS84<- as.numeric(faunal.cluster$Latitude_WGS84) 
faunal.cluster$Longitude_WGS84 <- as.numeric(faunal.cluster$Longitude_WGS84)
faunal.cluster$FaunalCluster <- as.factor(faunal.cluster$FaunalCluster)
str(faunal.cluster)
head(faunal.cluster)
#_______________________________________________________________________________
## Save modelling data for Anna (S, N sqrt - from object univmeascoord2.mod' and cluster from obect 'faunal.cluster' )

## Take clustering result (dropping numeric cluster groups)
faunal.cluster2 <- faunal.cluster[,c(1,3,2,5)]
head(faunal.cluster2)
faunal.cluster2$type <- 'categorical'
faunal.cluster2$paper <- 'vulnerability'
faunal.cluster2$metric <- 'assemblage'
colnames(faunal.cluster2) <- c('sample','x','y','value','type','paper','metric')
faunal.cluster2 <- faunal.cluster2[,c(6,1:3,7,5,4)]
head(faunal.cluster2)
dim(faunal.cluster2)

## Save file
write.csv(faunal.cluster2, "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/poseidon_assemblages_4_modelling.csv", row.names=FALSE)

#_______________________________________________________________________________
#### OUTPUT DATA FOR SIMPER ANALYSIS IN PRIMER ####

## Add Sample code to transformed faunal data matrix
data4simper=cbind(faunal.cluster$Sample,datat)

## Change name of col 1 from 'faunal.cluster$Sample' to 'Sample'
colnames(data4simper)[1] <- "Sample"

## Check df are same length
dim(data4simper)#22814 819
dim(faunal.cluster)#22814 5

## Check both dfs have a col named 'Sample'
colnames(data4simper)
colnames(faunal.cluster)

## Load packages
library(dplyr)

## Merge dataframes by col 'Sample
simper_merged_100 <- merge(faunal.cluster[,c(1,5)],data4simper, by='Sample')
dim(simper_merged_100)#22814   820
#View(simper_merged_100)

## Remove 'Sample' column
simper_merged_100 <- simper_merged_100%>%
  select(-Sample)

## Create df for treatment
simper_merged_100_cluster <- simper_merged_100[,1]

## Create df for fauna
simper_merged_100_fauna <- simper_merged_100%>%
  select(-FaunalCluster)

## Export both objects as .csv files for use with PRIMER6
write.csv(simper_merged_100_fauna,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/DATAFORSIMPER_100.csv",row.names=FALSE)
write.csv(simper_merged_100_cluster,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/FACTORFORSIMPER_100.csv",row.names=FALSE)
#_______________________________________________________________________________
# Randomly sample 50% of the rows by group
simper_merged_50 <- simper_merged_100 %>%
  group_by(FaunalCluster) %>%
  slice_sample(prop = 0.5)%>%#summarise(count = n())
  ungroup()
dim(simper_merged_50)#11406   820

## Create df for treatment
simper_merged_50_cluster <- simper_merged_50[,1]

## Create df for fauna
simper_merged_50_fauna <- simper_merged_50%>%
  select(-FaunalCluster)
## Export both objects as .csv files for use with PRIMER6
write.csv(simper_merged_50_fauna,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/DATAFORSIMPER_50.csv",row.names=FALSE)
write.csv(simper_merged_50_cluster,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/FACTORFORSIMPER_50.csv",row.names=FALSE)
#_______________________________________________________________________________
# Randomly sample 25% of the rows by group
simper_merged_25 <- simper_merged_100 %>%
  group_by(FaunalCluster) %>%
  slice_sample(prop = 0.25)%>%#summarise(count = n())
  ungroup()
dim(simper_merged_25)#5700  820

## Create df for treatment
simper_merged_25_cluster <- simper_merged_25[,1]

## Create df for fauna
simper_merged_25_fauna <- simper_merged_25%>%
  select(-FaunalCluster)

## Export both objects as .csv files for use with PRIMER6
write.csv(simper_merged_25_fauna,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/DATAFORSIMPER_25.csv",row.names=FALSE)
write.csv(simper_merged_25_cluster,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/FACTORFORSIMPER_25.csv",row.names=FALSE)
#_______________________________________________________________________________
# Randomly sample 10% of the rows by group
simper_merged_10 <- simper_merged_100 %>%
  group_by(FaunalCluster) %>%
  slice_sample(prop = 0.1)%>%#summarise(count = n())
  ungroup()
dim(simper_merged_10)#2275  820

## Create df for treatment
simper_merged_10_cluster <- simper_merged_10[,1]

## Create df for fauna
simper_merged_10_fauna <- simper_merged_10%>%
  select(-FaunalCluster)

## Export both objects as .csv files for use with PRIMER6
write.csv(simper_merged_10_fauna,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/DATAFORSIMPER_10.csv",row.names=FALSE)
write.csv(simper_merged_10_cluster,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/FACTORFORSIMPER_10.csv",row.names=FALSE)
#_______________________________________________________________________________
#### PLOT ASSEMBLAGE CLUSTERS USING LEAFLET ####

## Load packages
library(leaflet)
library(leafem)

# Use this to match cluster group number (popup) with chr label (see https://openscience.cefas.co.uk/ob_baseline/) 
head(faunal.cluster)

## Define cluster colours
pal <- colorFactor(
  palette = c(
              "#05aac1",#1
              "#ff0000",#2
              "#b40202",#3
              "#9a32cd",#4
              "#ffff00",#5
              "#9aff9a",#6
              "#ff8c00",#7
              "#00ffff",#8
              "#0000ee",#9
              "#00cd00",#10
              "#b4b404"#11
              ),domain = faunal.cluster$'results$cluster')

## Map
leaflet() %>%
  addProviderTiles(providers$Esri.OceanBasemap,options = providerTileOptions(noWrap = TRUE))%>%
  addCircleMarkers(data=faunal.cluster,~as.numeric(Longitude_WGS84), ~as.numeric(Latitude_WGS84), popup = ~as.character(results$cluster),radius = 3,stroke = F, color = "black",weight = 1,fill = TRUE, fillColor =~pal(results$cluster),fillOpacity = 1)%>%
  setView(-3,54.6,zoom=5)%>%
  addMouseCoordinates()
#_______________________________________________________________________________
#### GGPLOT MAP OF FAUNAL CLUSTERS ####

## Load Load packages
library(ggplot2)

## Produce map
p2= ggplot()+
  geom_point(data=faunal.cluster,aes(Longitude_WGS84,Latitude_WGS84,col=FaunalCluster), size=0.45,show.legend = TRUE)+
  scale_colour_manual(values = c("#0000ee",
                                 "#00ffff",
                                 "#05aac1",
                                 "#9a32cd",
                                 "#00cd00",
                                 "#9aff9a",
                                 "#b40202",
                                 "#ff0000",
                                 "#ff8c00",
                                 "#ffff00",
                                 "#b4b404" ),name="Cluster")+ guides(colour = guide_legend(override.aes = list(size=3)))+ # Change size of legend dots
  coord_map(xlim = c(-10.7, 6),ylim = c(48, 62))+ #set x,y limits of plot 
  theme_bw(base_size = 24)+
  labs(x="Longitude",y="Latitude")

fig4a=p2+theme(legend.key.size = unit(1, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))

## Save plot to an image file (png or tiff) 
#png("OUTPUTS/FIGURE 4a.png",width = 29.7,height = 42,units = "cm", res = 600, pointsize = 48) 
#tiff("OUTPUTS/FIGURE 4a.tiff",width = 29.7,height = 42,units = "cm",res = 600,pointsize = 48) 
fig4a
#dev.off() 
#_______________________________________________________________________________
#### GGPLOT MAP OF FAUNAL CLUSTERS (FACET BY CLUSTER) ####

## Produce map 
p6= ggplot()+
  geom_point(data=faunal.cluster,aes(Longitude_WGS84,Latitude_WGS84,col=FaunalCluster), size=0.15,show.legend = FALSE)+
  coord_map(xlim = c(-10.7, 4),ylim = c(48, 62)) +#set x,y limits of plot 
  theme_bw(base_size=24)+ 
  scale_colour_manual(values = c("#0000ee","#00ffff","#05aac1","#9a32cd","#00cd00","#9aff9a","#b40202","#ff0000","#ff8c00","#ffff00","#b4b404"),name="Cluster")+ 
  guides(colour = guide_legend(override.aes = list(size=5)))+ # Change size of legend dots #(were too small) 
  labs(x="Longitude",y="Latitude")+ 
  facet_wrap(~FaunalCluster) 

fig4b=p6+guides(fill=FALSE)# remove bathy legend 

## Save plot to an image file (png or tiff) 
#png("OUTPUTS/FIGURE 4b.png",width = 29.7,height = 42,units = "cm",res = 600,pointsize = 12) 
#tiff("OUTPUTS/FIGURE 4b.tiff",width = 29.7,height = 42,units = "cm",res = 600,pointsize = 12) 
fig4b 
#dev.off()
#_______________________________________________________________________________
#### CREATE TAXON (FAMILY) NAMES WITH PHYLUM CODE ####

## Create a df for cluster centers
cluster_centers <- as.data.frame(results$centers)
#View(cluster_centers)

## See taxon names from df 'cluster_centers'
names(cluster_centers)

## Identify 'family' and 'phylum' columns from raw data
names(data)

## Take just the 'family' and 'phylum' cols
fam_phy <-unique(data[,c(6,8)])
head(fam_phy)

## Check for presence of individual families
#fam_phy[fam_phy$family == "Aspidosiphonidae", ]

## Find Phyla
unique(fam_phy$phylum)

## Create new col for phylum code
fam_phy$phy_code <- fam_phy$phylum

## Create Phylum codes
fam_phy$phy_code[fam_phy$phy_code == 'Mollusca'] <- 'Mol'
fam_phy$phy_code[fam_phy$phy_code == 'Echinodermata'] <- 'Ech'
fam_phy$phy_code[fam_phy$phy_code == 'Annelida'] <- 'Ann'
fam_phy$phy_code[fam_phy$phy_code == 'Arthropoda'] <- 'Art'
fam_phy$phy_code[fam_phy$phy_code == 'Cnidaria'] <- 'Cni'
fam_phy$phy_code[fam_phy$phy_code == 'Sipuncula'] <- 'Si'
fam_phy$phy_code[fam_phy$phy_code == 'Nemertea'] <- 'N'
fam_phy$phy_code[fam_phy$phy_code == 'Phoronida'] <- 'Phor'
fam_phy$phy_code[fam_phy$phy_code == 'Hemichordata'] <- 'Hem'
fam_phy$phy_code[fam_phy$phy_code == 'Bryozoa'] <- 'Bry'
fam_phy$phy_code[fam_phy$phy_code == 'Platyhelminthes'] <- 'Pla'
fam_phy$phy_code[fam_phy$phy_code == 'Chordata'] <- 'Cho'
fam_phy$phy_code[fam_phy$phy_code == 'Porifera'] <- 'Por'
fam_phy$phy_code[fam_phy$phy_code == 'Entoprocta'] <- 'Ent'
fam_phy$phy_code[fam_phy$phy_code == 'Brachiopoda'] <- 'Bra'
fam_phy$phy_code[fam_phy$phy_code == 'Tracheophyta'] <- 'Tra'
fam_phy$phy_code[fam_phy$phy_code == 'Ochrophyta'] <- 'Och'
fam_phy$phy_code[fam_phy$phy_code == 'Myzozoa'] <- 'Myz'
fam_phy$phy_code[fam_phy$phy_code == 'Priapulida'] <- 'Pri'
fam_phy$phy_code[fam_phy$phy_code == 'Foraminifera'] <- 'For'
fam_phy$phy_code[fam_phy$phy_code == 'Ciliophora'] <- 'Cil'
fam_phy$phy_code[fam_phy$phy_code == 'Nemertina'] <- 'Ne'
fam_phy$phy_code[fam_phy$phy_code == 'Sipunculida'] <- 'Sip'
fam_phy$phy_code[fam_phy$phy_code == 'Pogonophora'] <- 'Pog'
fam_phy$phy_code[fam_phy$phy_code == 'Nematomorpha'] <- 'Nem'
fam_phy$phy_code[fam_phy$phy_code == 'Gastrotricha'] <- 'Gas'
fam_phy$phy_code[fam_phy$phy_code == 'Coelenterata'] <- 'Coe'

## Create mew 'taxon' col which is the family name and bracketed phylum code
fam_phy$taxon <- paste(fam_phy$family," (", fam_phy$phy_code, ")", sep = "")

# Remove duplicate rows
fam_phy2 <- fam_phy %>% distinct()
dim(fam_phy2)# 825   4
head(fam_phy2)

## Remove these records (where phyla has been revised and thus there are two values)
remove_values1 <- c('Aspidosiphonidae (Si)','Cnidaria (Coe)','Golfingiidae (Si)','Siboglinidae (Pog)','NA (Mol)')
fam_phy3 <- fam_phy2 %>% filter(!taxon %in% remove_values1)

## Remove these records (phylum not present in the data)
remove_values2 <- c('Nemertina','Sipunculida')
fam_phy3 <- fam_phy3 %>% filter(!phylum %in% remove_values2)
head(fam_phy3)

## check above values deleted
fam_phy3[fam_phy3$taxon == "Cnidaria (Coe)", ]

# Create a vector of column names from df 'cluster_centers'
column_names <- as.data.frame(colnames(cluster_centers))
dim(column_names)# 818
head(column_names)

## Rename column
colnames(column_names)[1] <- 'family'

# Replace dots with spaces in the 'family' column
column_names$family <- gsub("\\.", " ", column_names$family)

# Merge data frames by 'family', keeping only rows from df 'column_names'
merged_df <- left_join(column_names, fam_phy3, by = "family")
head(merged_df)
dim(merged_df)# 818

## Update column names (family and phy code)in df 'cluster_centers'
colnames (cluster_centers) <- merged_df$taxon

#write.csv(column_names,'C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/DATA/delete_column_names2.csv', row.names=F)
#write.csv(fam_phy3,'C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/DATA/delete_fam_phy2.csv', row.names=F)
#_______________________________________________________________________________
## Save df ' cluster_centres' as 'cluster_centres_4_simper' (see later in script). This is for simper based on vegan
cluster_centers_4_simper <- cluster_centers
#_______________________________________________________________________________
#### PREPARE DATA FOR PRODUCING PHYLA PIE CHARTS ####

## Inspect df 'cluster_centers'
head(cluster_centers)
#View(cluster_centers)
#class(cluster_centers)

## Transpose df 'cluster_centers' so taxa are rows
cluster_centers_t <- as.data.frame(t(cluster_centers))
head(cluster_centers_t)

## Make column of taxon names based on row names
cluster_centers_t$taxon <- rownames(cluster_centers_t)

## Remove rownames
rownames(cluster_centers_t) <- NULL
head(cluster_centers_t)

## Update column (cluster group) names
colnames(cluster_centers_t) <- c('A2b','D2a','D1','B1b','D2c','C1b','D2b','A2a','A1','C1a','D2d','taxon')
head(cluster_centers_t)
dim(cluster_centers_t)# 818 12

## Add family and phyla info
cluster_centers_t_wt_fam <- left_join(cluster_centers_t,fam_phy2, by='taxon')
dim(cluster_centers_t_wt_fam)# 818 15
head(cluster_centers_t_wt_fam)

## Drop cols you don't need and reorder
cluster_centers_t_wt_fam2 <- cluster_centers_t_wt_fam[,c(12,14,1:11)]
cluster_centers_t_wt_fam2

## Turn data from wide to long format
data_long <- gather(cluster_centers_t_wt_fam2, cluster, count, A2b:D2d, factor_key=TRUE)
data_long

## Reorder columns
data_long2 <- data_long[,c(3,1,2,4)]
data_long2

## Create a vector specifying the number of rows to select from each group
n_vector <- c('A2b'= 57,'D2a'=24 ,'D1'= 29,'B1b'= 56,'D2c'= 8,'C1b'= 42,'D2b'= 28,'A2a'= 52,'A1'= 65,'C1a'= 30,'D2d'= 17)

## Order data by group with count from high to low 
data_long3 <- data_long2%>%
  group_by(cluster,taxon)%>%
  arrange(desc(count))
head(data_long3)

## Select number of rows by group according to above vector
library(purrr)

selected_data <- data_long3 %>%
  group_by(cluster) %>%
  group_split() %>%
  map2_dfr(names(n_vector), ~ slice(.x, 1:n_vector[.y]))
head(selected_data)

## Update column names
colnames(selected_data)[1] <- 'Cluster'

# Create a directory for the images
image_dir <- "pie_images"
dir.create(image_dir, showWarnings = FALSE)
#_______________________________________________________________________________
## Identify colours for use in Pie charts

## Load package
library(ggplot2)
library(scales)

#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(11)
hex
# "#F8766D" "#DB8E00" "#AEA200" "#64B200" "#00BD5C" "#00C1A7" "#00BADE" "#00A6FF" "#B385FF" "#EF67EB" "#FF63B6"
#_______________________________________________________________________________
#### PRODUCE PIE CHARTS FOR THE 11 DIFFERENT CLUSTER GROUPS ####
#_______________________________________________________________________________
#### PIE A1 ####
selected_data_A1 <- selected_data[which(selected_data$Cluster=='A1'),]
selected_data_A1

## Identify phyla for cluster
unique(selected_data_A1$phylum)

p_A1 <- ggplot(selected_data_A1, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
  theme(legend.position="none")

p_A1

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_A1.png", plot = p_A1, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE A2a ####
selected_data_A2a <- selected_data[which(selected_data$Cluster=='A2a'),]
selected_data_A2a

## Identify phyla for cluster
unique(selected_data_A2a$phylum)

p_A2a <- ggplot(selected_data_A2a, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_A2a

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_A2a.png", plot = p_A2a, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE A2b ####
selected_data_A2b <- selected_data[which(selected_data$Cluster=='A2b'),]
selected_data_A2b

## Identify phyla for cluster
unique(selected_data_A2b$phylum)

p_A2b <- ggplot(selected_data_A2b, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_A2b

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_A2b.png", plot = p_A2b, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE B1b ####
selected_data_B1b <- selected_data[which(selected_data$Cluster=='B1b'),]
selected_data_B1b

## Identify phyla for cluster
unique(selected_data_B1b$phylum)

p_B1b <- ggplot(selected_data_B1b, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_B1b

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_B1b.png", plot = p_B1b, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE C1a ####
selected_data_C1a <- selected_data[which(selected_data$Cluster=='C1a'),]
selected_data_C1a

## Identify phyla for cluster
unique(selected_data_C1a$phylum)

p_C1a <- ggplot(selected_data_C1a, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_C1a

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_C1a.png", plot = p_C1a, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE C1b ####
selected_data_C1b <- selected_data[which(selected_data$Cluster=='C1b'),]
selected_data_C1b

## Identify phyla for cluster
unique(selected_data_C1b$phylum)

p_C1b <- ggplot(selected_data_C1b, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_C1b

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_C1b.png", plot = p_C1b, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE D1 ####
selected_data_D1 <- selected_data[which(selected_data$Cluster=='D1'),]
selected_data_D1

## Identify phyla for cluster
unique(selected_data_D1$phylum)

p_D1 <- ggplot(selected_data_D1, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_D1

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_D1.png", plot = p_D1, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE D2a ####
selected_data_D2a <- selected_data[which(selected_data$Cluster=='D2a'),]
selected_data_D2a

## Identify phyla for cluster
unique(selected_data_D2a$phylum)

p_D2a <- ggplot(selected_data_D2a, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_D2a

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_D2a.png", plot = p_D2a, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE D2b ####
selected_data_D2b <- selected_data[which(selected_data$Cluster=='D2b'),]
selected_data_D2b

## Identify phyla for cluster
unique(selected_data_D2b$phylum)

p_D2b <- ggplot(selected_data_D2b, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_D2b

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_D2b.png", plot = p_D2b, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE D2c ####
selected_data_D2c <- selected_data[which(selected_data$Cluster=='D2c'),]
selected_data_D2c

## Identify phyla for cluster
unique(selected_data_D2c$phylum)

p_D2c <- ggplot(selected_data_D2c, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_D2c

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_D2c.png", plot = p_D2c, width = 6, height = 4)
#_______________________________________________________________________________
#### PIE D2d ####
selected_data_D2d <- selected_data[which(selected_data$Cluster=='D2d'),]
selected_data_D2d

## Identify phyla for cluster
unique(selected_data_D2d$phylum)

p_D2d <- ggplot(selected_data_D2d, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#dd.col
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_D2d

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/pie_D2d.png", plot = p_D2d, width = 6, height = 4)
#_______________________________________________________________________________
#### CONVERT PIE CHARTS TO IMAGES FOR USE IN GT TABLE ####

library(base64enc)

# Set the directory containing the .png files
directory <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/"

# List all .png files in the directory
png_files <- list.files(directory, pattern = "\\.png$", full.names = TRUE)

# Create a dataframe from the list of .png files
pie_images<- data.frame(files = png_files)

## Add column for Cluster
pie_images$Cluster <- c('A1','A2a','A2b','B1b','C1a','C1b','D1','D2a','D2b','D2c','D2d')

##Update column order
pie_images <- pie_images[,2:1]

# Print the dataframe
print(pie_images)

# Set the directory containing the .png files
directory <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/R/pie_images/"

# List all .png files in the directory
png_files <- list.files(directory, pattern = "\\.png$", full.names = TRUE)

# Function to convert a file to base64
convert_to_base64 <- function(file) {
  base64encode(file)
}

# Apply the function to all .png files and create a dataframe
base64_df <- data.frame(
  file = png_files,
  base64 = sapply(png_files, convert_to_base64)
)

# Print the dataframe
View(base64_df)

## Create a df for rownames from df 'base64_df'. Names contain the cluster group name
Cluster <- rownames(base64_df)

## Remove rownames
rownames(base64_df) <- NULL

## Add cluster col to df 'base64_df'
pie_images <- cbind(Cluster,base64_df)

## Update column names
colnames(pie_images)[2] <- 'Image'
colnames(pie_images)[3] <- 'Image_base64'
#View(pie_images)

## Update cluster groups
pie_images$Cluster <- c('A1','A2a','A2b','B1b','C1a','C1b','D1','D2a','D2b','D2c','D2d')

# Convert images to base64
pie_images2 <- pie_images %>%
  rowwise() %>%
  mutate(Image_base64 = base64enc::dataURI(file = Image, mime = "image/png"))
head(pie_images2)
#_______________________________________________________________________________
#### IDENTIFY TAXA TO INCLUDE IN CHARACTERISTICS TABLE ####

## Analyse cluster centers
head(cluster_centers)
#View(cluster_centers)

# Function to return column names in order of values for each row
get_ordered_colnames <- function(row) {
  colnames(cluster_centers)[order(-row)]
}

# Apply the function to each row
result <- apply(cluster_centers, 1, get_ordered_colnames)

# Convert the result to a data frame for better readability
result_df <- data.frame(t(result))
colnames(result_df) <- paste0("Rank_", 1:ncol(result_df))
#View(result_df)

# Vector specifying the number of columns to keep for each row (this is the mean number of taxa)
N <- c(57,24,29,56,8,42,28,52,65,30,17)

# Function to return top N column names in order of values for each row, highest first
get_top_n_colnames <- function(row, n) {
  colnames(cluster_centers)[order(-row)][1:n]
}
## Convert df 'cluster_centers' to a matrix
cluster_centers <- as.matrix(cluster_centers)

# Apply the function to each row with different N values
result <- mapply(get_top_n_colnames, split(cluster_centers, row(cluster_centers)), N, SIMPLIFY = FALSE)

# Convert the list to a data frame for better readability
result_df <-as.data.frame( do.call(rbind, lapply(result, function(x) {
  length(x) <- max(N)
  return(x)
})))

# Label the rows with row numbers
rownames(result_df) <- c('A2b','D2a','D1','B1b','D2c','C1b','D2b','A2a','A1','C1a','D2d')

# Function to convert row values to a comma-separated string
row_to_string <- function(row) {
  paste(na.omit(row), collapse = ", ")
}

# Apply the function to each row
result2 <- apply(result_df, 1, row_to_string)

# Convert the result to a data frame for better readability
result_df2 <- data.frame(Row = paste0("Row_", 1:nrow(result_df)), Values = result2)
result_df2 <- data.frame(result2)
colnames(result_df2)[1] <- 'Taxa'

## Replace cluster labels
result_df2$Cluster <- c('A2b','D2a','D1','B1b','D2c','C1b','D2b','A2a','A1','C1a','D2d')

## Remove row names
rownames(result_df2) <- NULL

# Move the last column to the first position
result_df2 <- result_df2[, c(ncol(result_df2), 1:(ncol(result_df2)-1))]
#View(result_df2)

## Get rows in correct order
chr_taxa <- result_df2[order(result_df2$Cluster, decreasing = F),]
#_______________________________________________________________________________
## UNIVARIATE SUMMARY MEASURES ####

## Load packages
library(vegan)
library(dplyr)

## Calculate univariate summary measures based on faunal abundance data in df 'data5'
Richness = specnumber(data7) # Species Richness(S)
Abundance=rowSums(data7) # Abundance

## Need cluster results before this works
uni <- cbind(pos,Richness,Abundance,results$cluster)
colnames(uni)[6] <- "Cluster"

## Calculate metrics
uni2 <- uni %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n(),r_mean = mean(Richness),r_sd=sd(Richness),a_mean=mean(Abundance),a_sd=sd(Abundance) )

## Replace cluster numbers for names
uni2$Cluster[uni2$Cluster == "1"] <- "A2b"
uni2$Cluster[uni2$Cluster == "2"] <- "D2a"
uni2$Cluster[uni2$Cluster == "3"] <- "D1"
uni2$Cluster[uni2$Cluster == "4"] <- "B1b"
uni2$Cluster[uni2$Cluster == "5"] <- "D2c"
uni2$Cluster[uni2$Cluster == "6"] <- "C1b"
uni2$Cluster[uni2$Cluster == "7"] <- "D2b"
uni2$Cluster[uni2$Cluster == "8"] <- "A2a"
uni2$Cluster[uni2$Cluster == "9"] <- "A1"
uni2$Cluster[uni2$Cluster == "10"] <- "C1a"
uni2$Cluster[uni2$Cluster == "11"] <- "D2d"

## Update column names
colnames(uni2)[1] <- "Cluster"
head(uni2)

## Order by r_mean
uni2 <- uni2[order(-uni2$r_mean),]
head(uni2)
#_______________________________________________________________________________
#### ADD UNIVARIATE RESULTS TO TABLE OF CHARACTERISING TAXA ####

## Inspect existing table of characterising taxa
#View(chr_taxa)

## Merge dataframes
uni3 <- merge(uni2,chr_taxa,by="Cluster")
#View(uni3)

# Merge the base64 images data frame with the original data
merged_data <- left_join(uni3, pie_images2, by = "Cluster")
#View(merged_data)
#View(pie_images2)
#class(pie_images2)
#str(pie_images2)
#class(uni3)
#str(uni3)

## Rename
uni3 <- merged_data
names(uni3)

## Change colname from "Image_base64" to "Phylum" as this will be the col header in gt output
colnames(uni3)[9] <- 'Phylum'
#View(uni3)
#_______________________________________________________________________________
#### CREATE NICE TABLE FOR CLUSTER SUMMARY INFO ####

## Load packages
library(gt)
library(gtExtras)
library(dplyr)
library(stringr)
library(gtsummary)
library(webshot2)
library(chromote)

# Function to make selected text bold
make_bold <- function(text, bold_words) {
  for (word in bold_words) {
    text <- stringr::str_replace_all(text, word, paste0("<b>", word, "</b>"))
  }
  return(text)
}

# Apply the function to the 'Taxa' column
uni3$Taxa <- mapply(make_bold, uni3$Taxa, list(c("Styelidae",
                                                 "Balanidae",
                                                 "Spionidae",
                                                 "Terebellidae",
                                                 "Porcellanidae",
                                                 "Syllidae",
                                                 "Polynoidae",
                                                 "Sabellariidae",
                                                 "Capitellidae",
                                                 "Serpulidae",
                                                 "Phyllodocidae",
                                                 "Nemertea",
                                                 "Cirratulidae",
                                                 "Hiatellidae",
                                                 "Alcyonidiidae",
                                                 "Escharellidae",
                                                 "Trochidae",
                                                 "Pholoidae",
                                                 "Electridae",
                                                 "Amphiuridae"),
                                               c("Sabellariidae",
                                                 "Spionidae",
                                                 "Terebellidae",
                                                 "Polynoidae",
                                                 "Nemertea",
                                                 "Phyllodocidae",
                                                 "Lumbrineridae",
                                                 "Pholoidae",
                                                 "Syllidae",
                                                 "Capitellidae",
                                                 "Cirratulidae",
                                                 "Porcellanidae",
                                                 "Semelidae"),
                                               c("Serpulidae",
                                                 "Sabellariidae",
                                                 "Syllidae",
                                                 "Spionidae",
                                                 "Terebellidae",
                                                 "Polynoidae",
                                                 "Capitellidae",
                                                 "Lumbrineridae",
                                                 "Cirratulidae",
                                                 "Amphiuridae",
                                                 "Porcellanidae",
                                                 "Nemertea",
                                                 "Maldanidae",
                                                 "Phyllodocidae",
                                                 "Styelidae",
                                                 "Escharellidae"),
                                               c("Spionidae",
                                                 "Syllidae",
                                                 "Serpulidae",
                                                 "Phyllodocidae",
                                                 "Terebellidae",
                                                 "Glyceridae",
                                                 "Nemertea",
                                                 "Polynoidae",
                                                 "Bitectiporidae",
                                                 "Capitellidae",
                                                 "Amphiuridae",
                                                 "Fibulariidae",
                                                 "Escharellidae",
                                                 "Eunicidae",
                                                 "Lumbrineridae",
                                                 "Lichenoporidae"),
                                               c("Spionidae",
                                                 "Terebellidae",
                                                 "Capitellidae",
                                                 "Syllidae",
                                                 "Serpulidae",
                                                 "Lumbrineridae",
                                                 "Cirratulidae",
                                                 "Nemertea",
                                                 "Glyceridae"),
                                               c("Spionidae",
                                                 "Capitellidae",
                                                 "Terebellidae",
                                                 "Nemertea",
                                                 "Lumbrineridae",
                                                 "Cirratulidae",
                                                 "Ampeliscidae",
                                                 "Phyllodocidae",
                                                 "Pholoidae",
                                                 "Glyceridae",
                                                 "Polynoidae",
                                                 "Ampharetidae"),
                                               c("Amphiuridae",
                                                 "Lasaeidae",
                                                 "Spionidae",
                                                 "Nephtyidae",
                                                 "Nuculidae",
                                                 "Pectinariidae",
                                                 "Phoronidae",
                                                 "Semelidae",
                                                 "Pholoidae",
                                                 "Nemertea"),
                                               c("Spionidae",
                                                 "Nemertea",
                                                 "Glyceridae",
                                                 "Terebellidae",
                                                 "Fibulariidae",
                                                 "Capitellidae"),
                                               c("Spionidae",
                                                 "Oweniidae",
                                                 "Ampharetidae",
                                                 "Capitellidae",
                                                 "Amphiuridae",
                                                 "Nemertea",
                                                 "Thyasiridae",
                                                 "Lumbrineridae",
                                                 "Nephtyidae",
                                                 "Trichobranchidae"),
                                               c("Nephtyidae",
                                                 "Spionidae",
                                                 "Opheliidae"),
                                               c("Bathyporeiidae",
                                                 "Spionidae",
                                                 "Magelonidae",
                                                 "Nephtyidae")
                                               ))

## Check text has bold markers added
head(uni3)
names(uni3)
#View(uni3)

## Drop std dev cols and update column names
 uni3 <- uni3[,c(1,2,3,5,9,7,8)]
 colnames(uni3) <- c("Cluster","n","Richness","Abundance","Phylum","Taxa","Image")
#_______________________________________________________________________________
#### Remove names not included in SIMPER output ####
 
 ## Leave only simper taxa
 library(stringr)

 ##A1
 patterns1 <- c(" Verrucidae \\(Art\\),", " Molgulidae \\(Cho\\),",
               " Ammotheidae \\(Art\\), Semelidae \\(Mol\\),",
               " Cheirocratidae \\(Art\\), Galatheidae \\(Art\\), Ascidiidae \\(Cho\\), Mytilidae \\(Mol\\),",
               " Scalibregmatidae \\(Ann\\),",
               " Bugulidae \\(Bry\\),",
               " Rissoidae \\(Mol\\), Vesiculariidae \\(Bry\\), Flustridae \\(Bry\\), Candidae \\(Bry\\), Phoronidae \\(Phor\\), Sabellidae \\(Ann\\), Calloporidae \\(Bry\\), Sertulariidae \\(Cni\\), Lumbrineridae \\(Ann\\), Crisiidae \\(Bry\\), Leptochitonidae \\(Mol\\), Maldanidae \\(Ann\\), Dorvilleidae \\(Ann\\), Pedicellinidae \\(Ent\\), Corophiidae \\(Art\\), Nuculidae \\(Mol\\), Calyptraeidae \\(Mol\\), Ampeliscidae \\(Art\\), Anomiidae \\(Mol\\), Glyceridae \\(Ann\\), Celleporidae \\(Bry\\), Callipallenidae \\(Art\\), Golfingiidae \\(Ann\\), Campanulariidae \\(Cni\\), Porifera \\(Por\\), Bitectiporidae \\(Bry\\), Veneridae \\(Mol\\), Ampharetidae \\(Ann\\), Nephtyidae \\(Ann\\), Oweniidae \\(Ann\\), Tubuliporidae \\(Bry\\), Onchidorididae \\(Mol\\), Myidae \\(Mol\\), Thoridae \\(Art\\), Lasaeidae \\(Mol\\)")
 combined_pattern1 <- paste(patterns1, collapse = "|")
 uni3$Taxa[1] <- str_remove_all(uni3$Taxa[1], combined_pattern1)
 
## A2a 
 patterns2 <- c(" Ampeliscidae \\(Art\\),",
                " Actiniaria \\(Cni\\), Amphiuridae \\(Ech\\), Nereididae \\(Ann\\), Serpulidae \\(Ann\\), Nuculidae \\(Mol\\), Ampharetidae \\(Ann\\), Scalibregmatidae \\(Ann\\), Golfingiidae \\(Ann\\), Glyceridae \\(Ann\\), Sabellidae \\(Ann\\), Dorvilleidae \\(Ann\\), Lasaeidae \\(Mol\\), Electridae \\(Bry\\), Mytilidae \\(Mol\\), Photidae \\(Art\\), Unciolidae \\(Art\\), Ammotheidae \\(Art\\), Bugulidae \\(Bry\\), Hesionidae \\(Ann\\), Orbiniidae \\(Ann\\), Myidae \\(Mol\\), Nephtyidae \\(Ann\\), Sertulariidae \\(Cni\\), Vesiculariidae \\(Bry\\), Ophiuridae \\(Ech\\), Galatheidae \\(Art\\), Molgulidae \\(Cho\\), Alcyonidiidae \\(Bry\\), Campanulariidae \\(Cni\\), Bodotriidae \\(Art\\), Flustridae \\(Bry\\), Corophiidae \\(Art\\), Oweniidae \\(Ann\\), Rissoidae \\(Mol\\), Maldanidae \\(Ann\\), Calyptraeidae \\(Mol\\), Phoxichilidiidae \\(Art\\), Trochidae \\(Mol\\)")
 combined_pattern2 <- paste(patterns2, collapse = "|")
 uni3$Taxa[2] <- str_remove_all(uni3$Taxa[2], combined_pattern2)
 
 ## A2b 
 patternsA2b <- c(" Verrucidae \\(Art\\),",
                " Sabellidae \\(Ann\\), Scalibregmatidae \\(Ann\\), Ampeliscidae \\(Art\\),",
                " Glyceridae \\(Ann\\), Myidae \\(Mol\\), Bitectiporidae \\(Bry\\), Eunicidae \\(Ann\\), Balanidae \\(Art\\), Electridae \\(Bry\\), Golfingiidae \\(Ann\\), Sertulariidae \\(Cni\\), Leptochitonidae \\(Mol\\), Nuculidae \\(Mol\\), Pholoidae \\(Ann\\), Ampharetidae \\(Ann\\), Actiniaria \\(Cni\\), Calyptraeidae \\(Mol\\), Unciolidae \\(Art\\), Calloporidae \\(Bry\\), Alcyonidiidae \\(Bry\\), Celleporidae \\(Bry\\), Lichenoporidae \\(Bry\\), Tubuliporidae \\(Bry\\), Phoronidae \\(Phor\\), Veneridae \\(Mol\\), Porifera \\(Por\\), Trochidae \\(Mol\\), Campanulariidae \\(Cni\\), Mytilidae \\(Mol\\), Bugulidae \\(Bry\\), Paraonidae \\(Ann\\), Flustridae \\(Bry\\), Pyuridae \\(Cho\\), Chorizoporidae \\(Bry\\), Photidae \\(Art\\), Dorvilleidae \\(Ann\\), Plagioeciidae \\(Bry\\), Nephtyidae \\(Ann\\), Maeridae \\(Art\\), Corophiidae \\(Art\\)")
 combined_patternA2b <- paste(patternsA2b, collapse = "|")
 uni3$Taxa[3] <- str_remove_all(uni3$Taxa[3], combined_patternA2b) 
 
 ## B1b 
 patternsB1b <- c(" Scalibregmatidae \\(Ann\\), Cirratulidae \\(Ann\\),",
                  " Galatheidae \\(Art\\),",
                  " Adeonidae \\(Bry\\), Maeridae \\(Art\\), Ampharetidae \\(Ann\\), Leptochitonidae \\(Mol\\), Glycymerididae \\(Mol\\), Microporellidae \\(Bry\\), Chorizoporidae \\(Bry\\), Hippothoidae \\(Bry\\), Golfingiidae \\(Ann\\), Sabellidae \\(Ann\\), Tubuliporidae \\(Bry\\), Phidoloporidae \\(Bry\\), Porcellanidae \\(Art\\), Bryocryptellidae \\(Bry\\), Epizoanthidae \\(Cni\\), Veneridae \\(Mol\\), Ampeliscidae \\(Art\\), Sertulariidae \\(Cni\\), Celleporidae \\(Bry\\), Porifera \\(Por\\), Cribrilinidae \\(Bry\\), Trochidae \\(Mol\\), Smittinidae \\(Bry\\), Styelidae \\(Cho\\), Clionaidae \\(Por\\), Maldanidae \\(Ann\\), Flabelligeridae \\(Ann\\), Calliopiidae \\(Art\\), Plagioeciidae \\(Bry\\), Electridae \\(Bry\\), Nereididae \\(Ann\\), Cheiloporinidae \\(Bry\\), Myidae \\(Mol\\), Pholoidae \\(Ann\\), Psammobiidae \\(Mol\\), Cirolanidae \\(Art\\), Parechinidae \\(Ech\\)")
 combined_patternB1b <- paste(patternsB1b, collapse = "|")
 uni3$Taxa[4] <- str_remove_all(uni3$Taxa[4], combined_patternB1b) 
 
 ## C1a 
 patternsC1a <- c(" Sabellariidae \\(Ann\\),",
                  " Phyllodocidae \\(Ann\\), Polynoidae \\(Ann\\), Maldanidae \\(Ann\\), Electridae \\(Bry\\), Balanidae \\(Art\\), Amphiuridae \\(Ech\\), Ampeliscidae \\(Art\\), Porcellanidae \\(Art\\), Sertulariidae \\(Cni\\), Scalibregmatidae \\(Ann\\), Escharellidae \\(Bry\\), Pholoidae \\(Ann\\), Actiniaria \\(Cni\\), Styelidae \\(Cho\\), Eunicidae \\(Ann\\), Bitectiporidae \\(Bry\\), Alcyonidiidae \\(Bry\\), Semelidae \\(Mol\\), Nereididae \\(Ann\\), Ampharetidae \\(Ann\\)")
 combined_patternC1a <- paste(patternsC1a, collapse = "|")
 uni3$Taxa[5] <- str_remove_all(uni3$Taxa[5], combined_patternC1a)
 
 ## C1b 
 patternsC1b <- c(" Serpulidae \\(Ann\\), Scalibregmatidae \\(Ann\\), Semelidae \\(Mol\\), Maldanidae \\(Ann\\),",
                  " Lasaeidae \\(Mol\\), Syllidae \\(Ann\\), Oweniidae \\(Ann\\), Pectinariidae \\(Ann\\), Sabellariidae \\(Ann\\), Phoronidae \\(Phor\\), Amphiuridae \\(Ech\\), Ophiuridae \\(Ech\\), Goniadidae \\(Ann\\), Poecilochaetidae \\(Ann\\), Fibulariidae \\(Ech\\), Paraonidae \\(Ann\\), Nephtyidae \\(Ann\\), Urothoidae \\(Art\\), Actiniaria \\(Cni\\), Eunicidae \\(Ann\\), Photidae \\(Art\\), Nuculidae \\(Mol\\), Nereididae \\(Ann\\), Golfingiidae \\(Ann\\), Dorvilleidae \\(Ann\\), Sertulariidae \\(Cni\\), Veneridae \\(Mol\\), Sabellidae \\(Ann\\), Trichobranchidae \\(Ann\\), Hesionidae \\(Ann\\)")
 combined_patternC1b <- paste(patternsC1b, collapse = "|")
 uni3$Taxa[6] <- str_remove_all(uni3$Taxa[6], combined_patternC1b)
 
 ## D1 
 patternsD1 <- c(" Oweniidae \\(Ann\\),",
                 " Cirratulidae \\(Ann\\),",
                  " Magelonidae \\(Ann\\), Ampeliscidae \\(Art\\), Pharidae \\(Mol\\), Capitellidae \\(Ann\\), Veneridae \\(Mol\\), Lumbrineridae \\(Ann\\), Phyllodocidae \\(Ann\\), Corbulidae \\(Mol\\), Cylichnidae \\(Mol\\), Orbiniidae \\(Ann\\), Sigalionidae \\(Ann\\), Thyasiridae \\(Mol\\), Flabelligeridae \\(Ann\\), Terebellidae \\(Ann\\), Glyceridae \\(Ann\\), Ampharetidae \\(Ann\\), Loveniidae \\(Ech\\)")
 combined_patternD1 <- paste(patternsD1, collapse = "|")
 uni3$Taxa[7] <- str_remove_all(uni3$Taxa[7], combined_patternD1)
 
 ## D2a 
 patternsD2a <- c(" Syllidae \\(Ann\\), Phyllodocidae \\(Ann\\), Cirratulidae \\(Ann\\), Goniadidae \\(Ann\\), Opheliidae \\(Ann\\), Lumbrineridae \\(Ann\\), Polynoidae \\(Ann\\), Oweniidae \\(Ann\\), Sigalionidae \\(Ann\\), Nephtyidae \\(Ann\\), Urothoidae \\(Art\\), Ampeliscidae \\(Art\\), Dorvilleidae \\(Ann\\), Paraonidae \\(Ann\\), Veneridae \\(Mol\\), Amphiuridae \\(Ech\\), Scalibregmatidae \\(Ann\\), Semelidae \\(Mol\\)")
 combined_patternD2a <- paste(patternsD2a, collapse = "|")
 uni3$Taxa[8] <- str_remove_all(uni3$Taxa[8], combined_patternD2a)
 
 ## D2b 
 patternsD2b <- c(" Cirratulidae \\(Ann\\), Maldanidae \\(Ann\\),",
                 " Amphinomidae \\(Ann\\), Ampeliscidae \\(Art\\), Orbiniidae \\(Ann\\),",
                 " Flabelligeridae \\(Ann\\), Glyceridae \\(Ann\\), Paraonidae \\(Ann\\), Phoronidae \\(Phor\\), Semelidae \\(Mol\\), Goniadidae \\(Ann\\), Phoxocephalidae \\(Art\\), Terebellidae \\(Ann\\), Pectinariidae \\(Ann\\), Pharidae \\(Mol\\), Leuconidae \\(Art\\), Nuculidae \\(Mol\\), Magelonidae \\(Ann\\)")
 combined_patternD2b <- paste(patternsD2b, collapse = "|")
 uni3$Taxa[9] <- str_remove_all(uni3$Taxa[9], combined_patternD2b)
 
 ## D2c 
 patternsD2c <- c(" Glyceridae \\(Ann\\), Bathyporeiidae \\(Art\\), Nemertea \\(N\\), Orbiniidae \\(Ann\\), Urothoidae \\(Art\\)")
 combined_patternD2c <- paste(patternsD2c, collapse = "|")
 uni3$Taxa[10] <- str_remove_all(uni3$Taxa[10], combined_patternD2c)
 
 ## D2d 
 patternsD2d <- c(" Tellinidae \\(Mol\\), Cirratulidae \\(Ann\\), Semelidae \\(Mol\\), Orbiniidae \\(Ann\\), Sigalionidae \\(Ann\\), Nemertea \\(N\\), Goniadidae \\(Ann\\), Oweniidae \\(Ann\\), Urothoidae \\(Art\\), Veneridae \\(Mol\\), Opheliidae \\(Ann\\), Amphiuridae \\(Ech\\), Pharidae \\(Mol\\)")
 combined_patternD2d <- paste(patternsD2d, collapse = "|")
 uni3$Taxa[11] <- str_remove_all(uni3$Taxa[11], combined_patternD2d)
 #View(uni3)
#_______________________________________________________________________________
 
## Produce gt table
 uni_tab <-  uni3%>%gt()%>%
   # Reduce number of decimal places
   #fmt_number(columns = c(Richness,	Abundance), decimals = 0)%>%

   gt_plt_bar(column = Richness, scale_type = "number",  color = "grey",text_color = "black",width=60)%>% #keep_column = TRUE,,text_color = "black",,width = 80
   gt_plt_bar(column = Abundance, scale_type = "number", color = "grey",text_color = "black",width=60) %>%#keep_column = TRUE, scale_type = "number",


   cols_hide(columns = c(n))%>%
   ## Add pie charts for phyla
   text_transform(
     locations = cells_body(columns = c(Phylum)),
     fn = function(x) {
       web_image(url = x, height = 100)
     }
   ) %>%
   cols_hide(columns = c(Image))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#b4b404")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 11)
     ))  %>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#ffff00")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 10)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#ff8c00"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 9)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#ff0000")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 8)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#b40202"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 7)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#9aff9a")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 6)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#00cd00")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 5)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#9a32cd"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 4)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#05aac1"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 3)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#00ffff")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 2)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#0000ee"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 1)
     ))%>%
 # Apply HTML rendering (bold taxon names)
 fmt_markdown(columns = vars(Taxa))%>%
   ## Add a table caption
 #  tab_header(
#     title = md("<div style='text-align: left;'><b>Table 1</b>. Cluster group characteristics including mean richness, mean abundance, phylytic composition (pie charts) and characterising taxa (family level or above). Listed taxa are those with the highest mean centroid values, with the number of taxa reported based on the group mean richness value. Highlighted taxa are those identified by a SIMPER analysis as contributing to ~50% of the similarity between samples. Phyla codes are given in parenthesis (see Legend at foot of table).<br><br> </div>")
#   )%>%#: 'Ann' - Annelida, 'Art' - Arthropoda, 'Bra' - Brachiopoda, 'Bry' - Bryozoa, 'Cil' - Ciliophora, 'Cho' - Chordata, 'Coe' - Coelenterata, 'Cni' - Cnidaria, 'Ech' - Echinodermata, 'Ent' - Entoprocta, 'For' - Foraminifera, 'Gas' - Gastrotricha, 'Hem' - Hemichordata, 'Mol' - Mollusca, 'Myz' - Myzozoa, 'N' - Nemertea, 'Ne' - Nemertina, 'Nem' - Nematomorpha, 'Och' - Ochrophyta, 'Phor' - Phoronida, 'Pla' - Platyhelminthes, 'Pog' - Pogonophora, 'Por' - Porifera, 'Pri' - Priapulida, 'Si' - Sipuncula, 'Sip' - Sipunculida, 'Tra' - Tracheophyta)
   tab_options(
     table.border.top.style = "hidden"
     #heading.border.bottom.style = "hidden"
   )%>%tab_options(
     heading.title.font.size = px(16)  # Adjust the font size as needed
   )%>% 
   cols_label(
         Cluster = md("**Cluster**"),
         n = md("**n**"),
         Richness = md("**Richness (mean)**"),
         Abundance = md("**Abundance (mean)**"),
          Phylum = md("**Phyla**"),
         Taxa = md("**Taxa**"))%>%
   cols_align(
          align = "center",
     columns = everything())%>%
     ## Right align Taxa column
     cols_align(
      # uni_tab2,
       align = "left",
       columns = 'Taxa'
     )%>%
 ## Add pie legend
 
#tab_caption(
    #caption = md("**Legend:** 
    tab_footnote(
    footnote = md("**Pie Chart Legend:**
<span style='color:#F8766D;font-size: 25px;'>■</span> Annelida (Ann), 
<span style='color:#DB8E00;font-size: 25px;'>■</span> Arthropoda (Art), 
<span style='color:#AEA200;font-size: 25px;'>■</span> Bryozoa (Bry), 
<span style='color:#64B200;font-size: 25px;'>■</span> Chordata (Cho), 
<span style='color:#00BD5C;font-size: 25px;'>■</span> Cnidaria (Cni), 
<span style='color:#00C1A7;font-size: 25px;'>■</span> Echinodermata (Ech), 
<span style='color:#00BADE;font-size: 25px;'>■</span> Entoprocta (Ent), 
<span style='color:#00A6FF;font-size: 25px;'>■</span> Nemertea (N),
<span style='color:#B385FF;font-size: 25px;'>■</span> Mollusca (Mol),
<span style='color:#EF67EB;font-size: 25px;'>■</span> Phoronida (Phor),
<span style='color:#FF63B6;font-size: 25px;'>■</span> Porifera (Por)"
))
 
 ## Save Response traits table
 uni_tab %>%
   gt::gtsave(
     "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/Table_S2.png",vwidth = 1650, vheight = 1600)#,  zoom = 2, expand = 10


 #_______________________________________________________________________________  
#### SIMPER ####
 
 library(vegan)
 
 # Load example data
 data(dune)
 data(dune.env)
 
 head(dune)
 head(dune.env)
 
class(dune)
class(dune.env)
 # Perform SIMPER analysis
 sim <- with(dune.env, simper(dune, Management, permutations = 99))
 
 # View the results
 summary(sim)
 
View(cluster_centers_4_simper)
 class(cluster_centers_4_simper)
 cluster_centers_4_simper.env <-data.frame(cluster=c("A2b",
                         "D2a",
                         "D1",
                         "B1b",
                         "D2c",
                         "C1b",
                         "D2b",
                         "A2a",
                         "A1",
                         "C1a",
                         "D2d"))
 class(cluster_centers_4_simper.env)
 
 
 # Perform SIMPER analysis
 sim <- with(cluster_centers_4_simper.env, simper(cluster_centers_4_simper, cluster, permutations = 99))
 summary(sim)
 
 # Summarize the results
 sim_summary <- summary(sim)
 class(sim_summary)
 sim
 # Extract species with high contributions
 characterizing_species <- sim_summary[[44]] %>%
   filter(cumsum < 0.05)  # Select species contributing up to 70% of the dissimilarity
 
 print(characterizing_species)
#_______________________________________________________________________________
 #### RANDOM FOREST (BIO): QUICK LOOK ####
 #_______________________________________________________________________________
 #### RANDOM FOREST (BIO): PREPARE DATA ####

## Start with cluster results df
names(faunal.cluster)

## Get columns in correct order
#FaunalCluster=faunal.cluster[,c(1,3,2,4)]#model based on numeric
FaunalCluster=faunal.cluster[,c(1,3,2,5)]#model based on character
#View(FaunalCluster)

## Change names of cols
colnames(FaunalCluster)=c("Sample","lon","lat","cluster")
head(FaunalCluster)

## Number of samples
dim(FaunalCluster)# 22814     4

## Take only the coordinates
FaunalCluster2 <- FaunalCluster[,2:3]
head(FaunalCluster2)
#_______________________________________________________________________________
#### RANDOM FOREST: CREATE RASTER STACK FOR ENV PREDICTOR VARIABLES ####

## Load packages
library(raster)
library(rgdal)

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

# Simple names for predictor variables
names(predictors)=c("Bathymetry","Current","Gravel","Light","Mud","Oxygen","Phyto","Salinity","Silicate","SPM","Temperature","WOV")#"Sand",

## Plot raster stack
plot(predictors)
#_______________________________________________________________________________
#### 20. RANDOM FOREST (BIO): EXTRACT PREDICTOR VARIABLES FROM RASTER STACK ####

head(FaunalCluster)
## Unload extract function from tidyr package (otherwise it won't work)
#.rs.unloadPackage("tidyr")

## Create a df for predictor variables
sdata <- raster::extract(predictors, FaunalCluster[,2:3])

#View(sdata)
class(sdata)

## Change from matrix to df
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(FaunalCluster$Sample,FaunalCluster$cluster,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "Cluster"
head(sdata2)

## Remove rows with NA (those falling outside raster extent)
sdata2 <- sdata2[complete.cases(sdata2), ]

## Change cols to appropriate type
str(sdata2)
sdata2$Cluster=as.factor(sdata2$Cluster)

## First check cols of correct type
str(sdata2)
dim(sdata2)# 22318    14
#_______________________________________________________________________________
#### RANDOM FOREST: MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Load packages
#install.packages("caTools")
library(caTools)

## Vector for response variable
Y <- sdata2[,2]

## Take a random 90% of samples in proportion across the 11 groups
set.seed(2)
msk= sample.split( Y, SplitRatio = 9/10, group = NULL )

## Check it was 90% 10% split

## The training set
train = sdata2[ msk,] 
dim(train)# 20087    14
#View(train)

## Remove station labels for train
train2 =train[,2:14]#was 11
#View(train2)

## The test set
test  = sdata2[ !msk,]
dim(test)# 2231   14
#View(test)

## Remove station labels for test
test2 =test[,2:14]#was 11
#View(test2)
str(test2)

## Check number of observations for train (TRUE) and test (FALSE) sets
print(table(Y, msk)) 

## Check number of samples in train and test sets is equal to total
dim(sdata2) # 22318    14
dim(train)+dim(test)# 22318    28
#_______________________________________________________________________________
#### RANDOM FOREST: DO MODELLING ####

## Call library
#install.packages("randomForest")
library(randomForest)

## Model
model <- factor(Cluster) ~Bathymetry+Current+Gravel+Light+Mud+Oxygen+Phyto+Salinity+Silicate+SPM+Temperature+WOV

## Run model
rf2 <- randomForest(model, data=train2,na.action=na.exclude)
#_______________________________________________________________________________
#### RANDOM FOREST: EXAMINE HOW VARIABLES ARE INFLUENCING THE MODEL ####

## Produce plots 
varImpPlot(rf2)

png('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/variables_affecting_model_structure.png') #, height=nrow(pr), width=ncol(pr) EFFECTS TRAITS
varImpPlot(rf2)
dev.off()

preds <- names(rf2$forest$xlevels)

for (i in 1:length(preds)){
  partialPlot(rf2, train2, preds[i], which.class ='1')
  next
}
#_______________________________________________________________________________
#### RANDOM FOREST: EVALUATE THE MODEL PERFORMANCE ####

## Predict cluster group for test set
pred <- predict(rf2, newdata = test2)
table(pred, test2$Cluster)

## We can test the accuracy as follows:
(22+21+37+93+84+52+88+104+188+464+203)/ nrow(test2)#60.8

## Confusion matrix plot
#https://stackoverflow.com/questions/37897252/plot-confusion-matrix-in-r-using-ggplot

confusion_matrix <- as.data.frame(table(pred, test2$Cluster))

cm <- ggplot(confusion_matrix, aes(pred,sort(Var2,decreasing = T), fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#737373") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11")) +
  scale_y_discrete(labels=c("11","10","9","8","7","6","5","4","3","2","1"))

## Save confusion matrix
png('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/structure_confusion_matrix_plot.png')
#png('OUTPUTS/EFFECTS_TRAITS/et__confusion_matrix_plot.png')
cm
dev.off()
#_______________________________________________________________________________
#### RANDOM FOREST: PRODUCE FULL COVERAGE RASTER ####

## Use model to predict cluster group for each raster cell
pr <- predict(predictors, rf2)
plot(pr)
#_______________________________________________________________________________
#### 26. RANDOM FOREST (BIO): OUTPUT RASTER AS TIFF ####

## Save raster
writeRaster(pr,'C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/StructureCluster.tif',overwrite=TRUE,format = "GTiff")

#_______________________________________________________________________________
#### SIMPLE RASTER PLOT (NO BATHY) ####

## Colours in correct order (see results of cluster where nos are matched to groups)

colours <- c(
  "#0000ee",#9
  "#00ffff",#8
  "#05aac1",#1
  
  "#9a32cd",#4
  "#00cd00",#10
  "#9aff9a",#6
  "#b40202",#3
  "#ff0000",#2
   "#ff8c00",#7
  
  "#ffff00",#5
  
 
 
  
  
  "#b4b404")#11
  
## Plot raster minus legend
plot(pr,col=colours,legend = FALSE)

## Add category legend
legend("bottomright", legend = c("A1","A2a","A2b","B1b","C1a","C1b","D1","D2a","D2b","D2c","D2d"), fill = c("#0000ee","#00ffff","#05aac1","#9a32cd","#00cd00","#9aff9a","#b40202","#ff0000","#ff8c00","#ffff00","#b4b404"))
#_______________________________________________________________________________
#### FAUNAL CLUSTER MAP: GGPLOT ####

## Load packages
library(raster)
library(ggplot2)
library(scales)


## Load raster
pr = raster('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/StructureCluster.tif')

## Reduce size of raster using aggregate function (e.g. ggregate from 40x40 resolution to 120x120 (factor = 3))
pr.agg <- aggregate(pr, fact=4,fun = modal)
res(pr.agg)

## Plot raster
plot(pr.agg)

## Convert the reduced raster to points for plotting
s_raster.p <- rasterToPoints(pr.agg)

## Make the points a dataframe for ggplot
sdf <- data.frame(s_raster.p)

## Make appropriate column headings
colnames(sdf) <- c("Longitude", "Latitude", "FCluster")
#View(fdf)

## Check there are only 11 cluster groups
unique(sdf$FCluster)

## Change cluster group from numeric to character
sdf$FCluster <- as.character(sdf$FCluster)
str(sdf)

## Change numbers to codes
sdf$FCluster[sdf$FCluster == "1"] <- "A2b"
sdf$FCluster[sdf$FCluster == "2"] <- "D2a"
sdf$FCluster[sdf$FCluster == "3"] <- "D1"
sdf$FCluster[sdf$FCluster == "4"] <- "B1b"
sdf$FCluster[sdf$FCluster == "5"] <- "D2c"
sdf$FCluster[sdf$FCluster == "6"] <- "C1b"
sdf$FCluster[sdf$FCluster == "7"] <- "D2b"
sdf$FCluster[sdf$FCluster == "8"] <- "A2a"
sdf$FCluster[sdf$FCluster == "9"] <- "A1"
sdf$FCluster[sdf$FCluster == "10"] <- "C1a"
sdf$FCluster[sdf$FCluster == "11"] <- "D2d"

## Make Cluster a factor
sdf$FCluster=as.factor(sdf$FCluster)

## Now make the map
modfc=ggplot(data=sdf, aes(y=Latitude, x=Longitude)) +
  geom_raster(aes(fill=FCluster)) +
  scale_fill_manual(values = c(   "#0000ee",
                                  "#00ffff",
                                  "#05aac1",
                                  "#9a32cd",
                                  "#00cd00",
                                  "#9aff9a",
                                  "#b40202",
                                  "#ff0000",
                                  "#ff8c00",
                                  "#ffff00",
                                  "#b4b404"),name="Cluster")+
  labs(x="Longitude",y="Latitude")+
  theme_bw(base_size=12)+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  coord_quickmap(xlim = c(-10.5, 10.7),ylim=c(48.7, 60.9))+
  theme(legend.text=element_text(size=18))

modfc
#_______________________________________________________________________________
#### 2D RASTER CLUSTER PLOT WITH UNDERLYING BATHY (AS HILLSHADE) ####

## Load packages
library(raster)

## Load raster data
bathy <- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/bathy3.tif")
pr = raster('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/StructureCluster.tif')

## Exagerate the vertical scale (bathy not very clear unless you do this)
bathy2 <- bathy*bathy

## Generate 'slope' and 'aspect' using the terrain function
slope <- terrain(bathy2,opt='slope')
aspect <- terrain(bathy2,opt='aspect')

## Generate hillshade
hill <- hillShade(slope,aspect,20,0)

# Define a colour palatte (see https://mycolor.space/?hex=%23C2BF5E&sub=1)
colours <- c( "#05aac1",#1
              "#ff0000",#2
              "#b40202",#3
              "#9a32cd",#4
              "#ffff00",#5
              "#9aff9a",#6
              "#ff8c00",#7
              "#00ffff",#8
              "#0000ee",#9
              "#00cd00",#10
              "#b4b404")

## Produce plot
png('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/raster_structure_2Dhillshadebathy.png') 
plot(hill,col=grey(1:100/100),legend=FALSE)#,axes=F,box=FALSE
plot(pr, col=colours,add=TRUE,alpha=0.65,axes=F,box=FALSE,legend=FALSE)#
legend(x = 8.3, y = 53.9, legend = c("A1","A2a","A2b","B1b","C1a","C1b","D1","D2a","D2b","D2c","D2d"), fill = c("#0000ee","#00ffff","#05aac1","#9a32cd","#00cd00","#9aff9a","#b40202","#ff0000","#ff8c00","#ffff00","#b4b404"),cex = 0.85, inset = 0.9,bty = "n") # bty used to turn off legend border)
dev.off()
#_______________________________________________________________________________
#### FINAL 3D PLOT ####
#https://nils.ratnaweera.net/2020/06/06/using-rayshader-to-visualize-lake-zurich/

## Load packages
library(raster)
library(sf)
library(tidyverse)
library(rayshader)
#library(lazyraster)

## Bring in bathy data
bathy <- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/bathy3.tif")
bathy.agg <- aggregate(bathy, fact=8,fun = modal) # reduce the size of the raster
#plot(bathy.agg)

## Bring in Cluster data
pr = raster('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/StructureCluster.tif')
#pr = raster('Y:/C8210_OneBenthic/Working_Area/SDM/OUTPUT/AssemblageMaxClass.tif') #Anna's version
pr.agg <- aggregate(pr, fact=8,fun = modal)
#plot(pr.agg)

## Convert rater to matrix for use with rayshader 
elmat = raster_to_matrix(bathy.agg)

## define colours you need
colours <- c("#05aac1",#1
             "#ff0000",#2
             "#b40202",#3
             "#9a32cd",#4
             "#ffff00",#5
             "#9aff9a",#6
             "#ff8c00",#7
             "#00ffff",#8
             "#0000ee",#9
             "#00cd00",#10
             "#b4b404")
#colours <- c('#0000ee',  '#00ffff',  '#05aac1',  '#9a32cd',  '#00cd00',  '#9aff9a',  '#b40202',  '#ff0000',  '#ff8c00',  '#ffff00',  '#b4b404')#Annas version

## convert overlay to rgb file
x <- RGB(pr.agg,col=colours)
plot(x)

## Create objects for red, green and blue
cluster_r = x[[1]]
cluster_g = x[[2]]
cluster_b = x[[3]]
plot(cluster_r)

## Raster stack 
cluster_rbg = raster::stack(cluster_r, cluster_g, cluster_b)
raster::plotRGB(cluster_rbg, scale=255^2)

#cluster_rbg_corrected = sqrt(raster::stack(cluster_r, cluster_g, cluster_b))
#raster::plotRGB(cluster_rbg_corrected)

## convert to matrix for use with rayshader
cluster_r = rayshader::raster_to_matrix(cluster_rbg$red)
cluster_g = rayshader::raster_to_matrix(cluster_rbg$green)
cluster_b = rayshader::raster_to_matrix(cluster_rbg$blue)

## Create an create an 3-layer RGB array
cluster_rgb_array = array(0,dim=c(nrow(cluster_r),ncol(cluster_r),3))

## Populate array
cluster_rgb_array[,,1] = cluster_r/255 #Red layer
cluster_rgb_array[,,2] = cluster_g/255 #Blue layer
cluster_rgb_array[,,3] = cluster_b/255 #Green layer

## multi-dimensional transpose as rasters and arrays are orientated differently
cluster_rgb_array = aperm(cluster_rgb_array, c(2,1,3))

plot_map(cluster_rgb_array)

z_scale = 6

#----------------
### DEM PLOT ###
# Azimuth (phi) is the viewing angle (0 horizon, 90 looking down) 
# theta is the angle of approachangle (180 = looking north)
hill <- elmat %>%
  #height_shade(texture = rev(rainbow(256)))
  sphere_shade() 
dim(hill)
dim(cluster_rgb_array)
rgl::clear3d()

plot_3d(hillshade = hill,heightmap = elmat, windowsize = c(1000,600),zscale=z_scale,theta = 0,zoom=0.5,phi=35,solid=FALSE)

#----------------
## ADD IMAGE OVERLAY ###
hill <- add_overlay(hill, cluster_rgb_array)

rgl::clear3d()

plot_3d(hillshade = hill, heightmap = elmat, windowsize = c(1000, 600), zscale = z_scale, theta = 0, zoom = 0.5, phi = 35,solid = FALSE,shadowcolor = "grey80")
#----------------
## ADD SHADE ###
hill <- add_shadow(hillshade = hill, shadowmap = ambient_shade(elmat),max_darken = 0.1)

rgl::clear3d()

plot_3d(hillshade = hill, heightmap = elmat, windowsize = c(1000, 600), zscale = z_scale, theta = 0, zoom = 0.5, phi = 40,solid = FALSE,shadow = TRUE,background = "white",shadowcolor = "grey80")
#-----------------
## add water ##

#rgl::clear3d()

#plot_3d(hillshade = hill, heightmap = elmat, windowsize = c(1000, 600), zscale = z_scale, theta = 0, zoom = 0.5, phi = 40,solid = TRUE,shadow = TRUE,background = "white",shadowcolor = "grey80",water=TRUE, wateralpha = 0.7, watercolor = "lightblue", waterdepth = -50, waterlinecolor = "lightblue", waterlinealpha = 0.5)

#-----------------
## Capture images of map at different angles

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/')
#setwd('C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/OneBenthicBiodiversityLayersApp/SOCIAL MEDIA/')
## Generate png files
angles= seq(0,360,length.out = 300)[-1]
for(i in 1:300) {
  render_camera(theta=-45+angles[i])
  render_snapshot(filename = sprintf("cluster%i.png", i))
}
rgl::rgl.close()

## make a movie using https://gifmaker.me/ Animation speed: 150 milliseconds, canvas size 150%
#_______________________________________________________________________________
#### Try making video with R ####
library(av)

# Set the directory containing your PNG files
img_dir <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/"

# List all PNG files in the directory
png_files <- list.files(img_dir, pattern = "*.png", full.names = TRUE)

# Create the video
av_encode_video(png_files, output = "output_video.mp4", framerate = NULL)
#_______________________________________________________________________________
#### DRAIN THE SEA #####
# You'll need to run above code (to generate overlay) before you run this. If overlay not required then it's ok
#https://wcmbishop.github.io/rayshader-demo/

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicCertain/OUTPUTS/')

## Load packages
library(raster)
library(magrittr)
library(stringr)
library(rgl)
library(rayshader)
library(magick)
library(ggplot2)
library(rgdal)
library(webshot)

## Load raster data
bathy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/bathy3.tif")
crs(bathy) <- "+proj=longlat +datum=WGS84 +no_defs" # Update crs

## Produce simplified version of raster layers
bathy.agg <- aggregate(bathy, fact=6,fun = modal)

## Convert bathy to matrix
elmat = raster_to_matrix(bathy.agg)

## Create png files for different water levels
n_frames <- 82#21
zscale <- 6

# frame transition variables
#waterdepthvalues <- min(elev_matrix)/2 - min(elev_matrix)/2 * cos(seq(0,2*pi,length.out = n_frames))
#waterdepthvalues <- c(0,-62.5)
waterdepthvalues <- c(0,-2.5,-5,-7.5,-10,-12.5,-15,-17.5,-20,-22.5,-25,-27.5,-30,-32.5,-35,-37.5,-40,-42.5,-45,-47.5,-50,-52.5,-55,-57.5,-60,-62.5,-65,-67.5,-70,-72.5,-75,-77.5,-80,-82.5,-85,-87.5,-90,-92.5,-95,-97.5,-100,-105,-110,-115,-120,-125,-130,-135,-140,-145,-150,-155,-160,-165,-170,-175,-180,-185,-190,-195,-200,-205,-210,-215,-220,-225,-230,-235,-240,-245,-250,-260,-270,-280,-290,-300,-310,-320,-330,-340,-350,-360)
#thetavalues <- -90 + 45 * cos(seq(0, 2*pi, length.out = n_frames))
thetavalues <- 0 + 0 * cos(seq(0, 2*pi, length.out = n_frames))
# shadow layers
ambmat <- ambient_shade(elmat, zscale = zscale)
raymat <- ray_shade(elmat, zscale = zscale, lambert = TRUE)

## Generate .png frame images (wateralpha=0.8)
img_frames <- paste0("drain_structure", seq_len(n_frames), ".png")
for (i in seq_len(n_frames)) {
  message(paste(" - image", i, "of", n_frames))
  elmat %>%
    sphere_shade(texture = "desert") %>%
    add_overlay(cluster_rgb_array)%>%# REMOVE IF OVERLAY NOT REQUIRED
    add_shadow(ambmat, 0.5) %>%
    add_shadow(raymat, 0.5) %>%
    plot_3d(elmat, solid = TRUE, shadow = TRUE, zscale = zscale, 
            water = TRUE, watercolor = "imhof3", wateralpha = 0.94, #watercolor="imhof3" "dodgerblue"
            waterlinecolor = "#ffffff", waterlinealpha = 0.5,
            #waterdepth = waterdepthvalues[i]/zscale, 
            waterdepth = waterdepthvalues[i], 
            theta = thetavalues[i], phi = 45)
  render_snapshot(img_frames[i])
  rgl::clear3d()
}

## make a movie using https://gifmaker.me/

#____________________________________________________________________________________________________________________



