##Function for SQL extraction##
#Updated 6.17.24#

library(motus)
library(dplyr)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggmap)
library(DBI)
library(RSQLite)
library(sf)
library(viridis)
library(PROJ)
library(geos) 
Sys.setenv(TZ = "UTC")
source("~/Motus/Scripts/useful_functions.R") 


sqldb<-tagme(projRecv = 213, new = FALSE, update = FALSE, dir = "~/Project 213 Data/SQL") 
sp.list<-c(7720) 
##COPO = 7750, CONI = 7720, NSWO = 7680, LEOW = 7620, FLOW = 7290, LEWO = 10030, NOFL = 10460

Detection_Extraction<- function(sqldb, sp.list) {
  start_time <- Sys.time()
  #tag metadata
  tagmet<-tbl(sqldb, "tagDeps")
  tags.meta<- tagmet %>%
    select(deployID, tagID, projectID, tsStart, tsEnd, deferSec, speciesID, markerNumber, latitude, longitude, fullID, comments) %>%
    filter(projectID == 213 ) %>%
    filter(speciesID %in% sp.list) %>%
    collect() %>%
    as.data.frame() %>%
    mutate(tsStart = as_datetime(tsStart, tz = "UTC", origin = "1970-01-01"),
           tsEnd = as_datetime(tsEnd, tz = "UTC", origin = "1970-01-01")) %>%
    distinct() %>%
    filter(!is.na(longitude)) %>%
    group_by(tagID) %>%
    mutate(n = n()) %>%
    ungroup()
  #ambiguous tags
  allambigs<-tbl(sqldb, "allambigs")
  ambig<-allambigs %>%
    collect() %>%
    as.data.frame()
  tags.ambigs<-ambig %>%
    group_by(ambigID) %>%
    mutate(inds = ifelse(any(motusTagID %in% tags.meta$tagID), "yes", "no")) %>%
    filter(inds == "yes") %>%
    ungroup() %>%
    rename(tagID = motusTagID)
  #tag list
  tag.list<-tags.meta %>%
    select(tagID) %>%
    bind_rows(., tags.ambigs %>% select(tagID)) %>%
    distinct() %>% pull()
  #detection extraction
  alltag<- tbl(sqldb, "alltags")
  tags.detect<-alltag %>%
    filter(tagProjID == 213 & motusTagID %in% tag.list & speciesID %in% sp.list) %>%
    select(runID, ts, sig, port, noise, freqsd, motusTagID, 
           ambigID, runLen, tagProjID, tagDeployID, tagDeployStart, tagDeployEnd, 
           tagDepLat, tagDepLon, deviceID, recvDeployID, recv,
           speciesID, markerNumber, mfgID, motusFilter, antBearing, antType, antHeight) %>%
    collect() %>%
    as.data.frame() %>%
    mutate(ts = as_datetime(ts),
           tagDeployStart = as_datetime(tagDeployStart),
           tagDeployEnd = as_datetime(tagDeployEnd))
  #receiver metadata
  recvmet<-tbl(sqldb, "recvDeps")
  recv.meta<-recvmet %>%
    collect() %>%
    as.data.frame() %>%
    mutate(tsStart = as_datetime(tsStart),
           tsEnd = as_datetime(tsEnd))
  recv.meta<-recv.meta %>%
    rename(recvDeployID = deployID, recvProjectID = projectID) %>%
    distinct() %>%
    filter(!is.na(longitude)) %>%
    group_by(deviceID) %>%
    mutate(n = n()) %>%
    ungroup()
  #join tables
  recvs<-tags.detect %>%
    select(receiverID = recv, recvDeployID) %>%
    filter(!is.na(recvDeployID)) %>% 
    distinct() %>%
    left_join(., recv.meta) %>%
    select (receiverID, recvDeployID, recvProjectID, stationName, 
            latitude, longitude, receiverType, tsStart, tsEnd, isMobile) %>%
    filter(!is.na(longitude))
  tags.detect<-left_join(tags.detect, recvs %>%
                           select(recv = receiverID, recvDeployID, receiverType, 
                                  recvName = stationName, recvDeployLat = latitude, recvDeployLon = longitude,
                                  isMobile, recvProjID = recvProjectID) %>%
                           distinct())
  #preliminary filters
  tags.detect<-tags.detect %>%
    mutate(date = as.Date(ts), 
           year = as.numeric(format(ts, '%Y')),
           month = as.numeric(format(ts, '%m'))) %>%
    filter(ts>=tagDeployStart & ts <= tagDeployEnd,
           !is.na(recvDeployLat))
  tags.detect<-tags.detect %>%
    mutate(motusFilter = ifelse(freqsd > 0.1, 0, motusFilter))
  tags.detect<-tags.detect %>%
    group_by(motusTagID) %>%
    mutate(af = sum(motusFilter)) %>%
    filter(af > 0) %>%
    ungroup() %>%
    select(-af)
  tag.list.fil<-tags.detect %>%
    select(motusTagID) %>%
    distinct() %>% pull()
  tags.meta<-tags.meta %>%
    filter(tagID %in% tag.list)
  
  return(list(tags.detect = tags.detect, tags.meta = tags.meta, tag.list.fil = tag.list.fil, tags.ambig = tags.ambigs))
  end_time <- Sys.time()
  
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("Elapsed Time: ", format(round(elapsed_time, 2), nsmall = 2), " minutes\n")
}


Result<- Detection_Extraction(sqldb, sp.list)

Detections <- Result$tags.detect
Tag_metadata<- Result$tags.meta
Filtered_tag_list <- Result$tag.list.fil
Ambiguous_tags <- Result$tags.ambig


write.csv(Detections, file = "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_detections_UNFILTERED_12.28.24.csv", row.names = FALSE)
write.csv(Tag_metadata, file = "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_dep_metadata_12.28.24.csv", row.names = FALSE)

deployments<-Detections %>% select(motusTagID, tagDeployID) %>%
  distinct() %>%
  group_by(motusTagID) %>%
  mutate(n = n()) %>%
  ungroup()
any(deployments$n>1)
#any tags with more than one deployment, you will need to use tagDeployID instead of motusTagID



Secondary_Filter<- function(Detections) {
  tags.detpro<-Detections
  
  tags.detpro<-tags.detpro %>%
    group_by(date, recvName, motusTagID) %>%
    mutate(n = n(),
           ntrue = sum(motusFilter),
           ptrue = ntrue/n) %>%
    ungroup()
  
  tags.filt<-tags.detpro %>%
    mutate(ftemp = ifelse(ptrue > 0.30, 1, 0),
           motusFilter = ftemp) %>%
    filter(ftemp == 1)
  
  tags.filt.list <- unique(tags.filt$motusTagID)
  
  transit.check<- tags.filt %>%
    do(add_deploy_loc(.)) %>%
    do(site_transit_min(.))
  
  transit.check.suspect<-transit.check %>%
    filter(suspect.transit == "suspect")
  
  tc<-transit.check %>%
    mutate(state = ifelse(dist.min == 0 | rate >= 5, "connected", "not_connected")) %>%
    filter(state == "connected") %>%
    select(motusTagID, ts.x, lat.x, lon.x, lat.y, lon.y)
  
  return(list(tags.filt.list = tags.filt.list, transit.check.suspect = transit.check.suspect, tc = tc,
              tags.filt = tags.filt))
  
}

Filtered_Result<- Secondary_Filter(Detections)

Filtered_Detections <- Filtered_Result$tags.filt
Filtered_tag_list<- Filtered_Result$tags.filt.list
Transit_check<- Filtered_Result$transit.check.suspect
Connections<- Filtered_Result$tc

write.csv(Filtered_Detections, file = "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_detections_FILTERED_12.28.24.csv", row.names = FALSE)


##Plots##

#setting the theme for mapping with ggmap
theme_set(theme_bw())
world<-ne_countries(scale = "medium", returnclass = "sf")

#to download a shapefile for the first time use this code
#this directory is set to download into a folder named Map Data in your wd
#you can change or modify this directory as needed
lakes<-ne_download(scale = "medium", type = 'lakes', category = 'physical',
                   returnclass = "sf", destdir = "./Map Data/lakes")

#to load an already downloaded shapefile use this code
lakes<-ne_load(type = "lakes", scale = "medium" , category = 'physical',
               returnclass = "sf",
               destdir = paste0("~/Map Data/lakes"))


tags.path.filt<-fun.getpath(Filtered_Detections)
unique_tag_list <- unique(Filtered_tag_list)

for (current_tag in unique_tag_list) {
  bird <- tags.path.filt %>% filter(motusTagID == current_tag)
  pbird <- Connections %>% filter(motusTagID == current_tag) %>%
    distinct()
  
  xmin <- min(bird$recvDeployLon, bird$tagDepLon, na.rm = TRUE) - 2
  xmax <- max(bird$recvDeployLon, bird$tagDepLon, na.rm = TRUE) + 2
  ymin <- min(bird$recvDeployLat, bird$tagDepLat, na.rm = TRUE) - 1
  ymax <- max(bird$recvDeployLat, bird$tagDepLat, na.rm = TRUE) + 1
  
  filename <- paste("./Processed Detection Data/CONI/Diagnostic Figures/", current_tag, "-map.png", sep = "")
  
  png(filename = filename,
      width = 8, height = 8, units = "in", res = 600)
  
  print(ggplot(data = world) +
          geom_sf() +
          geom_sf(data = lakes, fill = "white") +
          coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
          geom_point(data = bird, aes(recvDeployLon, recvDeployLat, colour = month), alpha = 0.5, size = 5) +
          geom_segment(data = pbird, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y),
                       linetype = "dashed",
                       arrow = arrow(angle = 18, type = "closed", length = unit(0.1, "inches"))) +
          geom_point(data = bird,
                     aes(tagDepLon, tagDepLat), colour = "red",
                     shape = 8, size = 4) +
          scale_color_gradientn(colors = rev(viridis_pal()(12)), limits = c(1, 12)) +
          ggtitle(paste0("motusTagID:", current_tag)))
  
  dev.off()
}


##
for (i in 1:length(Filtered_tag_list)){
  png(filename = paste0("./Processed Detection Data/CONI/Diagnostic Figures/Signal Strength Plots", "sigstrength-",
                        Filtered_tag_list[i], ".png"),
      width = 11, height = 7.5, units = "in", res=600)
  print(plotTagSig_mod(Filtered_Detections, motusTagID = Filtered_tag_list[i]))
  dev.off()
}
rm(i)
