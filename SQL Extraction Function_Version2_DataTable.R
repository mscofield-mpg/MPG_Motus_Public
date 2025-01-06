###SQL Extraction with data.table###
#12.29.24#
library(motus)
library(RSQLite)
library(data.table)
Sys.setenv(TZ = "UTC")
source("./Scripts/useful_functions.R") 



sqldb<-tagme(projRecv = 213, new = FALSE, update = FALSE, dir = "~/Project 213 Data/SQL") 
sp.list<-c(7720) 
##COPO = 7750, CONI = 7720, NSWO = 7680, LEOW = 7620, FLOW = 7290, LEWO = 10030, NOFL = 10460



Detection_Extraction_and_Filter <- function(sqldb, sp.list) {
  start_time <- Sys.time()
  
  # Tag metadata
  tagmet <- tbl(sqldb, "tagDeps")
  tags.meta <- tagmet %>%
    select(deployID, tagID, projectID, tsStart, tsEnd, deferSec, speciesID, markerNumber, latitude, longitude, fullID, comments) %>%
    filter(projectID == 213) %>%
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
  
  # Ambiguous tags
  allambigs <- tbl(sqldb, "allambigs")
  ambig <- allambigs %>%
    collect() %>%
    as.data.frame()
  
  tags.ambigs <- ambig %>%
    group_by(ambigID) %>%
    mutate(inds = ifelse(any(motusTagID %in% tags.meta$tagID), "yes", "no")) %>%
    filter(inds == "yes") %>%
    ungroup() %>%
    rename(tagID = motusTagID)
  
  # Tag list
  tag.list <- tags.meta %>%
    select(tagID) %>%
    bind_rows(., tags.ambigs %>% select(tagID)) %>%
    distinct() %>%
    pull()
  
  # Detection extraction
  alltag <- tbl(sqldb, "alltags")
  tags.detect <- alltag %>%
    filter(tagProjID == 213 & motusTagID %in% tag.list & speciesID %in% sp.list) %>%
    select(ts, sig, port, noise, freqsd, motusTagID, ambigID, runLen, tagProjID, 
           tagDeployID, tagDeployStart, tagDeployEnd, tagDepLat, tagDepLon, deviceID, 
           recvDeployID, recv, speciesID, markerNumber, mfgID, motusFilter, 
           antBearing, antType, antHeight) %>%  # Removed hitID, batchID, runID
    collect() %>%
    as.data.frame() %>%
    mutate(ts = as_datetime(ts),
           tagDeployStart = as_datetime(tagDeployStart),
           tagDeployEnd = as_datetime(tagDeployEnd)) %>%
    distinct()
  
  # Receiver metadata
  recvmet <- tbl(sqldb, "recvDeps")
  recv.meta <- recvmet %>%
    collect() %>%
    as.data.frame() %>%
    mutate(tsStart = as_datetime(tsStart),
           tsEnd = as_datetime(tsEnd))
  
  recv.meta <- recv.meta %>%
    rename(recvDeployID = deployID, recvProjectID = projectID) %>%
    distinct() %>%
    filter(!is.na(longitude)) %>%
    group_by(deviceID) %>%
    mutate(n = n()) %>%
    ungroup()
  
  # Join tables
  recvs <- tags.detect %>%
    select(receiverID = recv, recvDeployID) %>%
    filter(!is.na(recvDeployID)) %>% 
    distinct() %>%
    left_join(., recv.meta) %>%
    select(receiverID, recvDeployID, recvProjectID, stationName, 
           latitude, longitude, receiverType, tsStart, tsEnd, isMobile) %>%
    filter(!is.na(longitude))
  
  tags.detect <- left_join(tags.detect, recvs %>%
                             select(recv = receiverID, recvDeployID, receiverType, 
                                    recvName = stationName, recvDeployLat = latitude, recvDeployLon = longitude,
                                    isMobile, recvProjID = recvProjectID) %>%
                             distinct())
  
  # Preliminary filters
  tags.detect <- tags.detect %>%
    mutate(date = as.Date(ts), 
           year = as.numeric(format(ts, '%Y')),
           month = as.numeric(format(ts, '%m')),
           motusFilter = ifelse(freqsd > 0.1, 0, motusFilter)) %>%
    filter(ts >= tagDeployStart & ts <= tagDeployEnd,
           !is.na(recvDeployLat)) %>%
    group_by(motusTagID) %>%
    mutate(af = sum(motusFilter)) %>%
    filter(af > 0) %>%
    ungroup() %>%
    select(-af)
  
  # Filter detections
  tags.detpro <- tags.detect %>%
    group_by(date, recvName, motusTagID) %>%
    mutate(n = n(),
           ntrue = sum(motusFilter),
           ptrue = ntrue / n) %>%
    ungroup()
  
  tags.filt <- tags.detpro %>%
    mutate(ftemp = ifelse(ptrue > 0.30, 1, 0),
           motusFilter = ftemp) %>%
    filter(ftemp == 1)
  
  tags.filt.list <- unique(tags.filt$motusTagID)
  
  transit.check <- tags.filt %>%
    do(add_deploy_loc(.)) %>%
    do(site_transit_min(.))
  
  transit.check.suspect <- transit.check %>%
    filter(suspect.transit == "suspect")
  
  tc <- transit.check %>%
    mutate(state = ifelse(dist.min == 0 | rate >= 5, "connected", "not_connected")) %>%
    filter(state == "connected") %>%
    select(motusTagID, ts.x, lat.x, lon.x, lat.y, lon.y)
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("Elapsed Time: ", format(round(elapsed_time, 2), nsmall = 2), " minutes\n")
  
  return(list(tags.detect = tags.detect, 
              tags.meta = tags.meta, 
              tags.filtered = tags.filt,
              tags.filt.list = tags.filt.list, 
              transit.check.suspect = transit.check.suspect, 
              connections = tc))
}

# Using the function
Result <- Detection_Extraction_and_Filter(sqldb, sp.list)

# Save outputs to CSV
fwrite(Result$tags.detect, "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_detections_UNFILTERED_12.31.24.csv")
fwrite(Result$tags.filtered, "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_detections_FILTERED_12.31.24.csv")
fwrite(Result$tags.meta, "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_dep_metadata_12.31.24.csv")


#####

