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
  
  # Load and process tag metadata
  tags.meta <- as.data.table(
    tbl(sqldb, "tagDeps") %>%
      select(deployID, tagID, projectID, tsStart, tsEnd, deferSec, speciesID, 
             markerNumber, latitude, longitude, fullID, comments) %>%
      filter(projectID == 213, speciesID %in% sp.list) %>%
      collect()
  )
  tags.meta <- tags.meta[
    !is.na(longitude),
    .(tagID, tsStart = as_datetime(tsStart, tz = "UTC", origin = "1970-01-01"),
      tsEnd = as_datetime(tsEnd, tz = "UTC", origin = "1970-01-01"), speciesID)
  ][, .SD[!duplicated(tagID)], by = tagID]
  
  # Ambiguous tags
  tags.ambigs <- as.data.table(
    tbl(sqldb, "allambigs") %>% collect()
  )[motusTagID %in% tags.meta$tagID, .(ambigID, tagID = motusTagID)]
  
  # Detection extraction
  tags.detect <- as.data.table(
    tbl(sqldb, "alltags") %>%
      filter(tagProjID == 213, motusTagID %in% c(tags.meta$tagID, tags.ambigs$tagID), 
             speciesID %in% sp.list) %>%
      select(runID, ts, sig, port, noise, freqsd, motusTagID, 
             ambigID, runLen, tagDeployID, tagDeployStart, 
             tagDeployEnd, tagDepLat, tagDepLon, deviceID, recvDeployID, recv,
             speciesID, mfgID, motusFilter, antBearing, antType, 
             antHeight) %>%
      collect()
  )
  tags.detect <- tags.detect[
    , `:=`(ts = as_datetime(ts), 
           tagDeployStart = as_datetime(tagDeployStart), 
           tagDeployEnd = as_datetime(tagDeployEnd))
  ]
  
  # Receiver metadata
  recv.meta <- as.data.table(
    tbl(sqldb, "recvDeps") %>% collect()
  )[!is.na(longitude), .(recvDeployID = deployID, recvProjectID = projectID, 
                         tsStart = as_datetime(tsStart), tsEnd = as_datetime(tsEnd),
                         longitude, latitude)]
  
  # Join receiver metadata
  recvs <- tags.detect[
    !is.na(recvDeployID),
    .(recv, recvDeployID)
  ][unique(recv.meta), on = "recvDeployID"]
  
  tags.detect <- tags.detect[
    recvs,
    on = .(recv, recvDeployID),
    `:=`(recvLat = i.latitude, recvLon = i.longitude, recvType = i.recvType)
  ][ts >= tagDeployStart & ts <= tagDeployEnd & !is.na(recvLat)]
  
  # Filter detections
  tags.detect <- tags.detect[
    freqsd <= 0.1, motusFilter := fifelse(motusFilter > 0, motusFilter, 0)
  ][, af := sum(motusFilter), by = motusTagID][af > 0]
  
  # Secondary filtering
  tags.detect <- tags.detect[
    , `:=`(n = .N, ntrue = sum(motusFilter), ptrue = ntrue / .N), 
    by = .(date(ts), recv, motusTagID)
  ][ptrue > 0.30]
  
  tags.filtered <- tags.detect[
    , `:=`(ftemp = 1, motusFilter = 1)
  ]
  
  # Transit check
  tags.filtered <- add_deploy_loc(tags.filtered)
  transit.check <- site_transit_min(tags.filtered)
  transit.suspect <- transit.check[suspect.transit == "suspect"]
  
  # Connection states
  tc <- transit.check[
    , state := fifelse(dist.min == 0 | rate >= 5, "connected", "not_connected")
  ][state == "connected", .(motusTagID, ts.x, lat.x, lon.x, lat.y, lon.y)]
  
  # Outputs
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat("Elapsed Time: ", format(round(elapsed_time, 2), nsmall = 2), " minutes\n")
  
  return(list(
    tags.detect = tags.detect,
    tags.filtered = tags.filtered,
    tags.meta = tags.meta,
    transit.suspect = transit.suspect,
    connections = tc
  ))
}

#Run function
Result <- Detection_Extraction_and_Filter(sqldb, sp.list)

# Save outputs to CSV
fwrite(Result$tags.detect, "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_detections_UNFILTERED_12.31.24.csv")
fwrite(Result$tags.filtered, "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_detections_FILTERED_12.31.24.csv")
fwrite(Result$tags.meta, "~/Project 213 Data/Processed Detection Data/CONI/CONI_tag_dep_metadata_12.31.24.csv")


#####

