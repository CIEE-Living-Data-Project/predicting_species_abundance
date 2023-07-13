#Title: Functions to calculate portion change for species pairs 
#Author: Alexandre Fuster
#Date: June 20th 2023






# Get portion change for abundance variables



get.log.prop.change_abund <-function(x,data,pairs.keep){
  
  #get pair
  genera.pair=pairs.keep[x,]
  
  #get studies
  study_1 <- data[which(data$ID==str_split(genera.pair$PAIR.ID,"_")[[1]][1]),] #get all timeseries for ID1
  study_2 <- data[which(data$ID==str_split(genera.pair$PAIR.ID,"_")[[1]][2]),] #get all timeseries for ID2
  
  #get years
  years1=study_1$YEAR[which(study_1$GENUS==genera.pair$Gn1)]
  years2=study_2$YEAR[which(study_2$GENUS==genera.pair$Gn2)]
  
  cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
  
  dat<-data.frame("Period"=NA,"value"=NA)
  for (y in 1:length(cont.overlapping)){
    dat[y,1]=y
    dat[y,2]=length(cont.overlapping[[y]])}
  
  years.overlap=cont.overlapping[[dat$Period[which(dat$value==max(dat$value))]]]
  
  n=dat$Period[which(dat$value>9)]
  
  if (length(n)>1){
    res=data.frame(matrix(ncol=21,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2", 
                                                                "log_prop_change_gn1_mean","log_prop_change_gn2_mean",
                                                                "log_prop_change_gn1_sd","log_prop_change_gn2_sd",
                                                                "log_prop_change_gn1_median","log_prop_change_gn2_median",
                                                                "log_prop_change_gn1_min","log_prop_change_gn2_min",
                                                                "log_prop_change_gn1_max","log_prop_change_gn2_max",
                                                                "log_prop_change_gn1_CoV","log_prop_change_gn2_CoV",
                                                                "PairID","Type","SERIES.n","SERIES.start","SERIES.end","YEAR.T","YEAR.T1")))) #makes an empty dataframe
    
    
    
    
    
    for (x in 1:length(n)){
      years.overlap=cont.overlapping[[n[x]]]
      
      start.year=years.overlap[1]
      end.year=years.overlap[length(years.overlap)]
      
      
      
      ## get measurements of abundance
      
      # Mean
      gen1.mean<-study_1$mean_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.mean<-study_2$mean_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # sd
      gen1.sd<-study_1$sd_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.sd<-study_2$sd_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # Median
      gen1.median<-study_1$median_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.median<-study_2$median_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # Min
      gen1.min<-study_1$min_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.min<-study_2$min_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # Max
      gen1.max<-study_1$max_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.max<-study_2$max_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # CoV
      gen1.CoV<-study_1$CoV_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.CoV<-study_2$CoV_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      
      
      
      
      #get year t and t+1
      year.Ts=years.overlap[1:(length(years.overlap)-1)]
      year.T1s=years.overlap[2:length(years.overlap)]
      
      
      
      
      
      #calculate  change in abundance measures
      
      # mean
      N1_gn1_mean = gen1.mean[2:length(gen1.mean)]
      N0_gn1_mean = gen1.mean[1:length(gen1.mean)-1]
      
      N1_gn2_mean = gen2.mean[2:length(gen2.mean)]
      N0_gn2_mean = gen2.mean[1:length(gen2.mean)-1]
      
      log_prop_change_gn1_mean = log(N1_gn1_mean/N0_gn1_mean)
      log_prop_change_gn2_mean = log(N1_gn2_mean/N0_gn2_mean)
      
      
      # sd
      N1_gn1_sd = gen1.sd[2:length(gen1.sd)]
      N0_gn1_sd = gen1.sd[1:length(gen1.sd)-1]
      
      N1_gn2_sd = gen2.sd[2:length(gen2.sd)]
      N0_gn2_sd = gen2.sd[1:length(gen2.sd)-1]
      
      log_prop_change_gn1_sd = log(N1_gn1_sd/N0_gn1_sd)
      log_prop_change_gn2_sd = log(N1_gn2_sd/N0_gn2_sd)
      
      
      # median
      N1_gn1_median = gen1.median[2:length(gen1.median)]
      N0_gn1_median = gen1.median[1:length(gen1.median)-1]
      
      N1_gn2_median = gen2.median[2:length(gen2.median)]
      N0_gn2_median = gen2.median[1:length(gen2.median)-1]
      
      log_prop_change_gn1_median = log(N1_gn1_median/N0_gn1_median)
      log_prop_change_gn2_median = log(N1_gn2_median/N0_gn2_median)
      
      
      # min
      N1_gn1_min = gen1.min[2:length(gen1.min)]
      N0_gn1_min = gen1.min[1:length(gen1.min)-1]
      
      N1_gn2_min = gen2.min[2:length(gen2.min)]
      N0_gn2_min = gen2.min[1:length(gen2.min)-1]
      
      log_prop_change_gn1_min = log(N1_gn1_min/N0_gn1_min)
      log_prop_change_gn2_min = log(N1_gn2_min/N0_gn2_min)
      
      
      # max
      N1_gn1_max = gen1.max[2:length(gen1.max)]
      N0_gn1_max = gen1.max[1:length(gen1.max)-1]
      
      N1_gn2_max = gen2.max[2:length(gen2.max)]
      N0_gn2_max = gen2.max[1:length(gen2.max)-1]
      
      log_prop_change_gn1_max = log(N1_gn1_max/N0_gn1_max)
      log_prop_change_gn2_max = log(N1_gn2_max/N0_gn2_max)
      
      
      # CoV
      N1_gn1_CoV = gen1.CoV[2:length(gen1.CoV)]
      N0_gn1_CoV = gen1.CoV[1:length(gen1.CoV)-1]
      
      N1_gn2_CoV = gen2.CoV[2:length(gen2.CoV)]
      N0_gn2_CoV = gen2.CoV[1:length(gen2.CoV)-1]
      
      log_prop_change_gn1_CoV = log(N1_gn1_CoV/N0_gn1_CoV)
      log_prop_change_gn2_CoV = log(N1_gn2_CoV/N0_gn2_CoV)
      
      
      
      
      
      
      
      res.add<-data.frame("Gn1"=genera.pair$Gn1,"Gn2"=genera.pair$Gn2,
                          
                          "log_prop_change_gn1_mean" = log_prop_change_gn1_mean,
                          "log_prop_change_gn2_mean" = log_prop_change_gn2_mean,
                          "log_prop_change_gn1_sd" = log_prop_change_gn1_sd,
                          "log_prop_change_gn2_sd" = log_prop_change_gn2_sd,
                          "log_prop_change_gn1_median" = log_prop_change_gn1_median,
                          "log_prop_change_gn2_median" = log_prop_change_gn1_median,
                          "log_prop_change_gn1_min" = log_prop_change_gn1_min,
                          "log_prop_change_gn2_min" = log_prop_change_gn1_min,
                          "log_prop_change_gn1_max" = log_prop_change_gn1_max,
                          "log_prop_change_gn2_max" = log_prop_change_gn1_max,
                          "log_prop_change_gn1_CoV" = log_prop_change_gn1_CoV,
                          "log_prop_change_gn2_CoV" = log_prop_change_gn1_CoV,
                          
                          "PairID"=genera.pair$PAIR.ID,
                          "Type"=genera.pair$Type,
                          "SERIES.n"=x,
                          "SERIES.start"=start.year,
                          "SERIES.end"=end.year,
                          "YEAR.T"=year.Ts,
                          "YEAR.T1"=year.T1s)
      
      res<-rbind(res,res.add)
    }
    
    res$SERIES.l=nrow(res)+1
    return(res)
    
    
  } #if multiple periods overlap at least 10 years calc for both!
  
  if (length(n)==1){
    years.overlap=cont.overlapping[[n]]
    
    start.year=years.overlap[1]
    end.year=years.overlap[length(years.overlap)]
    
    
    
    ## get measurements of abundance
    
    # Mean
    gen1.mean<-study_1$mean_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.mean<-study_2$mean_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # sd
    gen1.sd<-study_1$sd_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.sd<-study_2$sd_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # Median
    gen1.median<-study_1$median_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.median<-study_2$median_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # Min
    gen1.min<-study_1$min_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.min<-study_2$min_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # Max
    gen1.max<-study_1$max_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.max<-study_2$max_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # CoV
    gen1.CoV<-study_1$CoV_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.CoV<-study_2$CoV_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    
    
    #get year t and t+1
    year.Ts=years.overlap[1:(length(years.overlap)-1)]
    year.T1s=years.overlap[2:length(years.overlap)]
    
    
    
    
    #calculate  change in abundance measures
    
    # mean
    N1_gn1_mean = gen1.mean[2:length(gen1.mean)]
    N0_gn1_mean = gen1.mean[1:length(gen1.mean)-1]
    
    N1_gn2_mean = gen2.mean[2:length(gen2.mean)]
    N0_gn2_mean = gen2.mean[1:length(gen2.mean)-1]
    
    log_prop_change_gn1_mean = log(N1_gn1_mean/N0_gn1_mean)
    log_prop_change_gn2_mean = log(N1_gn2_mean/N0_gn2_mean)
    
    
    # sd
    N1_gn1_sd = gen1.sd[2:length(gen1.sd)]
    N0_gn1_sd = gen1.sd[1:length(gen1.sd)-1]
    
    N1_gn2_sd = gen2.sd[2:length(gen2.sd)]
    N0_gn2_sd = gen2.sd[1:length(gen2.sd)-1]
    
    log_prop_change_gn1_sd = log(N1_gn1_sd/N0_gn1_sd)
    log_prop_change_gn2_sd = log(N1_gn2_sd/N0_gn2_sd)
    
    
    # median
    N1_gn1_median = gen1.median[2:length(gen1.median)]
    N0_gn1_median = gen1.median[1:length(gen1.median)-1]
    
    N1_gn2_median = gen2.median[2:length(gen2.median)]
    N0_gn2_median = gen2.median[1:length(gen2.median)-1]
    
    log_prop_change_gn1_median = log(N1_gn1_median/N0_gn1_median)
    log_prop_change_gn2_median = log(N1_gn2_median/N0_gn2_median)
    
    
    # min
    N1_gn1_min = gen1.min[2:length(gen1.min)]
    N0_gn1_min = gen1.min[1:length(gen1.min)-1]
    
    N1_gn2_min = gen2.min[2:length(gen2.min)]
    N0_gn2_min = gen2.min[1:length(gen2.min)-1]
    
    log_prop_change_gn1_min = log(N1_gn1_min/N0_gn1_min)
    log_prop_change_gn2_min = log(N1_gn2_min/N0_gn2_min)
    
    
    # max
    N1_gn1_max = gen1.max[2:length(gen1.max)]
    N0_gn1_max = gen1.max[1:length(gen1.max)-1]
    
    N1_gn2_max = gen2.max[2:length(gen2.max)]
    N0_gn2_max = gen2.max[1:length(gen2.max)-1]
    
    log_prop_change_gn1_max = log(N1_gn1_max/N0_gn1_max)
    log_prop_change_gn2_max = log(N1_gn2_max/N0_gn2_max)
    
    
    # CoV
    N1_gn1_CoV = gen1.CoV[2:length(gen1.CoV)]
    N0_gn1_CoV = gen1.CoV[1:length(gen1.CoV)-1]
    
    N1_gn2_CoV = gen2.CoV[2:length(gen2.CoV)]
    N0_gn2_CoV = gen2.CoV[1:length(gen2.CoV)-1]
    
    log_prop_change_gn1_CoV = log(N1_gn1_CoV/N0_gn1_CoV)
    log_prop_change_gn2_CoV = log(N1_gn2_CoV/N0_gn2_CoV)
    
    
    
    
    res<-data.frame("Gn1"=genera.pair$Gn1,"Gn2"=genera.pair$Gn2,
                    
                    "log_prop_change_gn1_mean" = log_prop_change_gn1_mean,
                    "log_prop_change_gn2_mean" = log_prop_change_gn2_mean,
                    "log_prop_change_gn1_sd" = log_prop_change_gn1_sd,
                    "log_prop_change_gn2_sd" = log_prop_change_gn2_sd,
                    "log_prop_change_gn1_median" = log_prop_change_gn1_median,
                    "log_prop_change_gn2_median" = log_prop_change_gn2_median,
                    "log_prop_change_gn1_min" = log_prop_change_gn1_min,
                    "log_prop_change_gn2_min" = log_prop_change_gn2_min,
                    "log_prop_change_gn1_max" = log_prop_change_gn1_max,
                    "log_prop_change_gn2_max" = log_prop_change_gn2_max,
                    "log_prop_change_gn1_CoV" = log_prop_change_gn1_CoV,
                    "log_prop_change_gn2_CoV" = log_prop_change_gn2_CoV,
                    
                    "PairID"=genera.pair$PAIR.ID,
                    "Type"=genera.pair$Type,
                    "SERIES.n"=x,
                    "SERIES.start"=start.year,
                    "SERIES.end"=end.year,
                    "YEAR.T"=year.Ts,
                    "YEAR.T1"=year.T1s)
    
    
    
    res$SERIES.l=nrow(res)+1
    
    return(res)
  }}






# Get portion change for biomass variables


get.log.prop.change_biomass <-function(x,data,pairs.keep){
  
  #get pair
  genera.pair=pairs.keep[x,]
  
  #get studies
  study_1 <- data[which(data$ID==str_split(genera.pair$PAIR.ID,"_")[[1]][1]),] #get all timeseries for ID1
  study_2 <- data[which(data$ID==str_split(genera.pair$PAIR.ID,"_")[[1]][2]),] #get all timeseries for ID2
  
  #get years
  years1=study_1$YEAR[which(study_1$GENUS==genera.pair$Gn1)]
  years2=study_2$YEAR[which(study_2$GENUS==genera.pair$Gn2)]
  
  cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
  
  dat<-data.frame("Period"=NA,"value"=NA)
  for (y in 1:length(cont.overlapping)){
    dat[y,1]=y
    dat[y,2]=length(cont.overlapping[[y]])}
  
  years.overlap=cont.overlapping[[dat$Period[which(dat$value==max(dat$value))]]]
  
  n=dat$Period[which(dat$value>9)]
  
  if (length(n)>1){
    res=data.frame(matrix(ncol=21,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2", 
                                                                "log_prop_change_gn1_mean","log_prop_change_gn2_mean",
                                                                "log_prop_change_gn1_sd","log_prop_change_gn2_sd",
                                                                "log_prop_change_gn1_median","log_prop_change_gn2_median",
                                                                "log_prop_change_gn1_min","log_prop_change_gn2_min",
                                                                "log_prop_change_gn1_max","log_prop_change_gn2_max",
                                                                "log_prop_change_gn1_CoV","log_prop_change_gn2_CoV",
                                                                "PairID","Type","SERIES.n","SERIES.start","SERIES.end","YEAR.T","YEAR.T1")))) #makes an empty dataframe
    
    
    
    
    
    for (x in 1:length(n)){
      years.overlap=cont.overlapping[[n[x]]]
      
      start.year=years.overlap[1]
      end.year=years.overlap[length(years.overlap)]
      
      
      
      ## get measurements of abundance
      
      # Mean
      gen1.mean<-study_1$mean_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.mean<-study_2$mean_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # sd
      gen1.sd<-study_1$sd_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.sd<-study_2$sd_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # Median
      gen1.median<-study_1$median_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.median<-study_2$median_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # Min
      gen1.min<-study_1$min_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.min<-study_2$min_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # Max
      gen1.max<-study_1$max_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.max<-study_2$max_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      # CoV
      gen1.CoV<-study_1$CoV_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.CoV<-study_2$CoV_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      
      
      
      
      #get year t and t+1
      year.Ts=years.overlap[1:(length(years.overlap)-1)]
      year.T1s=years.overlap[2:length(years.overlap)]
      
      
      
      
      
      #calculate  change in biomass measures
      
      # mean
      N1_gn1_mean = gen1.mean[2:length(gen1.mean)]
      N0_gn1_mean = gen1.mean[1:length(gen1.mean)-1]
      
      N1_gn2_mean = gen2.mean[2:length(gen2.mean)]
      N0_gn2_mean = gen2.mean[1:length(gen2.mean)-1]
      
      log_prop_change_gn1_mean = log(N1_gn1_mean/N0_gn1_mean)
      log_prop_change_gn2_mean = log(N1_gn2_mean/N0_gn2_mean)
      
      
      # sd
      N1_gn1_sd = gen1.sd[2:length(gen1.sd)]
      N0_gn1_sd = gen1.sd[1:length(gen1.sd)-1]
      
      N1_gn2_sd = gen2.sd[2:length(gen2.sd)]
      N0_gn2_sd = gen2.sd[1:length(gen2.sd)-1]
      
      log_prop_change_gn1_sd = log(N1_gn1_sd/N0_gn1_sd)
      log_prop_change_gn2_sd = log(N1_gn2_sd/N0_gn2_sd)
      
      
      # median
      N1_gn1_median = gen1.median[2:length(gen1.median)]
      N0_gn1_median = gen1.median[1:length(gen1.median)-1]
      
      N1_gn2_median = gen2.median[2:length(gen2.median)]
      N0_gn2_median = gen2.median[1:length(gen2.median)-1]
      
      log_prop_change_gn1_median = log(N1_gn1_median/N0_gn1_median)
      log_prop_change_gn2_median = log(N1_gn2_median/N0_gn2_median)
      
      
      # min
      N1_gn1_min = gen1.min[2:length(gen1.min)]
      N0_gn1_min = gen1.min[1:length(gen1.min)-1]
      
      N1_gn2_min = gen2.min[2:length(gen2.min)]
      N0_gn2_min = gen2.min[1:length(gen2.min)-1]
      
      log_prop_change_gn1_min = log(N1_gn1_min/N0_gn1_min)
      log_prop_change_gn2_min = log(N1_gn2_min/N0_gn2_min)
      
      
      # max
      N1_gn1_max = gen1.max[2:length(gen1.max)]
      N0_gn1_max = gen1.max[1:length(gen1.max)-1]
      
      N1_gn2_max = gen2.max[2:length(gen2.max)]
      N0_gn2_max = gen2.max[1:length(gen2.max)-1]
      
      log_prop_change_gn1_max = log(N1_gn1_max/N0_gn1_max)
      log_prop_change_gn2_max = log(N1_gn2_max/N0_gn2_max)
      
      
      # CoV
      N1_gn1_CoV = gen1.CoV[2:length(gen1.CoV)]
      N0_gn1_CoV = gen1.CoV[1:length(gen1.CoV)-1]
      
      N1_gn2_CoV = gen2.CoV[2:length(gen2.CoV)]
      N0_gn2_CoV = gen2.CoV[1:length(gen2.CoV)-1]
      
      log_prop_change_gn1_CoV = log(N1_gn1_CoV/N0_gn1_CoV)
      log_prop_change_gn2_CoV = log(N1_gn2_CoV/N0_gn2_CoV)
      
      
      
      
      
      
      
      res.add<-data.frame("Gn1"=genera.pair$Gn1,"Gn2"=genera.pair$Gn2,
                          
                          "log_prop_change_gn1_mean" = log_prop_change_gn1_mean,
                          "log_prop_change_gn2_mean" = log_prop_change_gn2_mean,
                          "log_prop_change_gn1_sd" = log_prop_change_gn1_sd,
                          "log_prop_change_gn2_sd" = log_prop_change_gn2_sd,
                          "log_prop_change_gn1_median" = log_prop_change_gn1_median,
                          "log_prop_change_gn2_median" = log_prop_change_gn2_median,
                          "log_prop_change_gn1_min" = log_prop_change_gn1_min,
                          "log_prop_change_gn2_min" = log_prop_change_gn2_min,
                          "log_prop_change_gn1_max" = log_prop_change_gn1_max,
                          "log_prop_change_gn2_max" = log_prop_change_gn2_max,
                          "log_prop_change_gn1_CoV" = log_prop_change_gn1_CoV,
                          "log_prop_change_gn2_CoV" = log_prop_change_gn2_CoV,
                          
                          "PairID"=genera.pair$PAIR.ID,
                          "Type"=genera.pair$Type,
                          "SERIES.n"=x,
                          "SERIES.start"=start.year,
                          "SERIES.end"=end.year,
                          "YEAR.T"=year.Ts,
                          "YEAR.T1"=year.T1s)
      
      res<-rbind(res,res.add)
    }
    
    res$SERIES.l=nrow(res)+1
    return(res)
    
    
  } #if multiple periods overlap at least 10 years calc for both!
  
  if (length(n)==1){
    
    
    
    years.overlap=cont.overlapping[[n]]
    
    start.year=years.overlap[1]
    end.year=years.overlap[length(years.overlap)]
    
    
    
    ## get measurements of biomass
    
    # Mean
    gen1.mean<-study_1$mean_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.mean<-study_2$mean_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # sd
    gen1.sd<-study_1$sd_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.sd<-study_2$sd_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # Median
    gen1.median<-study_1$median_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.median<-study_2$median_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # Min
    gen1.min<-study_1$min_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.min<-study_2$min_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # Max
    gen1.max<-study_1$max_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.max<-study_2$max_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    # CoV
    gen1.CoV<-study_1$CoV_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.CoV<-study_2$CoV_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    
    
    #get year t and t+1
    year.Ts=years.overlap[1:(length(years.overlap)-1)]
    year.T1s=years.overlap[2:length(years.overlap)]
    
    
    
    
    #calculate  change in biomass measures
    
    # mean
    N1_gn1_mean = gen1.mean[2:length(gen1.mean)]
    N0_gn1_mean = gen1.mean[1:length(gen1.mean)-1]
    
    N1_gn2_mean = gen2.mean[2:length(gen2.mean)]
    N0_gn2_mean = gen2.mean[1:length(gen2.mean)-1]
    
    log_prop_change_gn1_mean = log(N1_gn1_mean/N0_gn1_mean)
    log_prop_change_gn2_mean = log(N1_gn2_mean/N0_gn2_mean)
    
    
    # sd
    N1_gn1_sd = gen1.sd[2:length(gen1.sd)]
    N0_gn1_sd = gen1.sd[1:length(gen1.sd)-1]
    
    N1_gn2_sd = gen2.sd[2:length(gen2.sd)]
    N0_gn2_sd = gen2.sd[1:length(gen2.sd)-1]
    
    log_prop_change_gn1_sd = log(N1_gn1_sd/N0_gn1_sd)
    log_prop_change_gn2_sd = log(N1_gn2_sd/N0_gn2_sd)
    
    
    # median
    N1_gn1_median = gen1.median[2:length(gen1.median)]
    N0_gn1_median = gen1.median[1:length(gen1.median)-1]
    
    N1_gn2_median = gen2.median[2:length(gen2.median)]
    N0_gn2_median = gen2.median[1:length(gen2.median)-1]
    
    log_prop_change_gn1_median = log(N1_gn1_median/N0_gn1_median)
    log_prop_change_gn2_median = log(N1_gn2_median/N0_gn2_median)
    
    
    # min
    N1_gn1_min = gen1.min[2:length(gen1.min)]
    N0_gn1_min = gen1.min[1:length(gen1.min)-1]
    
    N1_gn2_min = gen2.min[2:length(gen2.min)]
    N0_gn2_min = gen2.min[1:length(gen2.min)-1]
    
    log_prop_change_gn1_min = log(N1_gn1_min/N0_gn1_min)
    log_prop_change_gn2_min = log(N1_gn2_min/N0_gn2_min)
    
    
    # max
    N1_gn1_max = gen1.max[2:length(gen1.max)]
    N0_gn1_max = gen1.max[1:length(gen1.max)-1]
    
    N1_gn2_max = gen2.max[2:length(gen2.max)]
    N0_gn2_max = gen2.max[1:length(gen2.max)-1]
    
    log_prop_change_gn1_max = log(N1_gn1_max/N0_gn1_max)
    log_prop_change_gn2_max = log(N1_gn2_max/N0_gn2_max)
    
    
    # CoV
    N1_gn1_CoV = gen1.CoV[2:length(gen1.CoV)]
    N0_gn1_CoV = gen1.CoV[1:length(gen1.CoV)-1]
    
    N1_gn2_CoV = gen2.CoV[2:length(gen2.CoV)]
    N0_gn2_CoV = gen2.CoV[1:length(gen2.CoV)-1]
    
    log_prop_change_gn1_CoV = log(N1_gn1_CoV/N0_gn1_CoV)
    log_prop_change_gn2_CoV = log(N1_gn2_CoV/N0_gn2_CoV)
    
    
    
    
    res<-data.frame("Gn1"=genera.pair$Gn1,"Gn2"=genera.pair$Gn2,
                    
                    "log_prop_change_gn1_mean" = log_prop_change_gn1_mean,
                    "log_prop_change_gn2_mean" = log_prop_change_gn2_mean,
                    "log_prop_change_gn1_sd" = log_prop_change_gn1_sd,
                    "log_prop_change_gn2_sd" = log_prop_change_gn2_sd,
                    "log_prop_change_gn1_median" = log_prop_change_gn1_median,
                    "log_prop_change_gn2_median" = log_prop_change_gn2_median,
                    "log_prop_change_gn1_min" = log_prop_change_gn1_min,
                    "log_prop_change_gn2_min" = log_prop_change_gn2_min,
                    "log_prop_change_gn1_max" = log_prop_change_gn1_max,
                    "log_prop_change_gn2_max" = log_prop_change_gn2_max,
                    "log_prop_change_gn1_CoV" = log_prop_change_gn1_CoV,
                    "log_prop_change_gn2_CoV" = log_prop_change_gn2_CoV,
                    
                    "PairID"=genera.pair$PAIR.ID,
                    "Type"=genera.pair$Type,
                    "SERIES.n"=x,
                    "SERIES.start"=start.year,
                    "SERIES.end"=end.year,
                    "YEAR.T"=year.Ts,
                    "YEAR.T1"=year.T1s)
    
    
    
    res$SERIES.l=nrow(res)+1
    
    return(res)
  }}




















