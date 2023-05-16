## DATA PROCESSING FUNCTIONS
##created and updated by: AF
##############################################################


## Print progress bar for loops

set.prog.bar<-function(n_iter){
  #make progress bar
  progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                   total = n_iter,
                   complete = "=",
                   incomplete = "-",
                   current = ">",
                   clear = FALSE,
                   width = 100)
  
  
} 





#####################
# 3. OVERLAPPING YEARS
#####################

## Calculate years overlap of genera between timeseries

calc.overlap.between<-function(p,data){
  require(progress)
  require(dplyr)
  
  pair=bio.pairs[p,] #isolate which study pair we are focusing on
  
  
  #get timeseries
  timeseries_1 <- data[which(data$ID==pair$ID.1),] #get all timeseries for ID1
  timeseries_2 <- data[which(data$ID==pair$ID.2),] #get all timeseries for ID2
  
  #get genera lists
  genera1 <- timeseries_1$GENUS %>% unique()
  genera1 = genera1[!is.na(genera1)]
  
  genera2 <- timeseries_2$GENUS %>% unique()
  genera2 = genera2[!is.na(genera2)]
  
  if (length(genera1)==0){return(NULL)} #if no species found exit the function
  if (length(genera2)==0){return(NULL)} #if no species found exit the function
  
  #calc overlap for each genera-genera pairing
  mat<-data.frame("Gn1"=NA,"Gn2"=NA,"Max.Overlap"=NA) #make blank df
  mat<-mat[-1,] #remove 1st row
  for (i in 1:length(genera1)){
    gen1=genera1[i]
    #pb$tick()
    
    for (x in 1:length(genera2)){
      gen2=genera2[x]
      
      #get overlapping years
      years1=timeseries_1$YEAR[which(timeseries_1$GENUS==gen1)]
      years2=timeseries_2$YEAR[which(timeseries_2$GENUS==gen2)]
      
      cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
      
      dat<-data.frame("Period"=NA,"value"=NA)
      for (y in 1:length(cont.overlapping)){
        dat[y,1]=y
        dat[y,2]=length(cont.overlapping[[y]])}
      
      #take max
      mat<-rbind(mat,data.frame("Gn1"=gen1,"Gn2"=gen2,"Max.Overlap"=max(dat$value)))
      
    }
    
  } #for loop to look through genera
  mat$Type="Between" #add some stuff
  mat$PAIR.ID=paste0(pair$ID.1,"_",pair$ID.2) #add pair IDs
  
  return(mat) #exit the function
} #function to calc overlap between pairs BETWEEN two studies


## Calculate years overlap of genera within timeseries

calc.overlap.within<-function(ID,data){
  
  #get timeseries
  timeseries_1 <- data[which(data$ID==ID),]
  
  #get genera lists
  genera1 <- timeseries_1$GENUS %>% unique()
  genera1 = genera1[!is.na(genera1)]
  
  if (length(genera1)==0){return(NULL)}
  
  #then do ID1
  mat.id1<-data.frame("Gn1"=NA,"Gn2"=NA,"Max.Overlap"=NA)
  mat.id1<-mat.id1[-1,]
  for (i in 1:length(genera1)){
    #pb$tick()
    
    gen1=genera1[i]
    
    for (x in 1:length(genera1)){
      gen2=genera1[x]
      
      #get overlapping years
      years1=timeseries_1$YEAR[which(timeseries_1$GENUS==gen1)]
      years2=timeseries_1$YEAR[which(timeseries_1$GENUS==gen2)]
      
      cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
      
      dat<-data.frame("Period"=NA,"value"=NA)
      for (y in 1:length(cont.overlapping)){
        dat[y,1]=y
        dat[y,2]=length(cont.overlapping[[y]])}
      
      #take max
      mat.id1<-rbind(mat.id1,data.frame("Gn1"=gen1,"Gn2"=gen2,"Max.Overlap"=max(dat$value)))
      
      
    }
    
    
    #print(paste("DONE:",i,"out of",length(genera1)))
    
    
  }
  mat.id1$Type="Within"
  mat.id1$PAIR.ID=ID
  mat.id1<-mat.id1[-which(mat.id1$Gn1==mat.id1$Gn2),]
  
  return(mat.id1)
} #function to calc overlap between pairs WITHIN a single study



#####################
# 4. PROPORTION CHANGE
#####################

## Compute log proportion change

get.log.prop.change<-function(x,data,pairs.keep){
  
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
    res=data.frame(matrix(ncol=13,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2", "Log.prop.change.abun.Gn1","Log.prop.change.abun.Gn2","Log.prop.change.bio.Gn1","Log.prop.change.bio.Gn2","PairID","Type","SERIES.n","SERIES.start","SERIES.end","YEAR.T","YEAR.T1")))) #makes an empty dataframe
    for (x in 1:length(n)){
      years.overlap=cont.overlapping[[n[x]]]
      
      start.year=years.overlap[1]
      end.year=years.overlap[length(years.overlap)]
      
      #get abundances or biomass
      gen1.abun<-study_1$mean_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.abun<-study_2$mean_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      gen1.bio<-study_1$mean_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
      gen2.bio<-study_2$mean_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
      
      #get year t and t+1
      year.Ts=years.overlap[1:(length(years.overlap)-1)]
      year.T1s=years.overlap[2:length(years.overlap)]
      
      #calc  change in abun
      N1_gn1_abun = gen1.abun[2:length(gen1.abun)]
      N0_gn1_abun = gen1.abun[1:length(gen1.abun)-1]
      
      N1_gn2_abun = gen2.abun[2:length(gen2.abun)]
      N0_gn2_abun = gen2.abun[1:length(gen2.abun)-1]
      
      log_prop_change_gn1_abun = log(N1_gn1_abun/N0_gn1_abun)
      log_prop_change_gn2_abun = log(N1_gn2_abun/N0_gn2_abun)
      
      #calc  change in bio
      N1_gn1_bio = gen1.bio[2:length(gen1.bio)]
      N0_gn1_bio = gen1.bio[1:length(gen1.bio)-1]
      
      N1_gn2_bio = gen2.bio[2:length(gen2.bio)]
      N0_gn2_bio = gen2.bio[1:length(gen2.bio)-1]
      
      log_prop_change_gn1_bio = log(N1_gn1_bio/N0_gn1_bio)
      log_prop_change_gn2_bio = log(N1_gn2_bio/N0_gn2_bio)
      
      res.add<-data.frame("Gn1"=genera.pair$Gn1,"Gn2"=genera.pair$Gn2,
                          "Log.prop.change.abun.Gn1"=log_prop_change_gn1_abun,"Log.prop.change.abun.Gn2"=log_prop_change_gn2_abun,
                          "Log.prop.change.bio.Gn1"=log_prop_change_gn1_bio,"Log.prop.change.bio.Gn2"=log_prop_change_gn2_bio,
                          "PairID"=genera.pair$PAIR.ID,"Type"=genera.pair$Type,
                          "SERIES.n"=x,"SERIES.start"=start.year,"SERIES.end"=end.year,
                          "YEAR.T"=year.Ts,"YEAR.T1"=year.T1s)
      
      res<-rbind(res,res.add)
    }
    
    res$SERIES.l=nrow(res)+1
    return(res)
    
    
  } #if multiple periods overlap at least 10 years calc for both!
  
  if (length(n)==1){
    years.overlap=cont.overlapping[[n]]
    
    start.year=years.overlap[1]
    end.year=years.overlap[length(years.overlap)]
    
    #get abundances or biomass
    gen1.abun<-study_1$mean_abun_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.abun<-study_2$mean_abun_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    gen1.bio<-study_1$mean_bio_st[which(study_1$GENUS==genera.pair$Gn1 & study_1$YEAR%in%years.overlap)]
    gen2.bio<-study_2$mean_bio_st[which(study_2$GENUS==genera.pair$Gn2 & study_2$YEAR%in%years.overlap)]
    
    #get year t and t+1
    year.Ts=years.overlap[1:(length(years.overlap)-1)]
    year.T1s=years.overlap[2:length(years.overlap)]
    
    #calc  change in abun
    N1_gn1_abun = gen1.abun[2:length(gen1.abun)]
    N0_gn1_abun = gen1.abun[1:length(gen1.abun)-1]
    
    N1_gn2_abun = gen2.abun[2:length(gen2.abun)]
    N0_gn2_abun = gen2.abun[1:length(gen2.abun)-1]
    
    log_prop_change_gn1_abun = log(N1_gn1_abun/N0_gn1_abun)
    log_prop_change_gn2_abun = log(N1_gn2_abun/N0_gn2_abun)
    
    #calc  change in bio
    N1_gn1_bio = gen1.bio[2:length(gen1.bio)]
    N0_gn1_bio = gen1.bio[1:length(gen1.bio)-1]
    
    N1_gn2_bio = gen2.bio[2:length(gen2.bio)]
    N0_gn2_bio = gen2.bio[1:length(gen2.bio)-1]
    
    log_prop_change_gn1_bio = log(N1_gn1_bio/N0_gn1_bio)
    log_prop_change_gn2_bio = log(N1_gn2_bio/N0_gn2_bio)
    
    res<-data.frame("Gn1"=genera.pair$Gn1,"Gn2"=genera.pair$Gn2,
                    "Log.prop.change.abun.Gn1"=log_prop_change_gn1_abun,"Log.prop.change.abun.Gn2"=log_prop_change_gn2_abun,
                    "Log.prop.change.bio.Gn1"=log_prop_change_gn1_bio,"Log.prop.change.bio.Gn2"=log_prop_change_gn2_bio,
                    "PairID"=genera.pair$PAIR.ID,"Type"=genera.pair$Type,
                    "SERIES.n"=1,"SERIES.start"=start.year,"SERIES.end"=end.year,
                    "YEAR.T"=year.Ts,"YEAR.T1"=year.T1s)
    
    res$SERIES.l=nrow(res)+1
    
    return(res)}}


## Fetch meta data

make.meta<-function(data,meta){
  
  meta.data=data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("PairID", "ID1", "ID2","REALM1","REALM2","TAXA1","TAXA2","HABITAT1","HABITAT2","PROTECTED_AREA1","PROTECTED_AREA2")))) #makes an empty dataframe
  pb<-set.prog.bar(length(unique(data$PairID)))
  for (i in 1:length(unique(data$PairID))){
    pb$tick()
    
    PairID=unique(data$PairID)[i]
    
    ID1=str_split(PairID,"_")[[1]][1]
    ID2=str_split(PairID,"_")[[1]][2]
    
    met<-data.frame("PairID"=PairID,"ID1"=ID1,"ID2"=ID2,
                    "REALM1"=unique(meta$REALM[which(meta$ID==ID1)]),
                    "REALM2"=unique(meta$REALM[which(meta$ID==ID2)]),
                    "TAXA1"=unique(meta$TAXA[which(meta$ID==ID1)]),
                    "TAXA2"=unique(meta$TAXA[which(meta$ID==ID2)]),
                    "ORGANISMS1"=unique(meta$ORGANISMS[which(meta$ID==ID1)]),
                    "ORGANISMS2"=unique(meta$ORGANISMS[which(meta$ID==ID2)]),
                    "CLIMATE1"=unique(meta$CLIMATE[which(meta$ID==ID1)]),
                    "CLIMATE2"=unique(meta$CLIMATE[which(meta$ID==ID2)]))
    
    meta.data<-rbind(meta.data,met)
    
    
  }
  
  return(meta.data)
  
} #function to fetch meta data


