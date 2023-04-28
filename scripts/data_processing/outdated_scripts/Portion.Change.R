#Title: Calculate portion change for species pairs 
#Author: Isaac Eckert
#Date: April 26th, 2023

#Libraries
library(dplyr) 
library(tibble) 
library(progress)
library(stringr)

#Load data
load("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/preprocessing/cleaned_collated_standardized_MSF.Rdata")
between.studies.overlap<-readRDS("data/preprocessing/between.studies.overlap.RDS")
within.studies.overlap<-readRDS("data/preprocessing/within.studies.overlap.RDS")

bio.pairs <- read.csv("data/prep_biotime/bio_pairs_10km.csv")
meta.pairs <- read.csv("data/prep_biotime/meta_pairs_10km.csv")

pairs.keep=between.studies.overlap[which(between.studies.overlap$Max.Overlap>9),] #this is generated in the Overlap script
pairs.keep=within.studies.overlap[which(within.studies.overlap$Max.Overlap>9),] #this is generated in the Overlap script

#function to add progress bar
set.prog.bar<-function(n_iter){
  #make progress bar
  progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                   total = n_iter,
                   complete = "=",
                   incomplete = "-",
                   current = ">",
                   clear = FALSE,
                   width = 100)
  
} #function to make a progress bar

#function to calc prop change for genera pairing
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

results=data.frame(matrix(ncol=14,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2", "Log.prop.change.abun.Gn1","Log.prop.change.abun.Gn2","Log.prop.change.bio.Gn1","Log.prop.change.bio.Gn2","PairID","Type","SERIES.n","SERIES.start","SERIES.end","SERIES.l","YEAR.T","YEAR.T1")))) #makes an empty dataframe
pb<-set.prog.bar(nrow(pairs.keep))
for (i in 1:nrow(pairs.keep)){
  pb$tick()
  wk<-get.log.prop.change(i,collated.pairs_standardized_summary,pairs.keep)
  results<-rbind(results,wk)
  
} #get prop change for every genera pair

#add metadata
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
meta.data<-make.meta(results,collated.pairs_standardized_summary)
log.prop.change.with.meta.WITHIN<-left_join(results,meta.data) #join meta data with results df

#add unique genera ID col
log.prop.change.with.meta$UNIQUE.PAIR.ID=paste(log.prop.change.with.meta$Gn1,log.prop.change.with.meta$Gn2,log.prop.change.with.meta$PairID,sep="_")
log.prop.change.with.meta.WITHIN$UNIQUE.PAIR.ID=paste(log.prop.change.with.meta.WITHIN$Gn1,log.prop.change.with.meta.WITHIN$Gn2,log.prop.change.with.meta.WITHIN$PairID,sep="_")

#save
saveRDS(log.prop.change.with.meta,"data/log.prop.change.with.meta.RDS")
saveRDS(log.prop.change.with.meta.WITHIN,"data/preprocessing/log.prop.change.with.meta.WITHIN.RDS")

log.prop.change.with.meta<-readRDS("data/log.prop.change.with.meta.RDS")

#add better taxa cols
library(taxize)
taxa1<-tax_name(query = unique(log.prop.change.with.meta$Gn1), get = "class", db = "itis")
taxa2<-tax_name(query = unique(log.prop.change.with.meta$Gn2)[-which(unique(log.prop.change.with.meta$Gn2)%in%unique(log.prop.change.with.meta$Gn1))], get = "class", db = "itis")
taxa<-rbind(taxa1,taxa2)
#add new taxa col
taxa<-readRDS("data/class.assignment.RDS")
taxa$Gn1=taxa$query
taxa$Gn2=taxa$query
names(taxa)[which(names(taxa)=="class")]<-"RESOLVED.TAXA1"
taxa$RESOLVED.TAXA2=taxa$RESOLVED.TAXA1

log.prop.change.with.meta.taxa<-left_join(log.prop.change.with.meta,taxa[,c("Gn1","RESOLVED.TAXA1")])
log.prop.change.with.meta.taxa<-left_join(log.prop.change.with.meta.taxa,taxa[,c("Gn2","RESOLVED.TAXA2")])

log.prop.change.with.meta.WITHIN<-left_join(log.prop.change.with.meta.WITHIN,taxa[,c("Gn1","RESOLVED.TAXA1")])
log.prop.change.with.meta.WITHIN<-left_join(log.prop.change.with.meta.WITHIN,taxa[,c("Gn2","RESOLVED.TAXA2")])


saveRDS(log.prop.change.with.meta.taxa,"data/preprocessing/log.prop.change.with.meta.w.taxa.RDS")
saveRDS(log.prop.change.with.meta.WITHIN,"data/preprocessing/log.prop.change.with.meta.WITHIN.w.taxa.RDS")

table(log.prop.change.with.meta.WITHIN$RESOLVED.TAXA1)
table(log.prop.change.with.meta.WITHIN$RESOLVED.TAXA2)

#combine both and split for bio and abun
between<-readRDS("data/preprocessing/log.prop.change.with.meta.w.taxa.RDS")
witin<-readRDS("data/preprocessing/log.prop.change.with.meta.WITHIN.w.taxa.RDS")

all=rbind(between,witin)
remove(witin,between)

#split again for biomass and abun
abun=all[which(!is.na(all$Log.prop.change.abun.Gn1) & !is.na(all$Log.prop.change.abun.Gn2)),-which(names(all)%in%c("Log.prop.change.bio.Gn1","Log.prop.change.bio.Gn2"))] #keep only rows with non NA abundance values and get rid of bio cols
bio=all[which(!is.na(all$Log.prop.change.bio.Gn1) & !is.na(all$Log.prop.change.bio.Gn2)),-which(names(all)%in%c("Log.prop.change.abun.Gn1","Log.prop.change.abun.Gn2"))] #keep only rows with non NA biomass values and get rid of abun cols


names(abun)[which(names(abun)=="Log.prop.change.abun.Gn1")]<-"Prop.Change.Gn1"
names(abun)[which(names(abun)=="Log.prop.change.abun.Gn2")]<-"Prop.Change.Gn2"
names(bio)[which(names(bio)=="Log.prop.change.bio.Gn1")]<-"Prop.Change.Gn1"
names(bio)[which(names(bio)=="Log.prop.change.bio.Gn2")]<-"Prop.Change.Gn2"

abun$Metric="ABUNDANCE"
bio$Metric="BIOMASS"

full.data<-rbind(abun,bio)

saveRDS(full.data,"data/preprocessing/log.prop.change.full.data.RDS")





















