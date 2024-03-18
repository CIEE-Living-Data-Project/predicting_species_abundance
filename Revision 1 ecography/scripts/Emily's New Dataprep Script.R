#make new biotime data for the WG 
library(dplyr)
library(progress)

#load up biotime
biotime<-read.csv("~/Documents/Work and Career/LDP/Working Group/BioTIMEQuery_24_06_2021.csv")
metadata<-read.csv("~/Documents/Work and Career/LDP/Working Group/BioTIMEMetadata_24_06_2021.csv")

#remove marine and aquatic (keep only terr)
biotime.t<-biotime[which(biotime$STUDY_ID%in%metadata$STUDY_ID[which(metadata$REALM=="Terrestrial")]),]

#Fix study 39 before anything else
#Study 39
study_39 <- biotime.t %>%
  filter(STUDY_ID == "39")
unique(study_39$GENUS_SPECIES)
unique(study_39$GENUS)
study_39$GENUS[study_39$GENUS_SPECIES == "Yellow-throated Warbler"] <- "Setophaga"
study_39$GENUS[study_39$GENUS_SPECIES == "Black-capped chickadee"] <- "Poecile"
study_39$GENUS[study_39$GENUS_SPECIES == "Blue jay"] <- "Cyanocitta"
study_39$GENUS[study_39$GENUS_SPECIES == "Brown creeper"] <- "Certhia"
study_39$GENUS[study_39$GENUS_SPECIES == "Dark-eyed junco"] <- "Junco"
study_39$GENUS[study_39$GENUS_SPECIES == "Downy woodpecker"] <- "Picoides"
study_39$GENUS[study_39$GENUS_SPECIES == "Hairy woodpecker"] <- "Leuconotopicus"
study_39$GENUS[study_39$GENUS_SPECIES == "Ovenbird"] <- "Seiurus"
study_39$GENUS[study_39$GENUS_SPECIES == "Red-eyed vireo"] <- "Vireo"
study_39$GENUS[study_39$GENUS_SPECIES == "Rose-breasted grosbeak"] <- "Pheucticus"
study_39$GENUS[study_39$GENUS_SPECIES == "White-breasted nuthatch"] <- "Sitta"
study_39$GENUS[study_39$GENUS_SPECIES == "Winter wren"] <- "Troglodytes"
study_39$GENUS[study_39$GENUS_SPECIES == "Wood thrush"] <- "Hylocichla"
study_39$GENUS[study_39$GENUS_SPECIES == "Hermit Thrush"] <- "Catharus"
study_39$GENUS[study_39$GENUS_SPECIES == "Swainsons Thrush"] <- "Catharus"
study_39$GENUS[study_39$GENUS_SPECIES == "Yellow-bellied Sapsucker"] <- "Sphyrapicus"
study_39$GENUS[study_39$GENUS_SPECIES == "Solitary Vireo"] <- "Vireo"
study_39$GENUS[study_39$GENUS_SPECIES == "Red-breasted nuthatch"] <- "Sitta"
#Add back into biotime
biotime.t <- biotime.t %>%
  filter(!STUDY_ID == "39")
biotime.t <- bind_rows(biotime, study_39)
biotime.t <- biotime.t %>%
  arrange(STUDY_ID)

#Get study 221

study_221 <- biotime.t %>%
  filter(STUDY_ID == "221")
unique(study_221$GENUS)
#ChatGPT and research has provided genus guesses for each code, 
#in the context of this being a NA tree/shrub study 
#Reclassify the genera
map_genus <- function(code) {
  ifelse(grepl("SAL", code), "Salix", 
         ifelse(grepl("POPU", code), "Populus", 
                ifelse(grepl("BEL", code), "Betula", 
                       ifelse(grepl("ALNU", code), "Alnus", 
                              ifelse(grepl("BETU", code), "Betula", 
                                     ifelse(grepl("LARI", code), "Larix", 
                                            ifelse(grepl("MYRI", code), "Myrica", 
                                                   ifelse(grepl("PICE", code), "Picea",
                                                          ifelse(grepl("ROSA", code), "Rosa", 
                                                                 ifelse(grepl("RUBU", code), "Rubus", 
                                                                        ifelse(grepl("SHEP", code), "Shepherdia", 
                                                                               ifelse(grepl("VIBU", code), "Viburnum", 
                                                                                      code))))))))))))
}
study_221 <- study_221 %>% mutate(GENUS= map_genus(GENUS))
table(study_221$GENUS)
#Add back into biotime.t
biotime.t <- biotime.t %>%
  filter(!STUDY_ID=="221")
biotime.t <- bind_rows(biotime.t, study_221)
#Reorder
biotime.t <- biotime.t %>%
  arrange(STUDY_ID)





#any false 0s?
table(biotime.t$sum.allrawdata.ABUNDANCE[which(biotime.t$STUDY_ID%in%metadata$STUDY_ID[which(is.na(metadata$ABUNDANCE_TYPE))])])
#table(biotime.t$sum.allrawdata.BIOMASS[which(biotime.t$STUDY_ID%in%metadata$STUDY_ID[which(is.na(metadata$BIOMASS_TYPE))])])

#replace 0s with NAs if the study didnt record how they measured abudnance or biomass...
biotime.t$sum.allrawdata.ABUNDANCE[which(biotime.t$STUDY_ID%in%metadata$STUDY_ID[which(is.na(metadata$ABUNDANCE_TYPE))])]<-NA
#biotime.t$sum.allrawdata.BIOMASS[which(biotime.t$STUDY_ID%in%metadata$STUDY_ID[which(is.na(metadata$BIOMASS_TYPE))])]<-NA

#many plots are NAs which wont work during the aggregation step, lets make them a non-NA character
biotime.t$PLOT[which(is.na(biotime.t$PLOT))]<-"Was_NA"

#now we can aggregate to annual measures of abundance by averaging across plots in a study for each genera
biotime.t$UNIQUE_ID<-paste(biotime.t$STUDY_ID,biotime.t$PLOT,biotime.t$GENUS,biotime.t$YEAR,sep="~")
biotime.agg.abundance<-aggregate(cbind(sum.allrawdata.ABUNDANCE) ~ UNIQUE_ID, data = biotime.t[-which(is.na(biotime.t$sum.allrawdata.ABUNDANCE)),], mean)
#biotime.agg.biomass<-aggregate(cbind(sum.allrawdata.BIOMASS) ~ UNIQUE_ID, data = biotime.t[-which(is.na(biotime.t$sum.allrawdata.BIOMASS)),], mean)

#add back in some cols
biotime.agg.abundance$STUDY_ID<-biotime.t$STUDY_ID[match(biotime.agg.abundance$UNIQUE_ID,biotime.t$UNIQUE_ID)]
biotime.agg.abundance$PLOT<-biotime.t$PLOT[match(biotime.agg.abundance$UNIQUE_ID,biotime.t$UNIQUE_ID)]
biotime.agg.abundance$GENUS<-biotime.t$GENUS[match(biotime.agg.abundance$UNIQUE_ID,biotime.t$UNIQUE_ID)]
biotime.agg.abundance$YEAR<-biotime.t$YEAR[match(biotime.agg.abundance$UNIQUE_ID,biotime.t$UNIQUE_ID)]

# biotime.agg.biomass$STUDY_ID<-biotime.t$STUDY_ID[match(biotime.agg.biomass$UNIQUE_ID,biotime.t$UNIQUE_ID)]
# biotime.agg.biomass$PLOT<-biotime.t$PLOT[match(biotime.agg.biomass$UNIQUE_ID,biotime.t$UNIQUE_ID)]
# biotime.agg.biomass$GENUS<-biotime.t$GENUS[match(biotime.agg.biomass$UNIQUE_ID,biotime.t$UNIQUE_ID)]
# biotime.agg.biomass$YEAR<-biotime.t$YEAR[match(biotime.agg.biomass$UNIQUE_ID,biotime.t$UNIQUE_ID)]

#add in sample size (number of plots and number of species that went into calc'ing the mean)
ss_abun<-as.data.frame(table(biotime.t$UNIQUE_ID[-which(is.na(biotime.t$sum.allrawdata.ABUNDANCE))]))
#ss_bio<-as.data.frame(table(biotime.t$UNIQUE_ID[-which(is.na(biotime.t$sum.allrawdata.BIOMASS))]))
biotime.agg.abundance$SAMPLE_SIZE<-ss_abun$Freq[match(biotime.agg.abundance$UNIQUE_ID,ss_abun$Var1)]
#biotime.agg.biomass$SAMPLE_SIZE<-ss_bio$Freq[match(biotime.agg.biomass$UNIQUE_ID,ss_bio$Var1)]

head(ss_abun)

#remove all zero means
which(biotime.agg.abundance$sum.allrawdata.ABUNDANCE==0) #none
#length(which(biotime.agg.biomass$sum.allrawdata.BIOMASS==0)) #none

#biotime.agg.biomass<-biotime.agg.biomass[-which(biotime.agg.biomass$sum.allrawdata.BIOMASS==0),]
#dim(biotime.agg.biomass)

#add a new variable for study and plot
biotime.agg.abundance$STUDY_PLOT<-paste(biotime.agg.abundance$STUDY_ID,biotime.agg.abundance$PLOT,sep="~")
#biotime.agg.biomass$STUDY_PLOT<-paste(biotime.agg.biomass$STUDY_ID,biotime.agg.biomass$PLOT,sep="~")


########Emily's addition - fix genera names ##################


#Make a new frame with only study 39 and 221 
study_id_select <- c("39", "221")
study_221_39 <- biotime.agg.abundance %>%
  filter(STUDY_ID %in% study_id_select)

#find genera-genera pairs (within studies) that overlap for 10 years
find.genera.pairs.overlap.abundance<-function(study_plot){
  
  cropped<-study_221_39[which(study_221_39$STUDY_PLOT==study_plot),]

  #get genera list
  genera <- cropped$GENUS %>% unique()
  genera = genera[!is.na(genera)]
  
  if (length(genera)<2){return(NULL)} #if not enough genera found exit the function

  pairs<-as.data.frame(t(apply(combn(genera,2),2,paste))) #get all unique pairs
  names(pairs)=c("Gn1","Gn2")
  pairs$STUDY_PLOT=study_plot
  pairs$STUDY_ID=unique(cropped$STUDY)
  pairs$PLOT=unique(cropped$PLOT)
  pairs$Max.Overlap<-NA
  pairs$Type="Within"
  
  for (i in 1:nrow(pairs)){
    
    gen1=pairs$Gn1[i]
    gen2=pairs$Gn2[i]
    
    #get overlapping years
    years1=cropped$YEAR[which(cropped$GENUS==gen1)]
    years2=cropped$YEAR[which(cropped$GENUS==gen2)]
    
    cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
    
    dat<-data.frame("Period"=NA,"value"=NA)
    for (y in 1:length(cont.overlapping)){
      dat[y,1]=y
      dat[y,2]=length(cont.overlapping[[y]])}
    
    pairs$Max.Overlap[i]=max(dat$value)
    
  }
  
  return(pairs)
  
}#abundance one
# find.genera.pairs.overlap.biomass<-function(study_plot){
#   
#   cropped<-biotime.agg.biomass[which(biotime.agg.biomass$STUDY_PLOT==study_plot),]
#   
#   #get genera list
#   genera <- cropped$GENUS %>% unique()
#   genera = genera[!is.na(genera)]
#   
#   if (length(genera)<2){return(NULL)} #if not enough genera found exit the function
#   
#   pairs<-as.data.frame(t(apply(combn(genera,2),2,paste))) #get all unique pairs
#   names(pairs)=c("Gn1","Gn2")
#   pairs$STUDY_PLOT=study_plot
#   pairs$STUDY_ID=unique(cropped$STUDY)
#   pairs$PLOT=unique(cropped$PLOT)
#   pairs$Max.Overlap<-NA
#   pairs$Type="Within"
#   
#   for (i in 1:nrow(pairs)){
#     
#     gen1=pairs$Gn1[i]
#     gen2=pairs$Gn2[i]
#     
#     #get overlapping years
#     years1=cropped$YEAR[which(cropped$GENUS==gen1)]
#     years2=cropped$YEAR[which(cropped$GENUS==gen2)]
#     
#     cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
#     
#     dat<-data.frame("Period"=NA,"value"=NA)
#     for (y in 1:length(cont.overlapping)){
#       dat[y,1]=y
#       dat[y,2]=length(cont.overlapping[[y]])}
#     
#     pairs$Max.Overlap[i]=max(dat$value)
#     
#   }
#   
#   return(pairs)
#   
# }#biomass one (NOT RUN)

#run loop to calc all pairs
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

#Run abundance
pb<-set.prog.bar(length(unique(study_221_39$STUDY_PLOT))) #sets progress bar
within.studies.overlap.abundance=data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2","STUDY_PLOT","STUDY_ID","PLOT","Max.Overlap","Type")))) #makes an empty dataframe
for (h in 1:length(unique(study_221_39$STUDY_PLOT))){
  pb$tick()
  study_plot=unique(study_221_39$STUDY_PLOT)[h]
  data.pair<-find.genera.pairs.overlap.abundance(study_plot)
  within.studies.overlap.abundance<-rbind(within.studies.overlap.abundance,data.pair)
  
} #populates dataframe with pairs

# #Run biomass
# pb<-set.prog.bar(length(unique(biotime.agg.biomass$STUDY_PLOT))) #sets progress bar
# within.studies.overlap.biomass=data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2","STUDY_PLOT","STUDY_ID","PLOT","Max.Overlap","Type")))) #makes an empty dataframe
# for (h in 1:length(unique(biotime.agg.biomass$STUDY_PLOT))){
#   pb$tick()
#   study_plot=unique(biotime.agg.biomass$STUDY_PLOT)[h]
#   data.pair<-find.genera.pairs.overlap.biomass(study_plot)
#   within.studies.overlap.biomass<-rbind(within.studies.overlap.biomass,data.pair)
#   
# } #populates dataframe with pairs

#make lists of final genera pairs and their studies to use 
final.pairs.abundance<-within.studies.overlap.abundance[which(within.studies.overlap.abundance$Max.Overlap>9),c(3,1,2,4)]
# final.pairs.biomass<-within.studies.overlap.biomass[which(within.studies.overlap.biomass$Max.Overlap>9),c(3,1,2,4)]

dim(final.pairs.abundance)
# dim(final.pairs.biomass)

#add in some time series IDs
final.pairs.abundance$TS_ID<-paste(final.pairs.abundance$Gn1,final.pairs.abundance$Gn2,final.pairs.abundance$STUDY_PLOT,sep="~")
# final.pairs.biomass$TS_ID<-paste(final.pairs.biomass$Gn1,final.pairs.biomass$Gn2,final.pairs.biomass$STUDY_PLOT,sep="~")

#calculate prop change for final lists
get.log.prop.change.abundance<-function(ts){
  
  #get genera
  gen1=final.pairs.abundance$Gn1[which(final.pairs.abundance$TS_ID==ts)]
  gen2=final.pairs.abundance$Gn2[which(final.pairs.abundance$TS_ID==ts)]
  
  #get study
  study_plot=final.pairs.abundance$STUDY_PLOT[which(final.pairs.abundance$TS_ID==ts)]
  
  cropped=study_221_39[which(study_221_39$STUDY_PLOT==study_plot & study_221_39$GENUS%in%c(gen1,gen2)),]

  #get years
  years1=cropped$YEAR[which(cropped$GENUS==gen1)]
  years2=cropped$YEAR[which(cropped$GENUS==gen2)]
  
  cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
  
  dat<-data.frame("Period"=NA,"value"=NA)
  for (y in 1:length(cont.overlapping)){
    dat[y,1]=y
    dat[y,2]=length(cont.overlapping[[y]])}
  
  years.overlap=cont.overlapping[[dat$Period[which(dat$value==max(dat$value))]]] #these are the overlapping years
  
  n=dat$Period[which(dat$value>9)] #this is the period(s) to isolate
  
  if (length(n)>1){
    res=data.frame(matrix(ncol=20,nrow=0, dimnames=list(NULL, c("TS_ID","STUDY_PLOT","Gn1", "Gn2", "Log.prop.change.Gn1","Log.prop.change.Gn2","SERIES.n","SERIES.start","SERIES.end","YEAR.T","YEAR.T1","METRIC","Abs.change.Gn1","Abs.change.Gn2","Abun.T.Gn1","Abun.T1.Gn1","Abun.T.Gn2","Abun.T1.Gn2","Length.Unique.Values.Gn1","Length.Unique.Values.Gn2")))) #makes an empty dataframe
    for (x in 1:length(n)){
      
      years.overlap=cont.overlapping[[n[x]]]
      
      start.year=years.overlap[1]
      end.year=years.overlap[length(years.overlap)]
      
      #get measures
      gen1.abun<-cropped$sum.allrawdata.ABUNDANCE[which(cropped$GENUS==gen1 & cropped$YEAR%in%years.overlap)]
      gen2.abun<-cropped$sum.allrawdata.ABUNDANCE[which(cropped$GENUS==gen2 & cropped$YEAR%in%years.overlap)]
      
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
      
      abs_change_gn1_abun = N1_gn1_abun-N0_gn1_abun
      abs_change_gn2_abun = N1_gn2_abun-N0_gn2_abun
      
      res.add<-data.frame("TS_ID"=ts,"STUDY_PLOT"=study_plot,"Gn1"=gen1,"Gn2"=gen2,
                          "Log.prop.change.Gn1"=log_prop_change_gn1_abun,"Log.prop.change.Gn2"=log_prop_change_gn2_abun,
                          "SERIES.n"=x,"SERIES.start"=start.year,"SERIES.end"=end.year,
                          "YEAR.T"=year.Ts,"YEAR.T1"=year.T1s,"METRIC"="ABUNDANCE",
                          "Abs.change.Gn1"=abs_change_gn1_abun,"Abs.change.Gn2"=abs_change_gn2_abun,
                          "Abun.T.Gn1"=N0_gn1_abun,"Abun.T1.Gn1"=N1_gn1_abun,
                          "Abun.T.Gn2"=N0_gn2_abun,"Abun.T1.Gn2"=N1_gn2_abun,
                          "Length.Unique.Values.Gn1"=length(unique(abs(log_prop_change_gn1_abun))),
                          "Length.Unique.Values.Gn2"=length(unique(abs(log_prop_change_gn2_abun))))
      
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
    gen1.abun<-cropped$sum.allrawdata.ABUNDANCE[which(cropped$GENUS==gen1 & cropped$YEAR%in%years.overlap)]
    gen2.abun<-cropped$sum.allrawdata.ABUNDANCE[which(cropped$GENUS==gen2 & cropped$YEAR%in%years.overlap)]
    
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
    
    abs_change_gn1_abun = N1_gn1_abun-N0_gn1_abun
    abs_change_gn2_abun = N1_gn2_abun-N0_gn2_abun
    
    res<-data.frame("TS_ID"=ts,"STUDY_PLOT"=study_plot,"Gn1"=gen1,"Gn2"=gen2,
                    "Log.prop.change.Gn1"=log_prop_change_gn1_abun,"Log.prop.change.Gn2"=log_prop_change_gn2_abun,
                    "SERIES.n"=1,"SERIES.start"=start.year,"SERIES.end"=end.year,
                    "YEAR.T"=year.Ts,"YEAR.T1"=year.T1s,"METRIC"="ABUNDANCE",
                    "Abs.change.Gn1"=abs_change_gn1_abun,"Abs.change.Gn2"=abs_change_gn2_abun,
                    "Abun.T.Gn1"=N0_gn1_abun,"Abun.T1.Gn1"=N1_gn1_abun,
                    "Abun.T.Gn2"=N0_gn2_abun,"Abun.T1.Gn2"=N1_gn2_abun,
                    "Length.Unique.Values.Gn1"=length(unique(abs(log_prop_change_gn1_abun))),
                    "Length.Unique.Values.Gn2"=length(unique(abs(log_prop_change_gn2_abun))))
    
    res$SERIES.l=nrow(res)+1
    
    return(res)}}
#get.log.prop.change.biomass<-function(ts){
  # 
  # #get genera
  # gen1=final.pairs.biomass$Gn1[which(final.pairs.biomass$TS_ID==ts)]
  # gen2=final.pairs.biomass$Gn2[which(final.pairs.biomass$TS_ID==ts)]
  # 
  # #get study
  # study_plot=final.pairs.biomass$STUDY_PLOT[which(final.pairs.biomass$TS_ID==ts)]
  # 
  # cropped=biotime.agg.biomass[which(biotime.agg.biomass$STUDY_PLOT==study_plot & biotime.agg.biomass$GENUS%in%c(gen1,gen2)),]
  # 
  # #get years
  # years1=cropped$YEAR[which(cropped$GENUS==gen1)]
  # years2=cropped$YEAR[which(cropped$GENUS==gen2)]
  # 
  # cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
  # 
  # dat<-data.frame("Period"=NA,"value"=NA)
  # for (y in 1:length(cont.overlapping)){
  #   dat[y,1]=y
  #   dat[y,2]=length(cont.overlapping[[y]])}
  # 
  # years.overlap=cont.overlapping[[dat$Period[which(dat$value==max(dat$value))]]] #these are the overlapping years
  # 
  # n=dat$Period[which(dat$value>9)] #this is the period(s) to isolate
  # 
  # if (length(n)>1){
  #   res=data.frame(matrix(ncol=12,nrow=0, dimnames=list(NULL, c("TS_ID","STUDY_PLOT","Gn1", "Gn2", "Log.prop.change.Gn1","Log.prop.change.Gn2","SERIES.n","SERIES.start","SERIES.end","YEAR.T","YEAR.T1","METRIC")))) #makes an empty dataframe
  #   for (x in 1:length(n)){
  #     years.overlap=cont.overlapping[[n[x]]]
  #     
  #     start.year=years.overlap[1]
  #     end.year=years.overlap[length(years.overlap)]
  #     
  #     #get measures
  #     gen1.abun<-cropped$sum.allrawdata.BIOMASS[which(cropped$GENUS==gen1 & cropped$YEAR%in%years.overlap)]
  #     gen2.abun<-cropped$sum.allrawdata.BIOMASS[which(cropped$GENUS==gen2 & cropped$YEAR%in%years.overlap)]
  #     
  #     #get year t and t+1
  #     year.Ts=years.overlap[1:(length(years.overlap)-1)]
  #     year.T1s=years.overlap[2:length(years.overlap)]
  #     
  #     #calc  change in abun
  #     N1_gn1_abun = gen1.abun[2:length(gen1.abun)]
  #     N0_gn1_abun = gen1.abun[1:length(gen1.abun)-1]
  #     
  #     N1_gn2_abun = gen2.abun[2:length(gen2.abun)]
  #     N0_gn2_abun = gen2.abun[1:length(gen2.abun)-1]
  #     
  #     log_prop_change_gn1_abun = log(N1_gn1_abun/N0_gn1_abun)
  #     log_prop_change_gn2_abun = log(N1_gn2_abun/N0_gn2_abun)
  #     
  #     
  #     res.add<-data.frame("TS_ID"=ts,"STUDY_PLOT"=study_plot,"Gn1"=gen1,"Gn2"=gen2,
  #                         "Log.prop.change.Gn1"=log_prop_change_gn1_abun,"Log.prop.change.Gn2"=log_prop_change_gn2_abun,
  #                         "SERIES.n"=x,"SERIES.start"=start.year,"SERIES.end"=end.year,
  #                         "YEAR.T"=year.Ts,"YEAR.T1"=year.T1s,"METRIC"="BIOMASS")
  #     
  #     res<-rbind(res,res.add)
  #   }
  #   
  #   res$SERIES.l=nrow(res)+1
  #   return(res)
  #   
  #   
  # } #if multiple periods overlap at least 10 years calc for both!
  # 
  # if (length(n)==1){
  #   years.overlap=cont.overlapping[[n]]
  #   
  #   start.year=years.overlap[1]
  #   end.year=years.overlap[length(years.overlap)]
  #   
  #   #get abundances or biomass
  #   gen1.abun<-cropped$sum.allrawdata.BIOMASS[which(cropped$GENUS==gen1 & cropped$YEAR%in%years.overlap)]
  #   gen2.abun<-cropped$sum.allrawdata.BIOMASS[which(cropped$GENUS==gen2 & cropped$YEAR%in%years.overlap)]
  #   
  #   #get year t and t+1
  #   year.Ts=years.overlap[1:(length(years.overlap)-1)]
  #   year.T1s=years.overlap[2:length(years.overlap)]
  #   
  #   #calc  change in abun
  #   N1_gn1_abun = gen1.abun[2:length(gen1.abun)]
  #   N0_gn1_abun = gen1.abun[1:length(gen1.abun)-1]
  #   
  #   N1_gn2_abun = gen2.abun[2:length(gen2.abun)]
  #   N0_gn2_abun = gen2.abun[1:length(gen2.abun)-1]
  #   
  #   log_prop_change_gn1_abun = log(N1_gn1_abun/N0_gn1_abun)
  #   log_prop_change_gn2_abun = log(N1_gn2_abun/N0_gn2_abun)
  #   
  #   res<-data.frame("TS_ID"=ts,"STUDY_PLOT"=study_plot,"Gn1"=gen1,"Gn2"=gen2,
  #                   "Log.prop.change.Gn1"=log_prop_change_gn1_abun,"Log.prop.change.Gn2"=log_prop_change_gn2_abun,
  #                   "SERIES.n"=1,"SERIES.start"=start.year,"SERIES.end"=end.year,
  #                   "YEAR.T"=year.Ts,"YEAR.T1"=year.T1s,"METRIC"="BIOMASS")
  #   
  #   res$SERIES.l=nrow(res)+1
  #   
  #   return(res)}}

#loop through all abundance 
results.abundance=data.frame(matrix(ncol=21,nrow=0, dimnames=list(NULL, c("TS_ID","STUDY_PLOT","Gn1", "Gn2", "Log.prop.change.Gn1","Log.prop.change.Gn2","SERIES.n","SERIES.start","SERIES.end","YEAR.T","YEAR.T1","METRIC","Abs.change.Gn1","Abs.change.Gn2","Abun.T.Gn1","Abun.T1.Gn1","Abun.T.Gn2","Abun.T1.Gn2","Length.Unique.Values.Gn1","Length.Unique.Values.Gn2","SERIES.l")))) #makes an empty dataframe
pb<-set.prog.bar(nrow(final.pairs.abundance))
for (i in 1:nrow(final.pairs.abundance)){
  pb$tick()
  
  ts=final.pairs.abundance$TS_ID[i]
  wk<-get.log.prop.change.abundance(ts)
  results.abundance<-rbind(results.abundance,wk)
  
} #get prop change for every genera pair

write.csv(results.abundance,"~/Documents/Work and Career/LDP/Working Group/results.abundance_221_29.csv")
results.abundance=read.csv("C:\\Users\\isaac\\OneDrive - McGill University\\Sharepoint\\Datasets\\BioTIME\\results.abundance.csv")

#add in sp per genus/plot/year
b<-biotime.t
b$UNIQUE_ID_SP<-paste(b$UNIQUE_ID,b$SPECIES,sep="~")
b<-b[-which(duplicated(b$UNIQUE_ID_SP)),]

data<-as.data.frame(table(b$UNIQUE_ID[-which(is.na(b$sum.allrawdata.ABUNDANCE))]))
names(data)=c("UNIQUE_ID","Richness")

#add cols
results.abundance$Rich.T.Gn1=data$Richness[match(paste(results.abundance$STUDY_PLOT,results.abundance$Gn1,results.abundance$YEAR.T,sep="~"),data$UNIQUE_ID)]
results.abundance$Rich.T1.Gn1=data$Richness[match(paste(results.abundance$STUDY_PLOT,results.abundance$Gn1,results.abundance$YEAR.T1,sep="~"),data$UNIQUE_ID)]
results.abundance$Rich.T.Gn2=data$Richness[match(paste(results.abundance$STUDY_PLOT,results.abundance$Gn2,results.abundance$YEAR.T,sep="~"),data$UNIQUE_ID)]
results.abundance$Rich.T1.Gn2=data$Richness[match(paste(results.abundance$STUDY_PLOT,results.abundance$Gn2,results.abundance$YEAR.T1,sep="~"),data$UNIQUE_ID)]

#add total sample size
ss_abun<-as.data.frame(table(biotime.t$UNIQUE_ID[-which(is.na(biotime.t$sum.allrawdata.ABUNDANCE))]))
names(ss_abun)=c("UNIQUE_ID","Sample.Size")

#add cols
results.abundance$SS.T.Gn1=data$Richness[match(paste(results.abundance$STUDY_PLOT,results.abundance$Gn1,results.abundance$YEAR.T,sep="~"),data$UNIQUE_ID)]
results.abundance$SS.T1.Gn1=data$Richness[match(paste(results.abundance$STUDY_PLOT,results.abundance$Gn1,results.abundance$YEAR.T1,sep="~"),data$UNIQUE_ID)]
results.abundance$SS.T.Gn2=data$Richness[match(paste(results.abundance$STUDY_PLOT,results.abundance$Gn2,results.abundance$YEAR.T,sep="~"),data$UNIQUE_ID)]
results.abundance$SS.T1.Gn2=data$Richness[match(paste(results.abundance$STUDY_PLOT,results.abundance$Gn2,results.abundance$YEAR.T1,sep="~"),data$UNIQUE_ID)]



#loop through all biomass 
#results.biomass=data.frame(matrix(ncol=13,nrow=0, dimnames=list(NULL, c("TS_ID","STUDY_PLOT","Gn1", "Gn2", "Log.prop.change.Gn1","Log.prop.change.Gn2","SERIES.n","SERIES.start","SERIES.end","YEAR.T","YEAR.T1","METRIC","SERIES.l")))) #makes an empty dataframe
#pb<-set.prog.bar(nrow(final.pairs.biomass))
#for (i in 1:nrow(final.pairs.biomass)){
#   pb$tick()
#   
#   ts=final.pairs.biomass$TS_ID[i]
#   wk<-get.log.prop.change.biomass(ts)
#   results.biomass<-rbind(results.biomass,wk)
#   
# } #get prop change for every genera pair

#combine
#results=rbind(results.abundance,results.biomass)

#add in metadata
results.abundance$STUDY_ID<-gsub("\\~.*","",results.abundance$STUDY_PLOT)
results.abundance$REALM<-metadata$REALM[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]
results.abundance$CLIMATE<-metadata$CLIMATE[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]
results.abundance$HABITAT<-metadata$HABITAT[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]
results.abundance$PROTECTED_AREA<-metadata$PROTECTED_AREA[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]
results.abundance$TAXA<-metadata$TAXA[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]
results.abundance$AREA_SQ_KM<-metadata$AREA_SQ_KM[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]
results.abundance$ORGANISMS<-metadata$ORGANISMS[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]
results.abundance$LATITUDE<-metadata$CENT_LAT[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]
results.abundance$LONGITUDE<-metadata$CENT_LONG[match(results.abundance$STUDY_ID,metadata$STUDY_ID)]

dim(results.abundance)
length(unique(results.abundance$TS_ID)) #527874
length(unique(results.abundance$STUDY_ID)) #69

#save
write.csv(results.abundance,"~/Documents/Work and Career/LDP/Working Group/results.abundance_221_39.csv")
results<-read.csv("/Users/isaaceckert/Desktop/within.study.updated.data.csv")


table(results$METRIC)




