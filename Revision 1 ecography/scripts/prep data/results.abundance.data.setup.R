# Biotime cleaning, aggregating, and change in abundance calculations
# Studies 39 and 221 corrected (note that 225 still contains errors)
#Studies 39, 221, 225 ultimately removed from final analyses 
# Date created: Jan 24 2024
# Date updated: 18 April 2024

# Tested and checked by NC, 24.04.2024



#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_





#make new biotime data for the WG 
library(dplyr)
library(progress)

#load up biotime
biotime<-read.csv("~/Documents/Work and Career/LDP/Working Group/BioTIMEQuery_24_06_2021.csv")
metadata<-read.csv("~/Documents/Work and Career/LDP/Working Group/BioTIMEMetadata_24_06_2021.csv")

# # NC LOAD
# biotime<-read.csv('~/Desktop/Code/predicting_species_abundance/data/BioTIMECitations_24_06_2021.csv')
# metadata<-read.csv('data/BioTIMEMetadata_24_06_2021.csv')

# #NC CHECKS
# length(unique(biotime$STUDY_ID))
# unique(metadata$REALM)

#remove marine and aquatic (keep only terr)
biotime.t<-biotime[which(biotime$STUDY_ID %in% 
                           metadata$STUDY_ID[which(
                             metadata$REALM=="Terrestrial")]),]

# #NC CHECKS
# foo <- metadata$STUDY_ID[which(metadata$REALM=="Terrestrial")]
# length(metadata$STUDY_ID[which(metadata$REALM=="Terrestrial")])
# nrow(biotime.t)
# length(unique(biotime.t$STUDY_ID))


###### ENB addition - cleaning studies 39, 221 #######
#Fix study 39 before cleaning, aggregation
#Get only study 39
study_39 <- biotime.t %>%
  filter(STUDY_ID == "39")
#Summarize unique genera and species 
unique(study_39$GENUS_SPECIES)
unique(study_39$GENUS)
#Correct erroneous species names 
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
#Add back into full biotime dataset
biotime.t <- biotime.t %>%
  filter(!STUDY_ID == "39")
biotime.t <- bind_rows(biotime, study_39)
biotime.t <- biotime.t %>%
  arrange(STUDY_ID)

#Correct study 221
#Filter to only study 221
study_221 <- biotime.t %>%
  filter(STUDY_ID == "221")
#Summarize unique genera 
unique(study_221$GENUS)
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
                                                                               ifelse(grepl("VIBU", code), 
                                                                                      
                                                                                      #NC CHECK: should this also have 'ifelse'?
                                                                                      "Viburnum", 
                                                                                      code))))))))))))
}
study_221 <- study_221 %>% mutate(GENUS= map_genus(GENUS))
#Summarize unique and corrected genera
table(study_221$GENUS)
#Add back into biotime.t
biotime.t <- biotime.t %>%
  filter(!STUDY_ID=="221")
biotime.t <- bind_rows(biotime.t, study_221)
#Reorder to ascending study ID
biotime.t <- biotime.t %>%
  arrange(STUDY_ID)



##### Back to IE code ####

#any false 0s?
table(biotime.t$sum.allrawdata.ABUNDANCE[which(
  biotime.t$STUDY_ID%in%
    metadata$STUDY_ID[which(is.na(metadata$ABUNDANCE_TYPE))])])
#table(biotime.t$sum.allrawdata.BIOMASS[which(biotime.t$STUDY_ID%in%metadata$STUDY_ID[which(is.na(metadata$BIOMASS_TYPE))])])

#replace 0s with NAs if the study didnt record how they measured abudnance or biomass...
biotime.t$sum.allrawdata.ABUNDANCE[which(
  biotime.t$STUDY_ID%in%
    metadata$STUDY_ID[which(is.na(metadata$ABUNDANCE_TYPE))])] <- NA
#biotime.t$sum.allrawdata.BIOMASS[which(biotime.t$STUDY_ID%in%metadata$STUDY_ID[which(is.na(metadata$BIOMASS_TYPE))])]<-NA

#many plots are NAs which wont work during the aggregation step, lets make them a non-NA character
biotime.t$PLOT[which(is.na(biotime.t$PLOT))]<-"Was_NA"

#now we can aggregate to annual measures of abundance by averaging across plots in a study for each genera
biotime.t$UNIQUE_ID<-paste(biotime.t$STUDY_ID,biotime.t$PLOT,biotime.t$GENUS,biotime.t$YEAR,sep="~")
biotime.agg.abundance<-
  
  #NC CHECK: why using cbind here? and why not using the built-in 'na.action' for aggregate?
  aggregate(cbind(sum.allrawdata.ABUNDANCE) ~ UNIQUE_ID, 
            data = biotime.t[-which(is.na(biotime.t$sum.allrawdata.ABUNDANCE)),], mean)
#biotime.agg.biomass<-aggregate(cbind(sum.allrawdata.BIOMASS) ~ UNIQUE_ID, data = biotime.t[-which(is.na(biotime.t$sum.allrawdata.BIOMASS)),], mean)

#add back in identifiying cols
biotime.agg.abundance$STUDY_ID<-biotime.t$STUDY_ID[match(
  biotime.agg.abundance$UNIQUE_ID,biotime.t$UNIQUE_ID)]
biotime.agg.abundance$PLOT<-biotime.t$PLOT[match(
  biotime.agg.abundance$UNIQUE_ID,biotime.t$UNIQUE_ID)]
biotime.agg.abundance$GENUS<-biotime.t$GENUS[match(
  biotime.agg.abundance$UNIQUE_ID,biotime.t$UNIQUE_ID)]
biotime.agg.abundance$YEAR<-biotime.t$YEAR[match(
  biotime.agg.abundance$UNIQUE_ID,biotime.t$UNIQUE_ID)]

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

#NC CHECK: need to formally test this function if not already done so
#find genera-genera pairs (within studies) that overlap for 10 years
find.genera.pairs.overlap.abundance<-function(study_plot){
  
  cropped<-biotime.agg.abundance[which(biotime.agg.abundance$STUDY_PLOT==study_plot),]

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
    
    cont.overlapping<-split(years1[which(years1%in%years2)], 
                            cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
    
    dat<-data.frame("Period"=NA,"value"=NA)
    for (y in 1:length(cont.overlapping)){
      dat[y,1]=y
      dat[y,2]=length(cont.overlapping[[y]])}
    
    pairs$Max.Overlap[i]=max(dat$value)
    
  }
  
  return(pairs)
  
}#abundance one



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
pb<-set.prog.bar(length(unique(biotime.agg.abundance$STUDY_PLOT))) #sets progress bar
within.studies.overlap.abundance=data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2","STUDY_PLOT","STUDY_ID","PLOT","Max.Overlap","Type")))) #makes an empty dataframe
for (h in 1:length(unique(biotime.agg.abundance$STUDY_PLOT))){
  pb$tick()
  study_plot=unique(biotime.agg.abundance$STUDY_PLOT)[h]
  data.pair<-find.genera.pairs.overlap.abundance(study_plot)
  within.studies.overlap.abundance<-rbind(within.studies.overlap.abundance,data.pair)
  
} #populates dataframe with pairs



#make lists of final genera pairs and their studies to use 
final.pairs.abundance<-within.studies.overlap.abundance[which(within.studies.overlap.abundance$Max.Overlap>9),c(3,1,2,4)]

dim(final.pairs.abundance)

#add in some time series IDs
final.pairs.abundance$TS_ID<-paste(final.pairs.abundance$Gn1,final.pairs.abundance$Gn2,final.pairs.abundance$STUDY_PLOT,sep="~")

#calculate prop change for final lists
get.log.prop.change.abundance<-function(ts){
  
  #get genera
  gen1=final.pairs.abundance$Gn1[which(final.pairs.abundance$TS_ID==ts)]
  gen2=final.pairs.abundance$Gn2[which(final.pairs.abundance$TS_ID==ts)]
  
  #get study
  study_plot=final.pairs.abundance$STUDY_PLOT[which(final.pairs.abundance$TS_ID==ts)]
  
  cropped=biotime.agg.abundance[which(biotime.agg.abundance$STUDY_PLOT==study_plot & biotime.agg.abundance$GENUS%in%c(gen1,gen2)),]

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


#loop through all abundance 
results.abundance=data.frame(matrix(ncol=21,nrow=0, dimnames=list(NULL, c("TS_ID","STUDY_PLOT","Gn1", "Gn2", "Log.prop.change.Gn1","Log.prop.change.Gn2","SERIES.n","SERIES.start","SERIES.end","YEAR.T","YEAR.T1","METRIC","Abs.change.Gn1","Abs.change.Gn2","Abun.T.Gn1","Abun.T1.Gn1","Abun.T.Gn2","Abun.T1.Gn2","Length.Unique.Values.Gn1","Length.Unique.Values.Gn2","SERIES.l")))) #makes an empty dataframe
pb<-set.prog.bar(nrow(final.pairs.abundance))
for (i in 1:nrow(final.pairs.abundance)){
  pb$tick()
  
  ts=final.pairs.abundance$TS_ID[i]
  wk<-get.log.prop.change.abundance(ts)
  results.abundance<-rbind(results.abundance,wk)
  
} #get prop change for every genera pair

#Note: save locations should be updated if running code on own device 
#write.csv(results.abundance,"~/Documents/Work and Career/LDP/Working Group/results.abundance.csv")
#results.abundance=read.csv("~/Documents/Work and Career/LDP/Working Group/results.abundance.csv")







###### ENB addition to code: calculating richness for studies with incorrect names ####
#Note: right now only calculates richness for corrected studies 39 and 221
#Initially calculated richness for all studies, written by IE
#Can be easily adjusted to work for all studies
#Current data used for analyses has richness and sample size for all studies, 
#Including corrected 39 and 221 


#Read in full biotime data
biotime<-read.csv("~/Documents/Work and Career/LDP/Working Group/BioTIMEQuery_24_06_2021.csv")
#Not available in github - subset of complete data that was shared during analysis. 
#Corrected results available in Github
missing_data <- read.csv("~/Documents/Work and Career/LDP/Working Group/data_missing_richness.csv")
#Read in biotime metadata
metadata<-read.csv("~/Documents/Work and Career/LDP/Working Group/BioTIMEMetadata_24_06_2021.csv")

#remove marine and aquatic (keep only terr)
biotime.t<-biotime[which(biotime$STUDY_ID%in%metadata$STUDY_ID[which(metadata$REALM=="Terrestrial")]),]

#Filter biotime to only studies 39 and 221 (add 225 if using full data)
studies <- c("39", "221")
b_221_39 <- biotime.t %>%
  filter(STUDY_ID %in% studies)
#Get list of unique species names 
unique(b_221_39$GENUS_SPECIES)

#Correct the genus and species names for study 39, as above
name_mapping <- list(
  "SALIBEBB" = "Salix bebbiana",
  "POPUTRE2" = "Populus tremuloides",
  "SALIINTE" = "Salix interior",
  "POPUTRE" = "Populus tremuloides",
  "POPUTRE3" = "Populus tremuloides",
  "PICEGLA3" = "Picea glauca",
  "POPUBAL3" = "Populus balsamifera",
  "VIBUEDUL" = "Viburnum edule",
  "SALIMONT" = "Salix monticola",
  "ROSAACIC" = "Rosa acicularis",
  "BETUPAP3" = "Betula papyrifera",
  "BETUPAP2" = "Betula papyrifera",
  "RUBUIDAE" = "Rubus idaeus",
  "ALNUTENU" = "Alnus tenuifolia",
  "POPUBAL2" = "Populus balsamifera",
  "PICEGLA" = "Picea glauca",
  "SALINOVA" = "Salix novae-angliae",
  "SALISPE2" = "Salix petiolaris",
  "BETUPAP" = "Betula papyrifera",
  "PICEMAR3" = "Picea mariana",
  "ALNUCRIS" = "Alnus crispa",
  "RUBUSPEC" = "Rubus species",
  "SALIALAX" = "Salix alaxensis",
  "SALIGLAU" = "Salix glauca",
  "VIBUED" = "Viburnum edule",
  "SALILASI" = "Salix lasiandra",
  "PICEMAR2" = "Picea mariana",
  "LARILAR2" = "Larix laricina",
  "SALIPLAN" = "Salix planifolia",
  "SALIBRAC" = "Salix brachycarpa",
  "POPUTRE1" = "Populus tremuloides",
  "POPUBAL" = "Populus balsamifera",
  "LARILAR3" = "Larix laricina",
  "SALIARBU" = "Salix arbusculoides",
  "SALISCOU" = "Salix scouleriana",
  "PICEMAR" = "Picea mariana",
  "BETUGLAN" = "Betula glandulosa",
  "BETUPAP1" = "Betula papyrifera",
  "BETUNANA" = "Betula nana",
  "SALICOMM" = "Salix commutata",
  "SALIMYRT" = "Salix myrtillifolia",
  "PICEGLA1" = "Picea glauca",
  "ROSAASCI" = "Rosa acicularis",
  "BETUNAN" = "Betula nana",
  "BELNEO" = "Belangera neo-mexicana",
  "PICESPEC" = "Picea species",
  "MYRIGALE" = "Myrica gale",
  "PICEGLA2" = "Picea glauca",
  "SALINNTE" = "Salix interior",
  "SHEPCANA" = "Shepherdia canadensis",
  "SALBEBB" = "Salix bebbiana",
  "SALINIPH" = "Salix niphoclada"
)

# Function to replace unique names with correct genus and species names
convert_names <- function(x) {
  for (key in names(name_mapping)) {
    x <- gsub(key, name_mapping[[key]], x)
  }
  # Remove numbers from strings
  x <- gsub("[0-9]", "", x)
  # Remove extra spaces at the end of the string
  x <- trimws(x)
  return(x)
}

# Apply convert_names to studies
b_221_39$GENUS_SPECIES <- convert_names(b_221_39$GENUS_SPECIES)
#Check new list of unique names 
unique(b_221_39$GENUS_SPECIES)

#Convert bird names 
bird_mapping <- list(
  "Catharus guttatus" = "Catharus guttatus",
  "Catharus ustulatus" = "Catharus ustulatus",
  "Dendroica caerulescens" = "Setophaga caerulescens",
  "Catharus fuscescens" = "Catharus fuscescens",
  "Dendroica fusca" = "Setophaga fusca",
  "Empidonax minimus" = "Empidonax minimus",
  "Hylocichla mustelina" = "Catharus ustulatus",
  "Junco hyemalis" = "Junco hyemalis",
  "Pheucticus ludovicianus" = "Pheucticus ludovicianus",
  "Picoides pubescens" = "Picoides pubescens",
  "Picoides villosus" = "Picoides villosus",
  "Piranga olivacea" = "Piranga olivacea",
  "Seiurus aurocapilla" = "Seiurus aurocapilla",
  "Setophaga ruticilla" = "Setophaga ruticilla",
  "Setophaga virens" = "Setophaga virens",
  "Sitta carolinensis" = "Sitta carolinensis",
  "Sphyrapicus varius" = "Sphyrapicus varius",
  "Troglodytes hiemalis" = "Troglodytes hiemalis",
  "Vireo olivaceus" = "Vireo olivaceus",
  "Vireo philadelphicus" = "Vireo philadelphicus",
  "Vireo spec" = "Vireo NA",
  "Poecile atrcapillus" = "Poecile atricapillus",
  "Yellow-throated Warbler" = "Setophaga dominica",
  "Setophaga caerulescens" = "Setophaga caerulescens",
  "Setophaga fusca" = "Setophaga fusca",
  "Black-capped chickadee" = "Poecile atricapillus",
  "Blue jay" = "Cyanocitta cristata",
  "Brown creeper" = "Certhia americana",
  "Dark-eyed junco" = "Junco hyemalis",
  "Downy woodpecker" = "Dryobates pubescens",
  "Hairy woodpecker" = "Dryobates villosus",
  "Ovenbird" = "Seiurus aurocapilla",
  "Red-eyed vireo" = "Vireo olivaceus",
  "Rose-breasted grosbeak" = "Pheucticus ludovicianus",
  "White-breasted nuthatch" = "Sitta carolinensis",
  "Winter wren" = "Troglodytes hiemalis",
  "Wood thrush" = "Hylocichla mustelina",
  "Dryocopus pileatus" = "Dryocopus pileatus",
  "Hermit Thrush" = "Catharus guttatus",
  "Swainsons Thrush" = "Catharus ustulatus",
  "Yellow-bellied Sapsucker" = "Sphyrapicus varius",
  "Archilochus colubris" = "Archilochus colubris",
  "Contopus virens" = "Contopus virens",
  "Solitary Vireo" = "Vireo solitarius",
  "Bonasa umbellus" = "Bonasa umbellus",
  "Red-breasted nuthatch" = "Sitta canadensis",
  "Dendroica coronata" = "Setophaga coronata",
  "Certhia americana" = "Certhia americana",
  "Cyanocitta cristata" = "Cyanocitta cristata",
  "Carpodacus purpureus" = "Carpodacus purpureus",
  "Turdus migratorius" = "Turdus migratorius",
  "Mniotilta varia" = "Mniotilta varia"
)

convert_bird_names <- function(x) {
  for (key in names(bird_mapping)) {
    x <- gsub(key, bird_mapping[[key]], x)
  }
  return(x)
}

# Apply function to studies 39 and 221 to fix names 
b_221_39$GENUS_SPECIES <- convert_bird_names(b_221_39$GENUS_SPECIES)
unique(b_221_39$GENUS_SPECIES)

#Make new genus and species columns with corrected names 
genus_species <- strsplit(b_221_39$GENUS_SPECIES, " ")
b_221_39$GENUS <- sapply(genus_species, `[`, 1)
b_221_39$SPECIES <- sapply(genus_species, `[`, -1)
unique(b_221_39$SPECIES)

#add in sp per genus/plot/year
b<-b_221_39
b$UNIQUE_ID<-paste(b$STUDY_ID,b$PLOT,b$GENUS,b$YEAR,sep="~")
b$UNIQUE_ID_SP<-paste(b$UNIQUE_ID,b$SPECIES,sep="~")
b<-b[-which(duplicated(b$UNIQUE_ID_SP)),]

#Make table of unique IDs
data<-as.data.frame(table(b$UNIQUE_ID))
names(data)=c("UNIQUE_ID","Richness")

#Replace parts of study_plot names so data can talk to each other 
missing_data<- missing_data %>%
  mutate(STUDY_PLOT = str_replace_all(STUDY_PLOT, "Was_NA", "NA"))

#add richness cols
missing_data$Rich.T.Gn1=data$Richness[match(paste(missing_data$STUDY_PLOT,missing_data$Gn1,missing_data$YEAR.T,sep="~"),data$UNIQUE_ID)]
missing_data$Rich.T1.Gn1=data$Richness[match(paste(missing_data$STUDY_PLOT,missing_data$Gn1,missing_data$YEAR.T1,sep="~"),data$UNIQUE_ID)]
missing_data$Rich.T.Gn2=data$Richness[match(paste(missing_data$STUDY_PLOT,missing_data$Gn2,missing_data$YEAR.T,sep="~"),data$UNIQUE_ID)]
missing_data$Rich.T1.Gn2=data$Richness[match(paste(missing_data$STUDY_PLOT,missing_data$Gn2,missing_data$YEAR.T1,sep="~"),data$UNIQUE_ID)]

#add total sample size cols
ss_abun<-as.data.frame(table(b$UNIQUE_ID))
names(ss_abun)=c("UNIQUE_ID","Sample.Size")
missing_data$SS.T.Gn1=data$Richness[match(paste(missing_data$STUDY_PLOT,missing_data$Gn1,missing_data$YEAR.T,sep="~"),data$UNIQUE_ID)]
missing_data$SS.T1.Gn1=data$Richness[match(paste(missing_data$STUDY_PLOT,missing_data$Gn1,missing_data$YEAR.T1,sep="~"),data$UNIQUE_ID)]
missing_data$SS.T.Gn2=data$Richness[match(paste(missing_data$STUDY_PLOT,missing_data$Gn2,missing_data$YEAR.T,sep="~"),data$UNIQUE_ID)]
missing_data$SS.T1.Gn2=data$Richness[match(paste(missing_data$STUDY_PLOT,missing_data$Gn2,missing_data$YEAR.T1,sep="~"),data$UNIQUE_ID)]

#remove the unique ID column 
missing_data <- missing_data %>%
  select(-UNIQUE_ID)

#save - only for studies 221 and 39
#If running, update save location 
#write.csv(missing_data, "~/Documents/Work and Career/LDP/Working Group/data_missing_richness_ENB_032124.csv")


missing_data_check <- read.csv("~/Documents/Work and Career/LDP/Working Group/data_missing_richness.csv")

##### IE CODE - adding metadata into results.abundance #####
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

#save results.abundance with updated 
#write.csv(results.abundance,"~/Documents/Work and Career/LDP/Working Group/results.abundance_221_39.csv")
#results<-read.csv("/Users/isaaceckert/Desktop/within.study.updated.data.csv")


table(results$METRIC)




