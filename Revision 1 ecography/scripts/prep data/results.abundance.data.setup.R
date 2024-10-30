# Biotime cleaning, aggregating, and change in abundance calculations
# Date created: Jan 24 2024
# Date updated: 19 Sept 2024

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


######Note: Studies 39 and 221 contain errors, but were removed from analyses

##### Back to IE code ####

#any false 0s?
table(biotime.t$sum.allrawdata.ABUNDANCE[which(
  biotime.t$STUDY_ID%in%
    metadata$STUDY_ID[which(is.na(metadata$ABUNDANCE_TYPE))])])
#table(biotime.t$sum.allrawdata.BIOMASS[which(biotime.t$STUDY_ID%in%metadata$STUDY_ID[which(is.na(metadata$BIOMASS_TYPE))])])

#replace 0s with NAs if the study didnt record how they measured abundance or biomass...
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

#Note: save locations should be updated if running code on own device 
#save results.abundance 
#write.csv(results.abundance,"Revision 1 ecography/output/prep_data/results.abundance.csv")

table(results.abundance$METRIC)




