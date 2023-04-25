#calculate portion change for species pairs 

# # LIBRARIES # #
library(tidyverse) # a suite of tidy packages useful for data manipulation
library("dplyr") # data manipulation
library("tibble") # help to create simple dataframes
library("readr") # flexible package for reading in data files
library("ggplot2") # data visualization
library("magrittr") # aids in syntax and piping

rm(list=ls()) 

# # INPUT FILES # #
load("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/tidy/collated_pairs.RData")#collated pairs of overlapping studies
load("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/prep_biotime/bio_pairs_10km.RData") #metadata of overlapping studies

## So let's explore a particular pair of time series with lots of overlapping data
pair=164

pair_1_ID <- bio.pairs$ID.1[pair]
pair_2_ID <- bio.pairs$ID.2[pair]

timeseries_1 <- collated.pairs %>% dplyr::filter(ID == pair_1_ID)
timeseries_2 <- collated.pairs %>% dplyr::filter(ID == pair_2_ID)

years_overlap <- unique(timeseries_1$YEAR)[unique(timeseries_1$YEAR) %in% unique(timeseries_2$YEAR)] %>% sort()

## Let's make a list of the species in these two timeseries
timeseries_1_species <- timeseries_1$SPECIES %>% unique()
timeseries_1_species = timeseries_1_species[!is.na(timeseries_1_species)]
timeseries_1_species_length <- sapply(timeseries_1_species, function(x){
  timeseries_1 %>%
    dplyr::filter(SPECIES == x) %>%
    dplyr::select(YEAR) %>%
    unique() %>%
    unlist() %>%
    length()
})

timeseries_2_species <- timeseries_2$SPECIES %>% unique()
timeseries_2_species = timeseries_2_species[!is.na(timeseries_2_species)]
timeseries_2_species_length <- sapply(timeseries_2_species, function(x){
  timeseries_2 %>%
    dplyr::filter(SPECIES == x) %>%
    dplyr::select(YEAR) %>%
    unique() %>%
    unlist() %>%
    length()
})

sp1 = timeseries_1_species[sample(x = which(timeseries_1_species_length > 10), size = 1)]
sp2 = timeseries_2_species[sample(x = which(timeseries_2_species_length > 10), size = 1)]

#calc por change
calc.corr<-function(sp1,sp2,lag){
  
  sp1_data <- timeseries_1 %>%
    dplyr::filter(
      SPECIES == sp1,
      YEAR %in% years_overlap) %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(
      Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
        mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
        mean(sum.allrawdata.BIOMASS, na.rm = TRUE),
        na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(YEAR)%>%
    distinct(.)
  
  sp2_data <- timeseries_2 %>%
    dplyr::filter(
      SPECIES == sp2,
      YEAR %in% years_overlap) %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarize(
      Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
        mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
        mean(sum.allrawdata.BIOMASS, na.rm = TRUE),
        na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(YEAR)
  
  if (length(which(sp1_data$YEAR%in%sp2_data$YEAR))<4){print("Less than 4 overlapping years")}
  if (length(which(sp1_data$YEAR%in%sp2_data$YEAR))>=4){
  #only keep overlapping years
  sp1_data<-sp1_data[which(sp1_data$YEAR%in%sp2_data$YEAR),]
  sp2_data<-sp2_data[which(sp2_data$YEAR%in%sp1_data$YEAR),]
  
  
  N1_sp1 = sp1_data$Abundance[2:nrow(sp1_data)]
  N0_sp1 = sp1_data$Abundance[1:nrow(sp1_data)-1]
  
  N1_sp2 = sp2_data$Abundance[2:nrow(sp2_data)]
  N0_sp2 = sp2_data$Abundance[1:nrow(sp2_data)-1]
  
  log_prop_change_sp1 = log(N1_sp1/N0_sp1)
  log_prop_change_sp2 = log(N1_sp2/N0_sp2)
  
  length(log_prop_change_sp1)
  length(log_prop_change_sp2)
  
  cor.test(log_prop_change_sp1[1:(length(log_prop_change_sp1)-lag)],log_prop_change_sp2[(1+lag):length(log_prop_change_sp2)])}
  
}

mat<-as.data.frame(matrix(data=NA,ncol=4,nrow=0))
colnames(mat)=c("Species1","Species2","Corr","p")
for (i in 1:length(timeseries_1_species)){

  sp1=timeseries_1_species[i]
  
  for (x in 1:length(timeseries_2_species)){

    sp2=timeseries_2_species[x]
    
    c=calc.corr(sp1,sp2,lag=0)
    
    if (is.null(c)){
      
      mat<-rbind(mat,data.frame("Species1"=sp1,
                                "Species2"=sp2,
                                "Corr"=NA,
                                "p"=NA))
    }
    
    else {
      mat<-rbind(mat,data.frame("Species1"=sp1,
                                "Species2"=sp2,
                                "Corr"=c$estimate,
                                "p"=c$p.value))}

  }
  
  
  print(paste("DONE:",i,"out of",length(timeseries_1_species)))
  

  
}

mat.sig<-mat[which(mat$p<=0.05),]

collated.pairs.wking<-collated.pairs %>%
  group_by(ID,LATITUDE,LONGITUDE,SAMPLE_DESC,GENUS) %>%
  summarize(mean_abun=mean(sum.allrawdata.ABUNDANCE),
            mean_bio=mean(sum.allrawdata.BIOMASS))










