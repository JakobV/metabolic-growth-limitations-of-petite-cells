## The metabolic growth limitations of petite cells lacking the mitochondrial genome ##
## Vowinckel et al., 2021 #############################################################

## Analysis and Plotting scripts

## Environment----
source("Scripts/general_settings_v02.R")

library(gridExtra)
library(ggplot2)
library(png)
library(reshape)
library(RColorBrewer)
library(gplots)
library(grid)
library(dplyr)
library(stringr)
library(tidyr)
library(outliers)

Figures <- list()
Dimensions <- list()

## Data Analysis----

### Yeast Lifespan
{
  #lifespan experiment
  data <- read.csv("Data/150414_time_course_replicates.tsv", sep="\t")
  time <- data[which(data$Experiment == "cls"),]
  time$name <- factor(time$name)
  #add EC
  t <- cbind(time, ID=paste(time$StrainID, time$timestamp))
  t <- cast(t[which(t$name %in% levels(t$name)[c(8,6,7)]),], ID ~ name)
  t <- t[,c(1,4,2,3)]
  t <- cbind(t, EC=(t[,2]+t[,3]*0.5)/(t[,2]+t[,3]+t[,4]))
  t <- melt(t, id="ID")
  t <- cbind(t, t(as.data.frame(strsplit(as.character(t$ID), " "))))
  t <- t[which(t$variable == "EC"),-1]
  names(t) <- c("name", "value", "Strain", "Replicate", "timestamp")
  t <- cbind(t, StrainID=paste(t$Strain, t$Replicate))
  
  u <- time[which(time$name == "ATP [cytoplasm]"),]
  u <- merge(u, t, by=c("Strain", "StrainID", "timestamp", "Replicate"), all=T)[,c(-7,-8)]
  names(u)[10:11] <- c("name", "value")
  time <- rbind(time, u[,c(5,2,1,4,6,3,10,11,7,8,9)])
  
  #correct units
  #glucose: mM
  time$value[which(time$name == levels(time$name)[15])] <- time$value[which(time$name == levels(time$name)[15])]/1000
  #energy metabolites: mM
  time$value[which(time$name %in% levels(time$name)[6:8])] <- time$value[which(time$name %in% levels(time$name)[6:8])]/1000
  
  
  # summary
  {
  time_summary <- cbind(aggregate(time$value, by=list(time$Experiment, time$Strain, time$time, time$timestamp, time$name, time$unit), FUN=mean, na.rm=T), aggregate(time$value, by=list(time$Experiment, time$Strain, time$time, time$timestamp, time$name, time$unit), FUN=sd, na.rm=T)[7])
  names(time_summary) <- c("Experiment","Strain", "time", "timestamp", "name", "unit","value_mean", "value_sd")
  
  #correct timestamp unit
  time_summary$time <- as.POSIXlt(time_summary$time)
  time$time <- as.POSIXlt(time$time)
  time$timestamp_ <- time$timestamp/60/60
  write.table(time, "Data/time_course_viability_OD_EC.tsv", sep="\t", row.names=F)
  
  t <- aggregate(time_summary$timestamp, by=list(time_summary$name), FUN=max)
  time_summary$timestamp[which(time_summary$name %in% as.character(t[which(t$x > 200),1]))] <- time_summary$timestamp[which(time_summary$name %in% as.character(t[which(t$x > 200),1]))]/60/60
  
  #correct outliers
  time_summary$value_mean[which(time_summary$name == "EC" & round(time_summary$timestamp) == 39 & time_summary$Strain == "rho0_evo")] <- time_summary$value_mean[which(time_summary$name == "EC" & round(time_summary$timestamp) == 27 & time_summary$Strain == "rho0_evo")]
  time_summary$value_mean[which(time_summary$name == "EC" & round(time_summary$timestamp) == 39 & time_summary$Strain == "rho0")] <- time_summary$value_mean[which(time_summary$name == "EC" & round(time_summary$timestamp) == 27 & time_summary$Strain == "rho0")]
  time_summary$value_mean[which(time_summary$name == "biomass" & round(time_summary$timestamp) == 24 & time_summary$Strain == "rho+")] <- NA
  time_summary$value_mean[which(time_summary$name == "biomass" & round(time_summary$timestamp) == 54 & time_summary$Strain == "rho+")] <- NA
  time_summary$value_mean[which(time_summary$name == "biomass" & round(time_summary$timestamp) == 54 & time_summary$Strain == "rho0")] <- NA
  time_summary$value_mean[which(time_summary$name == "biomass" & round(time_summary$timestamp) == 54 & time_summary$Strain == "rho0_evo")] <- NA
  time_summary$value_mean[which(time_summary$name == "biomass" & round(time_summary$timestamp) == 48 & time_summary$Strain == "rho0")] <- NA
  time_summary$value_mean[which(time_summary$name == "biomass" & round(time_summary$timestamp) == 48 & time_summary$Strain == "rho0_evo")] <- NA
  time_summary$value_mean[which(time_summary$name == "biomass" & round(time_summary$timestamp) == 42 & time_summary$Strain == "rho0")] <- NA
  time_summary$value_mean[which(time_summary$name == "biomass" & round(time_summary$timestamp) == 42 & time_summary$Strain == "rho0_evo")] <- NA
  
  #correct different time scales
  time_summary <- cbind(time_summary, time_sec = as.numeric(time_summary$time))
  time_summary$time_sec <- time_summary$time_sec-min(time_summary$time_sec)
  time_summary$time_sec <- time_summary$time_sec/60/60
  
  t <- time_summary[which(time_summary$name %in% levels(time_summary$name)[c(8,9,56,57, 6,7,8,15)]),]

  write.table(t, "Data/time_course_viability_OD_EC_summary.tsv", sep="\t", row.names=F)
  }
  
  
}

### Metabolite Abundances
{
## get PPP data
code <- read.xlsx("Data/20131212_sample_code.xlsx", sheetIndex=1)
cdata <- read.xlsx("Data/20131217_TCAPPP_JH.xlsx", sheetIndex=2)

cdata <- cdata %>%
  dplyr::rename(Sample = `NA.`) %>%
  tidyr::pivot_longer(cols = -Sample, names_to = "Compound") %>%
  mutate(experiment = "PPP")

## get AA data

t <- read.csv("Data/2014-01-17_RESULT_concentration_in_µM.csv")

t <- t %>%
  dplyr::rename(Sample = sample) %>%
  tidyr::pivot_longer(cols = -Sample, names_to = "Compound") %>%
  mutate(experiment = "AA")

## combine data
cdata <- rbind(cdata, t)

cdata$Sample <- sub("Sample_", "", cdata$Sample)
cdata <- merge(cdata, code, by="Sample", all.x=T)

cdata <- cbind(cdata, ID = cdata$Strain)
for(i in c(" c1", " c2", " c3")){
  cdata$ID <- sub(i, "", cdata$ID)
}
cdata$Strain <- NULL

# clean up data
cdata$value[which(cdata$Compound == "Succinate")][c(3,6,7,14,15)] <- NA
cdata$value[which(cdata$Compound == "lysine" & cdata$ID == "rho0_313")][2] <- NA
cdata$value[which(cdata$Compound == "tyrosine" & cdata$ID == "rho0_313")][2] <- NA
cdata$value[which(cdata$Compound == "Fumarate" & cdata$ID == "rho0_ATP3_T911A")][2] <- NA
cdata$value[which(cdata$Compound == "Fumarate" & cdata$ID == "rho0_313")][2] <- NA
cdata$value[which(cdata$Compound == "Fumarate" & cdata$ID == "rho+_313")][3] <- NA
cdata$value[which(cdata$Compound == "Succinate" & cdata$ID == "rho+_313")][1] <- NA
cdata$value[which(cdata$Compound == "Succinate" & cdata$ID == "rho0_313")][3] <- NA
cdata$value[which(cdata$Compound == "X2.KG" & cdata$ID == "rho0_ATP3_T911A")][2] <- NA
cdata$value[which(cdata$Compound == "X2.KG" & cdata$ID == "rho+_313")][1] <- NA
cdata$experiment[cdata$Compound %in% c("Glucose", "G6P", "G3P", "F6P", "DHAP", "2.3.PG", "1.6.FP", "PEP")] <- "CCM"
cdata$experiment[cdata$Compound %in% c("Pyruvate", "Succinate", "Malate", "2.KG", "Fumarate", "Oxalacetate", "cis.Aconitate", "Citrate", "Acetyl.CoA")] <- "TCA"

write.table(cdata, "Data/20210816_Metabolite_Abundances.tsv", sep = "\t", row.names=F)
}

### GSH/GSSG Levels
{
  data <- read.csv("Data/140407_replicates.csv", stringsAsFactors = T)
  data <- data[,c(7,8,9,12,13,14,15,16,17,18,19)]
  data <- melt(data, id=c("dilution","strain", "genotype", "date", "replicate"))
  data <- cbind(data, t(as.data.frame(strsplit(as.character(data$variable), "_"))))  
  names(data)[8:9] <- c("compound", "transition")
  data <- aggregate(data$value, by=list(data$dilution, data$strain, data$genotype, data$date, data$replicate, data$compound), FUN=sum, na.rm=T)
  names(data)[1:7] <- c("dilution", "strain", "genotype", "date", "replicate", "compound", "sum")
  data <- data[which(data$date == "05/03/2014"),]
  data <- data[!(data$compound == "GSSG" & data$dilution == 20),]
  data <- data[!(data$compound == "GSH" & data$dilution == 1),]
  data <- cbind(data, condition=paste(data$strain, data$genotype, data$replicate))
  data$condition <- factor(data$condition)
  levels(data$condition)[1:3] <- sub(" +", "_+", levels(data$condition)[1:3], fixed=T)
  
  t <- cast(data, condition ~ compound, value="sum")
  t <- cbind(t, GSH.GSSG=t$GSH/t$GSSG)
  
  t <- melt(t)
  summary <- cbind(t, t(as.data.frame(strsplit(as.character(t$condition), " "))))
  summary <- cbind(summary, strainID=paste(summary[,4], summary[,5]))
  
  #p values
  pvalues <- data.frame()
  for(i in unique(summary$variable)){
    for(j in unique(summary$strainID)){
      pvalues <- rbind(pvalues, data.frame(compound=i, condition = j, pvalue=t.test(summary$value[which(summary$strainID == j & summary$variable == i)], summary$value[which(summary$strainID == unique(summary$strainID)[2] & summary$variable == i)])$p.value))
    }
  }
  
  sdata <- summary
  summary <- cbind(aggregate(summary$value, by=list(summary$variable, summary[,4], summary[,5]), FUN=mean), aggregate(summary$value, by=list(summary$variable, summary[,4], summary[,5]), FUN=sd)[4])
  names(summary) <- c("compound", "strain", "genotype", "mean", "sd")
  summary <- cbind(summary, condition=paste(summary$strain, summary$genotype))
  summary$mean[which(summary$compound %in% levels(summary$compound)[1:2])] <- summary$mean[which(summary$compound %in% levels(summary$compound)[1:2])]/1000
  summary$sd[which(summary$compound %in% levels(summary$compound)[1:2])] <- summary$sd[which(summary$compound %in% levels(summary$compound)[1:2])]/1000
  summary <- merge(summary, pvalues, by=c("condition", "compound"), all.x=T)
  summary <- cbind(summary, signif=NA)
  summary$signif[which(summary$pvalue < 0.05)] <- "*"
  
  sdata$value[sdata$variable %in% c("GSH", "GSSG")] <- sdata$value[sdata$variable %in% c("GSH", "GSSG")]/1000
  
  ggplot(summary, aes(x = condition, y = mean))+
    geom_bar(stat="identity")+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
    facet_wrap(~ compound, scales="free")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  write.table(summary, "Data/140407_replicates_summary.csv", sep=",")
  write.table(sdata, "Data/140407_replicates_data.csv", sep=",")
  
  
}

### DHE Fluorescence
{
  #split channels & re-format images to 8-bit TIF
  #analyse images using CP pipeline 151110_DHE.cpproj

  cells <- read.csv("Data/151110_DHE_cells.csv", stringsAsFactors = T)
  DHEcells <- cells[,c(10,8,12)]
  cells <- cells[,c(10,8,11)]
  names(DHEcells) <- unlist(DHEcells[1,])
  names(cells) <- unlist(cells[1,])
  DHEcells <- DHEcells[-1,]
  cells <- cells[-1,]
  DHEcells$Number_Object_Number <- as.numeric(as.character(DHEcells$Number_Object_Number))
  DHEcells <- DHEcells[!is.na(DHEcells$Number_Object_Number),]
  cells$Number_Object_Number <- as.numeric(as.character(cells$Number_Object_Number))
  
  summary <- aggregate(cells$Number_Object_Number, by=list(cells$Metadata_Strain, cells$Metadata_Image), FUN=length)
  names(summary) <- c("Strain", "Image", "Cells")
  t <- aggregate(DHEcells$Number_Object_Number, by=list(DHEcells$Metadata_Strain, DHEcells$Metadata_Image), FUN=length)
  names(t) <- c("Strain", "Image", "DHE_Cells")
  summary <- merge(summary, t, by=c("Strain", "Image"))
  summary <- cbind(summary, DHE_percent=100/summary$Cells*summary$DHE_Cells)
  
  ggplot(summary, aes(x = Strain, y = DHE_percent, group=Strain))+
    geom_boxplot()+
    theme_classic()
  
  summary$DHE_percent[which(summary$Strain == "wttBOOH")][6] <- NA
  
  #ttest
  out <- data.frame()
  for(i in unique(summary$Strain)){
    out <- rbind(out, data.frame(Strain=i, pvalue=t.test(summary$DHE_percent[which(summary$Strain == i)], summary$DHE_percent[which(summary$Strain == "wt")])$p.value))
  }
  
  summary2 <- cbind(aggregate(summary$DHE_percent, by=list(summary$Strain), FUN=mean, na.rm=T), aggregate(summary$DHE_percent, by=list(summary$Strain), FUN=sd, na.rm=T)[2])
  names(summary2) <- c("Strain", "mean", "sd")
  summary2 <- merge(summary2, out, by="Strain")
  summary2$Strain <- factor(summary2$Strain, levels = levels(summary2$Strain)[c(4,2,3,5)])
  levels(summary2$Strain)[4] <- "wt+tertBOOH"
  summary2 <- cbind(summary2, signif=NA)
  summary2$signif[which(summary2$pvalue < 0.01)] <- "*"

  #write to xls file
  wb <- createWorkbook()
  
  sheet1 <- createSheet(wb, sheetName="raw data")
  sheet2 <- createSheet(wb, sheetName="percentages")
  sheet3 <- createSheet(wb, sheetName="summary")
  addDataFrame(cells, sheet1)
  addDataFrame(summary, sheet2)
  addDataFrame(summary2, sheet3)
  
  saveWorkbook(wb, "Data/151110_DHE_Summary.xlsx")
}

### PTS1-GFP Microscopy
{
  #read CellProfiler output data
  #Rho0
  data <- read.csv("Data/Rho0/Cells.csv")
  data <- data[,-4]
  names(data)[2] <- "Cells"
  cells <- read.csv("Data/Rho0/FilteredCells.csv")
  cells <- cells[,c(1,8,2,3,5,6,7)]
  names(cells)[2:3] <- c("Cells", "FilteredCells")
  data <- merge(data, cells, by=c("ImageNumber", "Cells", "Metadata_cell", "Metadata_genotype"), all=T)
  
  #WT
  t <- read.csv("Data/WT/Cells.csv")
  t <- t[,-4]
  names(t)[2] <- "Cells"
  cells <- read.csv("Data/WT/FilteredCells.csv")
  cells <- cells[,c(1,8,2,3,5,6,7)]
  names(cells)[2:3] <- c("Cells", "FilteredCells")
  t <- merge(t, cells, by=c("ImageNumber", "Cells", "Metadata_cell", "Metadata_genotype"), all=T)
  
  data <- rbind(data, t)
  
  #remove filtered cells
  data <- data[which(data$Children_FilteredCells_Count > 0),]
  
  #remove cells with FormFactor <= 0.6 & Area > 10000
  data <- data[which(data$AreaShape_FormFactor > 0.6 & data$AreaShape_Area < 10000),]
  
  summary <- cbind(aggregate(data$Children_PXS_Count, by=list(data$Metadata_genotype), FUN=mean), aggregate(data$Children_PXS_Count, by=list(data$Metadata_genotype), FUN=median)[2],aggregate(data$Children_PXS_Count, by=list(data$Metadata_genotype), FUN=sd)[2], aggregate(data$Children_PXS_Count, by=list(data$Metadata_genotype), FUN=se)[2])
  names(summary) <- c("Metadata_genotype", "mean", "median", "stdev", "sem")
  summary$Metadata_genotype <- factor(summary$Metadata_genotype, levels=c("WT_313", "Rho_313", "Rho_ATP3T911A", "Rho_ATP3G919C"))
  data$Metadata_genotype <- factor(data$Metadata_genotype, levels=c("WT_313", "Rho_313", "Rho_ATP3T911A", "Rho_ATP3G919C"))
  
  #calculate n and significance
  
  n <- aggregate(data$FilteredCells, by=list(data$Metadata_genotype), FUN=length)
  names(n) <- c("Metadata_genotype", "n")
  t <- t.test(data$Children_PXS_Count[which(data$Metadata_genotype == "WT_313")], data$Children_PXS_Count[which(data$Metadata_genotype == "Rho_313")])[3]
  ttest <- data.frame(group1= "WT_313", group2="Rho_313", p=t)
  t <- data.frame(group1= "WT_313", group2="Rho_ATP3T911A", p=t.test(data$Children_PXS_Count[which(data$Metadata_genotype == "WT_313")], data$Children_PXS_Count[which(data$Metadata_genotype == "Rho_ATP3T911A")])[3])
  ttest <- rbind(ttest, t)
  t <- data.frame(group1= "WT_313", group2="WT_313", p=t.test(data$Children_PXS_Count[which(data$Metadata_genotype == "WT_313")], data$Children_PXS_Count[which(data$Metadata_genotype == "WT_313")])[3])
  ttest <- rbind(ttest, t)

  names(ttest)[1:2] <- c("Reference", "Metadata_genotype")
  summary <- merge(summary, ttest[,2:3], by="Metadata_genotype")
  summary <- merge(summary, n, by="Metadata_genotype")
  write.table(summary, "Data/140227_PTS1_results.csv", sep=",", quote=F, row.names=F)
  write.table(data, "Data/140227_PTS1_results_raw.csv", sep=",", quote=F, row.names=F)
}

## Proteomics data (General)
{
  pdata <- data.table::fread("Data/201229_YSBN11_MR_150429_Report.xls")
  pdata <- pdata %>% filter(!str_detect(R.FileName, "-Fe"))
  
  pdata$SampleID <- sub("11_", "", stringr::str_extract(string = pdata$R.FileName, pattern = "11_.+"))
  
  pdata$Genotype <- NA
  pdata$Genotype[grepl("wt", pdata$SampleID, ignore.case = T)] <- "WT"
  pdata$Genotype[grepl("rho", pdata$SampleID, ignore.case = T)] <- "Rho0"
  
  pdata$Plasmid <- NA
  pdata$Plasmid[grepl("T911A", pdata$SampleID, ignore.case = T)] <- "ATP3-T911A"
  pdata$Plasmid[grepl("G919C", pdata$SampleID, ignore.case = T)] <- "ATP3-G919C"
  pdata$Plasmid[is.na(pdata$Plasmid) & grepl("ATP3", pdata$SampleID, ignore.case = T)] <- "ATP3"
  pdata$Plasmid[is.na(pdata$Plasmid)] <- "Control"
  
  pdata$Replicate <- str_extract(string = pdata$SampleID, pattern = "\\d$")
  pdata$Experiment <- "Fe"
  
  pdata %>%
    distinct(R.FileName, .keep_all=T) %>%
    select(SampleID, Genotype, Plasmid, Replicate, Experiment)
  
  pdata <- pdata %>%
    distinct(R.FileName, PG.ProteinAccessions, .keep_all = T) %>%
    mutate(StrainID = paste(Genotype, Plasmid)) %>%
    filter(StrainID %in% c("WT Control", "Rho0 Control", "Rho0 ATP3-T911A")) %>%
    mutate(StrainID = factor(StrainID, levels = c("WT Control", "Rho0 Control", "Rho0 ATP3-T911A")))
  
  pdata$ID <- pdata$StrainID
  levels(pdata$ID) <- c("rho+_313", "rho0_313", "rho0_ATP3_T911A")
  
  odata <- pdata %>%
    filter(Experiment == "Fe") %>%
    mutate(experiment = "proteomics") %>%
    select(R.FileName, PG.Genes, PG.Quantity, experiment, ID) %>%
    dplyr::rename(Sample = R.FileName, Compound = PG.Genes, value = PG.Quantity)
  
  data.table::fwrite(odata, "Data/20210816_Proteomics_General.tsv", sep="\t", row.names=F)
}

## Proteomics data (Iron)
{
  pdata <- data.table::fread("Data/201229_YSBN11_MR_150429_Report.xls")
  pdata$SampleID <- sub("11_", "", stringr::str_extract(string = pdata$R.FileName, pattern = "11_.+"))
  
  pdata$Genotype <- NA
  pdata$Genotype[grepl("wt", pdata$SampleID, ignore.case = T)] <- "WT"
  pdata$Genotype[grepl("rho", pdata$SampleID, ignore.case = T)] <- "Rho0"
  
  pdata$Plasmid <- NA
  pdata$Plasmid[grepl("T911A", pdata$SampleID, ignore.case = T)] <- "ATP3-T911A"
  pdata$Plasmid[grepl("G919C", pdata$SampleID, ignore.case = T)] <- "ATP3-G919C"
  pdata$Plasmid[is.na(pdata$Plasmid) & grepl("ATP3", pdata$SampleID, ignore.case = T)] <- "ATP3"
  pdata$Plasmid[is.na(pdata$Plasmid)] <- "Control"
  
  pdata$Replicate <- str_extract(string = pdata$SampleID, pattern = "\\d$")
  
  pdata$Iron <- NA
  pdata$Iron[grepl("-Fe", pdata$SampleID, ignore.case = T, fixed=T)] <- "- Iron"
  pdata$Iron[grepl("+Fe", pdata$SampleID, ignore.case = T, fixed=T)] <- "+ Iron"
  
  pdata$Experiment <- "Fe"
  
  pdata %>%
    distinct(R.FileName, .keep_all=T) %>%
    select(SampleID, Genotype, Plasmid, Replicate, Experiment, Iron)
  
  pdata <- pdata %>%
    distinct(R.FileName, PG.ProteinAccessions, .keep_all = T) %>%
    mutate(StrainID = paste(Genotype, Plasmid, Iron)) 
  
  pdata <- pdata %>%
    mutate(StrainID = factor(StrainID, levels = c("WT Control - Iron", "Rho0 Control - Iron", "Rho0 ATP3-T911A - Iron", "WT Control + Iron", "Rho0 Control + Iron", "Rho0 ATP3-T911A + Iron")))
  
  pdata$ID <- pdata$StrainID
  
  odata <- pdata %>%
    filter(Experiment == "Fe") %>%
    mutate(experiment = "proteomics") %>%
    select(R.FileName, PG.Genes, PG.Quantity, experiment, ID) %>%
    dplyr::rename(Sample = R.FileName, Compound = PG.Genes, value = PG.Quantity)
  odata$sgdName <- paste0(str_to_title(tolower(odata$Compound)), "p")
  cdata <- odata
  
  # fold change
  normdata <- cdata %>%
    group_by(Compound, ID, experiment) %>%
    summarise(mean = mean(value, na.rm=T), sd = sd(value, na.rm=T)) %>%
    dplyr::rename(Strain = ID, Experiment = experiment) %>%
    mutate(Compound = sub("^X", "", Compound)) %>%
    ungroup() %>%
    filter(Strain == "WT Control + Iron") %>%
    select(Compound, mean) %>%
    dplyr::rename(norm = mean)
  
  fcdata <- cdata %>%
    dplyr::rename(Strain = ID, Experiment = experiment) %>%
    mutate(Compound = sub("^X", "", Compound)) %>%
    left_join(normdata) %>%
    mutate(foldchange = value/norm) %>%
    select(-norm)
  
  fcsummary <- fcdata %>%
    group_by(Compound, Strain, Experiment) %>%
    summarise(mean = mean(foldchange, na.rm=T), sd = sd(foldchange, na.rm=T))
  
  # zscore
  zdata <- cdata %>%
    dplyr::group_by(Compound, experiment) %>%
    mutate(zscore = scale(value)) %>%
    dplyr::rename(Strain = ID, Experiment = experiment) %>%
    mutate(Compound = sub("^X", "", Compound)) %>%
    ungroup()
  
  zsummary <- zdata %>%
    group_by(Compound, Strain, Experiment) %>%
    summarise(mean = mean(zscore, na.rm=T), sd = sd(zscore, na.rm=T))
  
  
  # significance
  signif <- data.frame()
  for(i in unique(fcsummary$Compound)){
    t <- fcdata[fcdata$Compound == i,]
    for(j in unique(t$Strain)){
      signif <- rbind(signif, data.frame(Compound = i, Strain = j, p.value = t.test(t$foldchange[t$Strain == j], t$foldchange[t$Strain == "WT Control + Iron"], var.equal = T)$p.value))
      
    }
  }
  signif$signif <- F
  signif$signif[signif$p.value < 0.01] <- T
  
  write.table(signif, "Data/20210816_IronData_Signif.tsv", sep="\t", row.names=F)
  write.table(fcsummary, "Data/20210816_IronData_FoldChange_Summarized.tsv", sep="\t", row.names=F)
  write.table(fcdata, "Data/20210816_IronData_FoldChange_Data.tsv", sep="\t", row.names=F)
  write.table(zsummary, "Data/20210816_IronData_ZScores_Summarized.tsv", sep="\t", row.names=F)
  write.table(zdata, "Data/20210816_IronData_ZScores_Data.tsv", sep="\t", row.names=F)
  write.table(pdata, "Data/20210816_IronData_Raw_Data.tsv", sep="\t", row.names=F)
  
}

## Proteomics data (aco)
{
  pdata <- data.table::fread("Data/20210111_211757_210111_BY4741_aco1_MR_Report.xls")
  
  pdata$SampleID <- sub("9_aco1_", "", stringr::str_extract(string = pdata$R.FileName, pattern = "9_aco1_.+"))
  
  pdata$Genotype <- NA
  pdata$Genotype[grepl("rho", pdata$SampleID, ignore.case = T)] <- "Rho0"
  pdata$Genotype[grepl("wt", pdata$SampleID, ignore.case = T)] <- "WT"
  pdata$Genotype[grepl("aco1", pdata$SampleID, ignore.case = T)] <- "aco1"
  
  pdata$Plasmid <- NA
  pdata$Plasmid[grepl("T911A", pdata$SampleID, ignore.case = T)] <- "ATP3-T911A"
  pdata$Plasmid[is.na(pdata$Plasmid)] <- "Control"
  
  pdata$Replicate <- str_extract(string = pdata$SampleID, pattern = "\\d$")
  
  pdata %>%
    distinct(R.FileName, .keep_all=T) %>%
    select(SampleID, Genotype, Plasmid, Replicate)
  
  pdata <- pdata %>%
    distinct(R.FileName, PG.ProteinAccessions, .keep_all = T) %>%
    mutate(StrainID = paste(Genotype, Plasmid)) 
  
  pdata <- pdata %>%
    mutate(StrainID = factor(StrainID, levels = c("WT Control",      "WT ATP3-T911A", "Rho0 Control",    "Rho0 ATP3-T911A", "aco1 Control",    "aco1 ATP3-T911A"  )))
  
  pdata$ID <- pdata$StrainID
  #levels(pdata$ID) <- c("rho+_313", "rho0_313", "rho0_ATP3_T911A")
  
  odata <- pdata %>%
    select(R.FileName, PG.Genes, PG.Quantity, ID) %>%
    dplyr::rename(Sample = R.FileName, Compound = PG.Genes, value = PG.Quantity)
  odata$sgdName <- paste0(str_to_title(tolower(odata$Compound)), "p")
  cdata <- odata
  
  # fold change
  normdata <- cdata %>%
    group_by(Compound, ID) %>%
    summarise(mean = mean(value, na.rm=T), sd = sd(value, na.rm=T)) %>%
    dplyr::rename(Strain = ID) %>%
    mutate(Compound = sub("^X", "", Compound)) %>%
    ungroup() %>%
    filter(Strain == "WT Control") %>%
    select(Compound, mean) %>%
    dplyr::rename(norm = mean)
  
  fcdata <- cdata %>%
    dplyr::rename(Strain = ID) %>%
    mutate(Compound = sub("^X", "", Compound)) %>%
    left_join(normdata) %>%
    mutate(foldchange = value/norm) %>%
    select(-norm)
  
  fcsummary <- fcdata %>%
    group_by(Compound, Strain) %>%
    summarise(mean = mean(foldchange, na.rm=T), sd = sd(foldchange, na.rm=T))
  
  # zscore
  zdata <- cdata %>%
    dplyr::group_by(Compound) %>%
    mutate(zscore = scale(value)) %>%
    dplyr::rename(Strain = ID) %>%
    mutate(Compound = sub("^X", "", Compound)) %>%
    ungroup()
  
  zsummary <- zdata %>%
    group_by(Compound, Strain) %>%
    summarise(mean = mean(zscore, na.rm=T), sd = sd(zscore, na.rm=T))
  
  # significance
  # signif <- data.frame()
  # for(i in unique(fcsummary$Compound)){
  #   t <- fcdata[fcdata$Compound == i,]
  #   for(j in unique(t$Strain)){
  #     signif <- rbind(signif, data.frame(Compound = i, Strain = j, p.value = t.test(t$foldchange[t$Strain == j], t$foldchange[t$Strain == "WT Control + Iron"], var.equal = T)$p.value))
  #     
  #   }
  # }
  # signif$signif <- F
  # signif$signif[signif$p.value < 0.01] <- T  
  
  write.table(fcsummary, "Data/20210816_AcoData_FoldChange_Summarized.tsv", sep="\t", row.names=F)
  write.table(fcdata, "Data/20210816_AcoData_FoldChange_Data.tsv", sep="\t", row.names=F)
  write.table(zsummary, "Data/20210816_AcoData_ZScores_Summarized.tsv", sep="\t", row.names=F)
  write.table(zdata, "Data/20210816_AcoData_ZScores_Data.tsv", sep="\t", row.names=F)
  write.table(pdata, "Data/20210816_AcoData_Raw_Data.tsv", sep="\t", row.names=F)
}

## Combined proteomics and metabolite data
{
  cdata <- data.table::fread("Data/20210816_Metabolite_Abundances.tsv")
  odata <- data.table::fread("Data/20210816_Proteomics_General.tsv")
  cdata <- rbind(cdata, odata)
  
  ## calculate fold change and zscore
  {
    # fold change
    normdata <- cdata %>%
      group_by(Compound, ID, experiment) %>%
      dplyr::summarise(mean = mean(value, na.rm=T), sd = sd(value, na.rm=T)) %>%
      dplyr::rename(Strain = ID, Experiment = experiment) %>%
      mutate(Compound = sub("^X", "", Compound)) %>%
      ungroup() %>%
      filter(Strain == "rho+_313") %>%
      select(Compound, mean) %>%
      dplyr::rename(norm = mean)
    
    fcdata <- cdata %>%
      dplyr::rename(Strain = ID, Experiment = experiment) %>%
      mutate(Compound = sub("^X", "", Compound)) %>%
      left_join(normdata) %>%
      mutate(foldchange = value/norm) %>%
      select(-norm)
    
    fcsummary <- fcdata %>%
      group_by(Compound, Strain, Experiment) %>%
      dplyr::summarise(mean = mean(foldchange, na.rm=T), sd = sd(foldchange, na.rm=T))
    
    # zscore
    zdata <- cdata %>%
      dplyr::group_by(Compound, experiment) %>%
      dplyr::mutate(zscore = scale(value, center = T, scale = T)) %>%
      dplyr::rename(Strain = ID, Experiment = experiment) %>%
      dplyr::mutate(Compound = sub("^X", "", Compound)) %>%
      ungroup()
    
    zsummary <- zdata %>%
      group_by(Compound, Strain, Experiment) %>%
      dplyr::summarise(mean = mean(zscore, na.rm=T), sd = sd(zscore, na.rm=T))
    
    # significance
    signif <- data.frame()
    for(i in unique(fcsummary$Compound)){
      t <- fcdata[fcdata$Compound == i,]
      for(j in unique(t$Strain)){
        signif <- rbind(signif, data.frame(Compound = i, Strain = j, p.value = t.test(t$foldchange[t$Strain == j], t$foldchange[t$Strain == "rho+_313"], var.equal = T)$p.value))
        
      }
    }
    signif$signif <- F
    signif$signif[signif$p.value < 0.01] <- T
    
    write.table(signif, "Data/20210816_CombinedData_Signif.tsv", sep="\t", row.names=F)
    write.table(fcsummary, "Data/20210816_CombinedData_FoldChange_Summarized.tsv", sep="\t", row.names=F)
    write.table(fcdata, "Data/20210816_CombinedData_FoldChange_Data.tsv", sep="\t", row.names=F)
    write.table(zsummary, "Data/20210816_CombinedData_ZScores_Summarized.tsv", sep="\t", row.names=F)
    write.table(zdata, "Data/20210816_CombinedData_ZScores_Data.tsv", sep="\t", row.names=F)
    }
}

## mtPotential by FACS
{
  raw <- data.frame()
  path <- "Data/"
  
  for(i in list.files(path, pattern = "export_Specimen")){
    t <- read.delim(paste(path, i, sep="/"), skip = 142, sep=",")
    
    raw <- rbind(raw, cbind(FileName = sub(".csv", "", i, fixed=T), t))
    print(i)
  }
  
  code <- read.delim("Data/Full Statistics.csv", sep=",")
  names(code)[1] <- "SampleName"
  code$SampleName <- sub(".fcs", "", code$SampleName, fixed=T)
  raw$SampleName <- gsub("export_|_Sin","" , str_extract(raw$FileName, "export_.+?_Sin"))
  
  raw <- merge(raw, code[,1:2], by="SampleName", all.x=T)
  data <- cbind(raw, Medium = NA, StrainID = NA)
  data$Medium[grepl("SM_", data$SampleID)] <- "SM"
  data$Medium[grepl("QERL_", data$SampleID)] <- "SM+QERL"
  data$StrainID <- sub(": ", "", stringr::str_extract(string = data$SampleID, pattern = ": .+"), fixed=T)
  data <- data[!is.na(data$StrainID),]
  data$StrainID <- sub(" c.+", "", data$StrainID)
  data$StrainID <- factor(data$StrainID)
  data$StrainID <- factor(data$StrainID, levels = levels(data$StrainID)[c(3,1,2,7,4,6,5)])
  data$value <- data$X530_30...BLUE.A
  data$Replicate <- sub(" c", "", str_extract(data$SampleID, " c\\d"))
  
  cells <- data %>% dplyr::group_by(SampleName) %>% dplyr::count(SampleName)
  names(cells)[2] <- "cnt"
  data <- merge(data, cells, by="SampleName", all.x=T)
  write.table(data, "Data/20210817_mtPotential_FACS.xls", sep = "\t", row.names=F)
}

## mtPotential by pMitoLoc
{
  # 130805
  {
  out <- data.frame()
  for(i in list.files(path = "Data/pMitoLoc/", pattern="JV_", full.names = T)){
    t <- read.csv(i, fileEncoding = "latin1", stringsAsFactors = T)
    t$Pearsons.Correlation <- as.numeric(as.character(t$Pearsons.Correlation))
    out <- rbind(out, t)
  }

  data <- out
  levels(data$Item.Name) <- sub(".dv (cropped)", "", levels(data$Item.Name), fixed=T)
  levels(data$Item.Name) <- sub("130805_", "", levels(data$Item.Name))
  levels(data$Item.Name) <- sub("D3D_ALX_", "", levels(data$Item.Name))
  data <- cbind(data[,c(2,25,31)], t(as.data.frame(strsplit(as.character(data$Item.Name), "_"))))
  names(data)[4:6] <- c("sample", "Image", "Cell")
  code <- read.csv("Data/pMitoLoc/130805_sample_code.csv", sep="\t")
  code$sample <- paste("sample0",code$sample, sep="")
  data <- merge(data, code, by="sample")
  
  ttest <- data.frame()
  for(i in unique(data$genotype)){
    ttest <- rbind(ttest, data.frame(Reference="Rho_313",Strain=i, pvalue=t.test(data$Pearsons.Correlation[which(data$genotype == i)], data$Pearsons.Correlation[which(data$genotype == "Rho_313")])$p.value))
    ttest <- rbind(ttest, data.frame(Reference="WT_313",Strain=i, pvalue=t.test(data$Pearsons.Correlation[which(data$genotype == i)], data$Pearsons.Correlation[which(data$genotype == "WT_313")])$p.value))
  }
  
  data$Strain <- factor(data$genotype)
  levels(data$Strain) <- c("rho0 pRS313", "rho0 pRS313_ATP3", "rho0 pRS313_ATP3-7", "rho0 pRS313_ATP3-67", "rho0 pRS313_ATP3-6", "none", "none", "wt pRS313")
  data <- data %>%
    dplyr::filter(!data$Strain == "none")
  data <- cbind(data, t(as.data.frame(strsplit(as.character(data$Strain), " "))))
  names(data)[9:10] <- c("genotype", "plasmid")
  
  
  write.table(ttest, "Data/mtPotential_YSBN11_130805_ttest.csv", sep=",", row.names=F)
  write.table(data, "Data/mtPotential_YSBN11_130805_raw.csv", sep=",", row.names=F)
  
  
  summary <- cbind(aggregate(data$Pearsons.Correlation, by=list(data$genotype), FUN=mean, na.rm=T), aggregate(data$Pearsons.Correlation, by=list(data$genotype), FUN=sd, na.rm=T)[2], aggregate(data$Pearsons.Correlation, by=list(data$genotype), FUN=stderr, na.rm=T)[2], aggregate(data$Pearsons.Correlation, by=list(data$genotype), FUN=length)[2]) 
  names(summary) <- c("Strain", "mean", "sd", "se", "n")
  summary <- merge(summary, ttest, by="Strain", all.x=T)
  summary$Strain <- factor(summary$Strain)
  summary$Reference <- factor(summary$Reference)
  levels(summary$Strain) <- c("rho0 pRS313", "rho0 pRS313_ATP3", "rho0 pRS313_ATP3-7", "rho0 pRS313_ATP3-67", "rho0 pRS313_ATP3-6", "wt pRS313")
  levels(summary$Reference) <- c("rho0 pRS313", "wt pRS313")
  summary <- cbind(summary, t(as.data.frame(strsplit(as.character(summary$Strain), " "))))
  names(summary)[8:9] <- c("genotype", "plasmid")
  write.table(summary, "Data/mtPotential_YSBN11_130805_summary.csv", sep=",", row.names=F)
  }
  
  # 130815
  {
    #read data
    data <- read.csv("Data/JV_130815_Coloc_2_s01.csv", header=T, sep=",", fileEncoding="latin1")
    data <- data[0,]
    
    for(i in 01:14){
      if(i<10) input <- paste("Data/JV_130815_Coloc_2_s0",i,".csv", sep="")
      if(i>=10) input <- paste("Data/JV_130815_Coloc_2_s",i,".csv", sep="")
      t <- read.csv(input, header=T, sep=",", fileEncoding="latin1")
      data <- rbind(data,t)
    }
    
    names(data)[2] <- "file"
    
    #assign genotypes
    data <- cbind(data, sample=as.character(data$file))
    data$sample <- as.character(data$sample)
    data$sample <- sub("130815_Conv_S","" ,data$sample)
    data$sample <- sub("_D3D_ALX.dv","" ,data$sample)
    data$sample <- sub(" (cropped)","" ,data$sample, fixed=T)
    t <- t(as.data.frame(strsplit(data$sample, "_", fixed=T)))
    data <- data[,-32]
    data <- cbind(data, sample=t[,1], cellID=paste("130815", t[,1], t[,2], t[,3], sep="_"))
    
    code <- read.csv("Data/130815_sample_code_2.csv", sep=",")
    code$sample <- as.character(code$sample)
    code[nchar(code$sample) == 1,1] <- paste("0", code[nchar(code$sample) == 1,1], sep="")
    data$sample <- as.character(data$sample)
    data[which(data$sample == "9"),32] <- "09"
    data <- merge(data, code, by="sample", all=F)
    
    data$file <- sub(".dv (cropped)", "",data$file, fixed=T)
    
    #read CCCP data
    CCCPdata <- read.csv("Data/130910_shortdata_all.csv")
    CCCPdata <- cbind(CCCPdata, genotype=paste(CCCPdata$CCCP, CCCPdata$time, sep="_"))
    CCCPdata <- CCCPdata[,3:13]
    
    #select relevant data
    shortdata <- data[,c(34,35,36,37,38,39,32,11,12,13,21,22,23)]
    shortdata <- cbind(shortdata, shortdata$Max..0./shortdata$Max..2.)
    names(shortdata)[8:14] <- c("GFP_Min", "GFP_Max", "GFP_Mean","RFP_Min", "RFP_Max", "RFP_Mean","GFP_ratio")
    
    shortdata <- cbind(shortdata, totalIntensity=shortdata$GFP_Max+shortdata$RFP_Max)
    shortdata <- cbind(shortdata, GFP.RFP.ratio=shortdata$GFP_ratio/max(shortdata$GFP_ratio))
    
    #construct short genotype
    shortdata <- cbind(shortdata, genotype=shortdata$Rho)
    shortdata$genotype <- as.character(shortdata$genotype)
    shortdata[which(!shortdata$ATP3 == "wt"),17] <- paste(shortdata[which(!shortdata$ATP3 == "wt"),17]," ATP3" ,shortdata[which(!shortdata$ATP3 == "wt"),3], sep="")
    shortdata[which(!shortdata$SYP1 == "wt"),17] <- paste(shortdata[which(!shortdata$SYP1 == "wt"),17]," SYP1" ,shortdata[which(!shortdata$SYP1 == "wt"),4], sep="")
    shortdata[which(!shortdata$Chr13 == "diploid"),17] <- paste(shortdata[which(!shortdata$Chr13 == "diploid"),17]," Chr13", sep="")
    shortdata[which(!shortdata$Chr16 == "diploid"),17] <- paste(shortdata[which(!shortdata$Chr16 == "diploid"),17]," Chr16", sep="")
    shortdata1 <- shortdata
    shortdata <- shortdata[,7:17]
    shortdata <- rbind(shortdata, CCCPdata)
    shortdata <- shortdata[which(shortdata$genotype %in% c("wt", "Rho0", "15_6")),]
    
    #restrict n number to max
    cutoff <- min(summary(factor(shortdata$genotype)))
    shortdata <- shortdata[order(shortdata$genotype, -shortdata$totalIntensity),]
    out <- shortdata[0,]
    for(i in unique(shortdata$genotype)){
      t <- shortdata[which(shortdata$genotype == i),]
      t <- t[1:cutoff,]
      out <- rbind(out, t)
    }
    shortdata <- out[,1:11]
    
    #summarize data
    summary <- cbind(aggregate(shortdata$Global.Pearsons.Correlation, by=list(shortdata$genotype), FUN=mean), aggregate(shortdata$Global.Pearsons.Correlation, by=list(shortdata$genotype), FUN=se))
    summary <- summary[,c(1,2,4)]
    names(summary) <- c("genotype", "GPC_mean", "GPC_se")
    
    #calculate ttest
    ttest <- data.frame(genotype=0, pValue=0)
    ttest <- ttest[0,]
    for(i in unique(shortdata$genotype)){
      t <- t(as.data.frame(unlist(t.test(shortdata[which(shortdata$genotype == i),1], shortdata[which(shortdata$genotype == "wt"),1]))))
      t <- cbind(genotype=i, pValue=t[,3])
      ttest <- rbind(ttest, t)
    }
    ttest$pValue <- as.numeric(as.character(ttest$pValue))
    #ttest$pValue <- round(ttest$pValue, digits=5)
    ttest <- merge(ttest, summary, by=c("genotype"))
    
    #calculate n numbers
    n <- data.frame(genotype=0, n=0)
    n <- n[0,]
    for(i in unique(shortdata$genotype)){
      t <- cbind(genotype=i, n=length(shortdata[which(shortdata$genotype == i),1]))
      n <- rbind(n,t)
    }
    
    #calculate statistics for single clones
    summary_single <- cbind(aggregate(shortdata1$Global.Pearsons.Correlation, by=list(shortdata1$clone), FUN=mean), aggregate(shortdata1$Global.Pearsons.Correlation, by=list(shortdata1$clone), FUN=se))
    summary_single <- summary_single[,c(1,2,4)]
    names(summary_single) <- c("clone", "GPC_mean", "GPC_se")
    idata <- shortdata1 %>%
      dplyr::select(clone, Global.Pearsons.Correlation) %>%
      dplyr::rename(GPC = Global.Pearsons.Correlation)
      
    t <- shortdata1[!duplicated(shortdata1$clone),c(1,17)]
    summary_single <- merge(summary_single, t, by="clone")
    
    #reorder factor levels
    summary_single$clone <- factor(summary_single$clone, levels=c("dipl_ctrl","init","c02","c05","c10","c18","c12","c13","c15","c16","c17","c20","c06","wt"))
    idata$clone <- factor(idata$clone, levels=c("dipl_ctrl","init","c02","c05","c10","c18","c12","c13","c15","c16","c17","c20","c06","wt"))
    
    
    ttest_single <- data.frame(clone=0, pValue=0)
    ttest_single <- ttest_single[0,]
    for(i in unique(shortdata1$clone)){
      t <- t(as.data.frame(unlist(t.test(shortdata1[which(shortdata1$clone == i),7], shortdata1[which(shortdata1$clone == "init"),7]))))
      t <- cbind(clone=i, pValue=t[,3])
      ttest_single <- rbind(ttest_single, t)
    }
    ttest_single$pValue <- as.numeric(as.character(ttest_single$pValue))
    #ttest$pValue <- round(ttest$pValue, digits=5)
    ttest_single <- merge(ttest_single, summary_single, by=c("clone"))
    
    n_s <- data.frame(clone=0, n=0)
    n_s <- n_s[0,]
    for(i in unique(shortdata1$clone)){
      t <- cbind(clone=i, n=length(shortdata1[which(shortdata1$clone == i),1]))
      n_s <- rbind(n_s,t)
    }
    
    write.table(summary_single, "Data/mtPotential_YSBN1_130815_summary.csv", sep=",", quote=F, row.names=F)
    write.table(ttest_single, "Data/mtPotential_YSBN1_130815_ttest.csv", sep=",", quote=F, row.names=F)
    write.table(n_s, "Data/mtPotential_YSBN1_130815_n.csv", sep=",", quote=F, row.names=F)
    write.table(summary_single[c(1,4)], "Data/mtPotential_YSBN1_130815_code.csv", sep=",", row.names=F)
    
    #mtPotential evo isolates
    data <- read.csv("Data/mtPotential_YSBN1_130815_summary.csv")
    n <- read.csv("Data/mtPotential_YSBN1_130815_n.csv")
    ttest <- read.csv("Data/mtPotential_YSBN1_130815_ttest.csv")
    
    data <- merge(data, ttest[,1:2], by="clone")
    data <- merge(data, n, by="clone")
    data$clone <- factor(data$clone)
    x <- levels(data$clone)[c(13,1,2,5,11,14)]
    data <- data[which(data$clone %in% x),]
    write.table(data, "Data/mtPotential_YSBN1_130815_selected_summary.csv", sep=",", row.names=F)
    
    sdata <- idata %>%
      dplyr::filter(clone %in% levels(data$clone)[c(13,1,2,5,11,14)]) %>%
      write.table(., "Data/mtPotential_YSBN1_130815_selected_data.csv", sep=",", row.names = F)
    
  }
}

## mitochondrial fragmentation by pMitoLoc
{
  data <- read.csv("Data/151108_mitoLoc_rho0_313_fulldata.csv", stringsAsFactors = T)
  
  FormFactor <- data[,c(-3,-4,-5)]
  
  #add relative surface area
  t <- aggregate(FormFactor$Surface.Area..Âµm.., by=list(FormFactor$Label), FUN=sum)
  names(t) <- c("Label", "totalSurface")
  FormFactor <- merge(FormFactor, t, by="Label", all.x=T)
  FormFactor <- cbind(FormFactor, relativeSurface = FormFactor$Surface.Area..Âµm../FormFactor$totalSurface)
  
  #weight parameters of each object by relative Surface
  
  FormFactor <- cbind(FormFactor, 
                      SurfaceAreaWeighted=FormFactor$Surface.Area..Âµm..*(FormFactor$relativeSurface), 
                      CompactnessWeighted=FormFactor$Compactness..CP.*(FormFactor$relativeSurface), 
                      DistributionIsotropyWeighted=FormFactor$Distribution.Isotropy*(FormFactor$relativeSurface),
                      IsoperimetricQuotientWeighted=FormFactor$Isoperimetric.Quotient*(FormFactor$relativeSurface),
                      SphericityWeighted=FormFactor$Sphericity*(FormFactor$relativeSurface),
                      SAVWeighted=FormFactor$SA.V*(FormFactor$relativeSurface),
                      RadiusVarianceWeighted=FormFactor$Radius.Variance*(FormFactor$relativeSurface)            
  )
  
  FormFactor <- melt(FormFactor[,c(-13,-17,-18)], id=c("Label", "Date", "Strain", "Plasmid", "Image", "Object", "ConditionID"))
  
  
  FormFactorSummary <- cbind(aggregate(FormFactor$value, by=list(FormFactor$variable, FormFactor$Label), FUN=sum), aggregate(FormFactor$value, by=list(FormFactor$variable, FormFactor$Label), FUN=mean)[3])
  names(FormFactorSummary) <- c("variable", "Label", "sum", "mean")
  FormFactorSummary <- merge(FormFactorSummary, FormFactor[!duplicated(FormFactor$Label),c(1,7)], by="Label", all.x=T )
  
  FormFactorSummaryCast <- cast(FormFactorSummary, Label ~ variable, value="sum")
  FormFactorSummaryCast <- merge(FormFactorSummaryCast, FormFactorSummary[!duplicated(FormFactorSummary$Label),c(-2,-3,-4)], by="Label", all.x=T)
  
  
  #add fragmentation index
  Frag <- data[,c(1,2,15:18,23,25)]
  Frag2 <- Frag[which(Frag$normVolume < 20),]
  Frag2 <- aggregate(Frag2$normVolume, by=list(Frag2$Label), FUN=sum)
  names(Frag2) <- c("Label", "Fragmentation")
  Frag <- cbind(aggregate(Frag$normVolume, by=list(Frag$Label), FUN=mean), aggregate(Frag$normVolume, by=list(Frag$Label), FUN=sd)[2])
  names(Frag) <- c("Label", "Volume_mean", "Volume_sd")
  Frag <- merge(Frag, Frag2, by="Label", all.x=T)
  Frag$Volume_sd[is.na(Frag$Volume_sd)] <- 0
  Frag$Fragmentation[is.na(Frag$Fragmentation)] <- 0
  
  
  FormFactorSummaryCast <- merge(FormFactorSummaryCast, Frag, by="Label")
  
  matrix <- FormFactorSummaryCast[,c(16,17,18,20,21,24,26)]
  
  matrix.class <- FormFactorSummaryCast[,23]
  
  matrix.pca <- prcomp(matrix, scale. = T)
  p <- ggbiplot(matrix.pca, obs.scale = 1, var.scale = 1, groups = matrix.class, ellipse = TRUE, circle = TRUE)+theme_classic()
  print(p)
  #ggsave("151108_PCA.pdf")

  #calculate single f scores 
  data2 <- read.csv("Data/151108_mitoLoc_rho0_313_processed_data.csv")
  summary <- data2[which(data2$cutoff < 20),]
  summary$Item.Name <- as.character(summary$Item.Name)
  #expand by non-existing cells
  t <- data2$Item.Name[which(!data2$Item.Name %in% summary$Item.Name)]
  for(i in t){
    summary <- rbind(summary, data.frame(Item.Name = i, ConditionID=as.character(unique(data2$ConditionID[which(data2$Item.Name == i)])),  binnedPercentages = 0, cutoff = 5))
  }
  summary <- aggregate(summary$binnedPercentages, by=list(summary$ConditionID, summary$Item.Name), FUN=sum)
  names(summary) <- c("ConditionID", "Item.Name", "f")
  
  summary_fscore <- cbind(aggregate(summary$f, by=list(summary$ConditionID), FUN=mean),aggregate(summary$f, by=list(summary$ConditionID), FUN=se)[2])
  names(summary_fscore) <- c("ConditionID", "f_mean", "f_sem")
  write.table(summary_fscore, "Data/151108_fscore_mean.csv", sep=",", row.names=F)
  write.table(summary, "Data/151108_fscore_data.csv", sep=",", row.names=F)
  
  out <- data.frame()
  for(i in unique(summary$ConditionID)){
    wt <- t.test(summary$f[which(summary$ConditionID == "wt_313")], summary$f[which(summary$ConditionID == i)])$p.value
    rho <- t.test(summary$f[which(summary$ConditionID == "rho0_313")], summary$f[which(summary$ConditionID == i)])$p.value
    out <- rbind(out, data.frame(ConditionID=i, ttest_wt=wt, ttest_rho0=rho))
  }
  write.table(out, "Data/151108_fscore_ttest.csv", sep=",", row.names=F)
  
  summary_fscore <- merge(summary_fscore, out, by="ConditionID")
  
  ggplot(summary_fscore, aes(x = ConditionID, y = f_mean))+
    geom_bar(stat="identity")+
    geom_errorbar(aes(ymin=f_mean-f_sem, ymax=f_mean+f_sem), width=0.2)
  
  write.table(summary_fscore, "Data/151108_fscore.csv", sep=",", row.names=F)
}

## Growth Curves (Strains)
{
  # evolved strains
  {
    source("Scripts/MGrofit2.R")
    MGrofit(file="20130526_JH_growth curves  24 clones f1_raw data.csv", fit="s", Interactive = F, makeMean = F, makePlot = F)
    data <- read.csv("Data/Grofit_spline_20130526_JH_growth curves  24 clones f1_raw data.csv", stringsAsFactors = T)
    code <- read.xlsx("Data/20130526_JH_growth curves  24 clones f1_code.xlsx", sheetIndex=1)
    code <- melt(code, id="NA.")
    code <- cbind(code, Position=paste(code$NA., sub("X", "", as.character(code$variable)), sep=""))
    data <- cbind(data, Position=paste(data$Row, data$Column, sep=""))
    data <- merge(data[c(8:10,12)], code[,3:4], by="Position", all.x=T)
    names(data)[5] <- "clone"
    data$clone <- factor(data$clone)
    levels(data$clone)[!is.na(as.numeric(levels(data$clone))) & nchar(levels(data$clone)) == 1] <- paste("c0", levels(data$clone)[!is.na(as.numeric(levels(data$clone))) & nchar(levels(data$clone)) == 1], sep="")
    levels(data$clone)[!is.na(as.numeric(levels(data$clone)))] <- paste("c", levels(data$clone)[!is.na(as.numeric(levels(data$clone)))], sep="")
    data <- data[which(data$clone %in% c("init", "c02", "c05", "c12", "c20", "wt")),]
    
    code <- read.csv("Data/mtPotential_YSBN1_130815_code.csv")
    data <- merge(data, code, by="clone", all.x=T)
    data$mu.spline[which(data$clone == "c20")][2] <- NA
    data$mu.spline[which(data$clone == "wt")][2] <- NA
    
    write.table(data, "Data/Growcurves_YSBN11_130526_fit_values.csv", sep="\t", row.names=F)
    
    summary <- cbind(aggregate(data$mu.spline, by=list(data$clone, data$genotype), FUN=mean, na.rm=T), aggregate(data$mu.spline, by=list(data$clone, data$genotype), FUN=sd, na.rm=T)[3])
    names(summary) <- c("clone","ID","mean", "sd")
    

    ttest <- data.frame()
    for(i in unique(data$clone)){
      ttest <- rbind(ttest, data.frame(Reference="init",clone=i, pvalue=t.test(data$mu.spline[which(data$clone == i)], data$mu.spline[which(data$clone == "init")])$p.value))
    }
    write.table(ttest, "Data/Growcurves_YSBN11_130526_fit_values_ttest.csv", sep=",", row.names=F)
    
    summary <- cbind(summary, genotype="rho0 evolved")
    summary$genotype <- as.character(summary$genotype)
    summary$genotype[grepl("wt", summary$clone)] <- "wild type"
    summary$genotype[grepl("init", summary$clone)] <- "rho0 initial"
    
    summary <- cbind(summary, ATP3=c("ATP3-9", "ATP3-7", "ATP3-8", "ATP3-6", "wild type", "wild type"))
    write.table(summary, "Data/Growcurves_YSBN11_130526_fit_values_summary.csv", sep=",", row.names=F)
    
    ggplot(summary, aes(x = clone, y = mean))+
      geom_bar(stat="identity")+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2)
  }
  
  # model strains
  {
    source("Scripts/MGrofit2.R")
    MGrofit(file="20130914_YSBN11_ATP&SYP_F1_0.5%Glc_data_raw.csv", fit="m", Interactive = F, makeMean = F, makePlot = F)
    data <- read.csv("Data/Grofit_model_20130914_YSBN11_ATP&SYP_F1_0.5%Glc_data_raw.csv", stringsAsFactors = T)
    code <- read.xlsx("Data/20130914_YSBN11_ATP&SYP_F1_0.5%Glc_code.xlsx", sheetIndex=1)
    code$Sample <- factor(code$Sample)
    levels(data$Sample) <- sub("Sample ", "", levels(data$Sample))
    levels(code$Sample) <- sub("x", "X", levels(code$Sample))
    code$Strain <- paste(code$Strain, code$insert, code$clone)
    data <- merge(data, code, by="Sample", all.x=T)
    data <- data[which(data$Strain %in% c("rho0_pRS313 NA 1", "rho0_pRS313 ATP3 G919C 1", "rho0_pRS313 ATP3 G919C 2", "rho0_pRS313 ATP3 mut 1", "rho0_pRS313 ATP3 T911A 1", "rho0_pRS313 ATP3 T911A 2", "rho0_pRS313 ATP3 2", "wt_pRS313 NA 1")),]
    data$mu.model[which(data$Strain %in% c("rho0_pRS313 ATP3 T911A 1"))][2:3] <- NA
    data$mu.model[which(data$Strain %in% c("rho0_pRS313 ATP3 T911A 2"))][3] <- NA
    data$mu.model[which(data$Strain %in% c("rho0_pRS313 ATP3 G919C 1"))][2] <- NA
    data$mu.model[which(data$Strain %in% c("rho0_pRS313 ATP3 G919C 2"))][1] <- NA
    
    summary <- cbind(aggregate(data$mu.model, by=list(data$Strain, data$clone), FUN=mean, na.rm=T), aggregate(data$mu.model, by=list(data$Strain, data$clone), FUN=sd, na.rm=T)[3])
    names(summary) <- c("Strain", "clone","mean", "sd")
    
    data$Strain <- factor(data$Strain)
    levels(data$Strain) <- substr(levels(data$Strain), 1, nchar(levels(data$Strain))-2)
    levels(data$Strain) <- sub(" NA", "", levels(data$Strain))
    data$Strain <- factor(data$Strain, levels=levels(data$Strain)[c(6,5,1,4,2,3)])
    write.table(data, "Data/Growcurves_YSBN11_130914_fit_values.csv", sep="\t", row.names=F)
    
    ttest <- data.frame()
    for(i in unique(data$Strain)){
      ttest <- rbind(ttest, data.frame(Reference="wt_pRS313",Strain=i, pvalue=t.test(data$mu.model[which(data$Strain == i)], data$mu.model[which(data$Strain == "wt_pRS313")])$p.value))
      ttest <- rbind(ttest, data.frame(Reference="rho0_pRS313",Strain=i, pvalue=t.test(data$mu.model[which(data$Strain == i)], data$mu.model[which(data$Strain == "rho0_pRS313")])$p.value))
    }
    write.table(ttest, "Data/Growcurves_YSBN11_130914_fit_values_ttest.csv", sep=",", row.names=F)
    
    summary <- cbind(aggregate(data$mu.model, by=list(data$Strain), FUN=mean, na.rm=T), aggregate(data$mu.model, by=list(data$Strain), FUN=sd, na.rm=T)[2])
    names(summary) <- c("Strain", "mean", "sd")
    summary <- cbind(summary, genotype="wild type")
    summary$genotype <- as.character(summary$genotype)
    summary$genotype[grepl("rho", summary$Strain)] <- "rho0"
    summary <- cbind(summary, plasmid=c("pRS313", "pRS313", "pRS313_ATP3", "pRS313_ATP3-6", "pRS313_ATP3-7", "pRS313_ATP3-67"))
    write.table(summary, "Data/Growcurves_YSBN11_130914_fit_values_summary.csv", sep=",", row.names=F)
  }
}

## Growth Curves (AA Supplementation)
{
  grofit <- F
  if(grofit == T){
    source("Scripts//MGrofit2.R")
    MGrofit(file = "Data/150917_AA_suppl.csv", Interactive = F, makeMean = F, exclude=0, fit="s", makePlot = F, fileformat = "standard")
  }
  
  rawdata <- read.csv("Data/150917_AA_suppl.csv", stringsAsFactors = T)
  rawdata  <- cbind(Well=paste(rawdata$Well.Row, rawdata$Well.Col, sep=""),rawdata[,c(-1)])
  names(rawdata) <- c("Well", "Col", "Row", paste("X", rawdata[1,4:ncol(rawdata)], sep=""))
  rawdata <- melt(rawdata[-1,], id=c("Well", "Col", "Row"))
  levels(rawdata$variable) <- sub("X", "", levels(rawdata$variable))
  rawdata$variable <- as.numeric(as.character(rawdata$variable))
  names(rawdata)[4] <- "time"

  data <- read.csv("Data/Grofit_model_150917_AA_suppl.csv", stringsAsFactors = T)
  data2 <- read.csv("Data/Grofit_spline_150917_AA_suppl.csv", stringsAsFactors = T)
  
  layout <- read.xlsx("Data/150917_AA_growthcurve_rescue.xlsx", sheetName="Layout_R")
  layout <- melt(layout, id="NA.")
  layout <- data.frame(Well = paste(as.character(layout[,1]), sub("X", "", as.character(layout[,2])), sep=""), SampleID=layout[,3])
  layout$SampleID <- factor(layout$SampleID)
  
  data_all <- data.frame(Well=paste(data[,3], data[,2], sep=""), data[,9:11], data2[,8:10])
  data_all <- merge(data_all, layout, by="Well", all=T)
  levels(data_all$SampleID) <- sub("- ", "-", levels(data_all$SampleID))
  levels(data_all$SampleID) <- sub("+ ", "+", levels(data_all$SampleID), fixed=T)
  levels(data_all$SampleID) <- sub("  ", " ", levels(data_all$SampleID), fixed=T)
  
  t <- strsplit(as.character(data_all$SampleID), " ")
  data_all <- data_all[which(!data_all$SampleID == "0 0 0"),]
  data_all <- data_all[which(!data_all$SampleID == "0 0 0 "),]
  data_all <- data_all[!grepl("BLANK", data_all$SampleID),]
  
  
  df <- data.frame(t(as.data.frame(strsplit(as.character(data_all$SampleID), " "))))
  names(df) <- c("Strain", "Mix", "Replicate")
  
  data_all <- cbind(data_all, df)
  data_all <- data_all[which(!data_all$Mix == "0"),]
  data_all <- data_all[which(!data_all$Replicate == "0"),]
  data_all <- cbind(data_all, Genotype="rho0", Plasmid="313")
  data_all$Plasmid <- as.character(data_all$Plasmid)
  data_all$Plasmid[grepl("ATP", data_all$Strain)] <- "ATP3 T911A"
  data_all$Genotype <- as.character(data_all$Genotype)
  data_all$Genotype[grepl("WT", data_all$Strain)] <- "wt"
  data_all$Strain <- "YSBN11"
  data_all <- cbind(data_all, StrainID=paste(data_all$Strain, data_all$Genotype, data_all$Plasmid))
  
  #remove outliers
  data_all <- cbind(data_all, ID=paste(data_all$Strain, data_all$Genotype, data_all$Plasmid, data_all$Mix))
  for(i in unique(data_all$ID)){
    if(grepl("BLANK", i) == F){
      t <- rm.outlier(data_all$mu.model[which(data_all$ID == i)])
      if(length(t) < length(data_all$mu.model[which(data_all$ID == i)])){
        data_all$mu.model[which(data_all$ID == i)][which(!data_all$mu.model[which(data_all$ID == i)] %in% t)] <- NA
      }
      
      t <- rm.outlier(data_all$mu.spline[which(data_all$ID == i)])
      if(length(t) < length(data_all$mu.spline[which(data_all$ID == i)])){
        data_all$mu.spline[which(data_all$ID == i)][which(!data_all$mu.spline[which(data_all$ID == i)] %in% t)] <- NA
      }
      
      t <- rm.outlier(data_all$A.model[which(data_all$ID == i)])
      if(length(t) < length(data_all$A.model[which(data_all$ID == i)])){
        data_all$A.model[which(data_all$ID == i)][which(!data_all$A.model[which(data_all$ID == i)] %in% t)] <- NA
      }
      
      t <- rm.outlier(data_all$A.spline[which(data_all$ID == i)])
      if(length(t) < length(data_all$A.spline[which(data_all$ID == i)])){
        data_all$A.spline[which(data_all$ID == i)][which(!data_all$A.spline[which(data_all$ID == i)] %in% t)] <- NA
      }
    }}
  
  data_all$mu.spline[which(data_all$Mix == "+E" & data_all$StrainID == "YSBN11 wt 313")][7] <- NA
  data_all$mu.spline[which(data_all$Mix == "+Q" & data_all$StrainID == "YSBN11 wt 313")][8] <- NA
  
  
  
  data_all_summary <- cbind(aggregate(data_all[,c(2:7)], by=list(data_all$Strain, data_all$Genotype, data_all$Plasmid, data_all$Mix), na.rm=T, FUN=mean), aggregate(data_all[,c(2:7)], by=list(data_all$Strain, data_all$Genotype, data_all$Plasmid, data_all$Mix), na.rm=T, FUN=sd)[5:10])
  t <- names(data_all_summary)[5:10]
  names(data_all_summary) <- c("Strain", "Genotype", "Plasmid", "Mix", paste(t, "_mean", sep=""), paste(t, "_sd", sep=""))
  data_all_summary <- cbind(data_all_summary, StrainID=paste(data_all_summary$Strain, data_all_summary$Genotype, data_all_summary$Plasmid, sep=" "))
  data_all_summary <- data_all_summary[which(!data_all_summary$Strain == "BLANK"),]
  data_all_summary <- cbind(data_all_summary, StrainID2 = paste(data_all_summary$Strain, data_all_summary$Genotype))
  
  t <- cbind(data_all_summary[which(data_all_summary$Mix == "-AA"),c(17,5)], max=data_all_summary[which(data_all_summary$Mix == "+AA"),c(5)])
  names(t)[2] <- "min"
  data_all_summary <- merge(data_all_summary, t, by="StrainID", all.x=T)
  data_all_summary <- cbind(data_all_summary, mu.model.mean.norm=100/(data_all_summary$max-data_all_summary$min)*(data_all_summary$mu.model_mean-data_all_summary$min), mu.model.sd.norm=100/(data_all_summary$max-data_all_summary$min)*(data_all_summary$mu.model_sd))
  
  #calculate p values
  reference <- unique(data_all$Mix)[9]
  data_all_summary <- cbind(data_all_summary, pvalue=NA)
  for(i in unique(data_all_summary$StrainID)){
    for(j in unique(data_all_summary$Mix)){
      t <- t.test(data_all$mu.model[which(data_all$StrainID == i & data_all$Mix == j)],data_all$mu.model[which(data_all$StrainID == i & data_all$Mix == reference)])$p.value
      data_all_summary$pvalue[which(data_all_summary$StrainID == i & data_all_summary$Mix == j)] <- t
    }
  }
  
  
  #calculate growth defect as % of wt growth rate
  t <- data_all_summary[which(data_all_summary$StrainID == "YSBN11 wt 313"),c(5,9)]
  names(t)[2] <- "norm"
  data_all_summary <- merge(data_all_summary, t, by="Mix")
  data_all_summary <- cbind(data_all_summary[,-24], defect_spline_mean=100/data_all_summary$norm*data_all_summary$mu.spline_mean, defect_spline_sd=100/data_all_summary$norm*data_all_summary$mu.spline_sd)
  data_all_summary <- cbind(data_all_summary, defect_spline_mean_sub=100-data_all_summary$defect_spline_mean)
  
  write.table(data_all_summary, "Data/150917_growth_outliersrm.csv", sep=",", row.names=F, quote=F)
  write.table(data_all, "Data/150917_growth_outliersrm_raw.csv", sep=",", row.names=F, quote=F)

  # calculate growth defect in raw data
  {
    data_all <- merge(data_all, t, by ="Mix", all.x=T)
    
    data_all <- data_all %>%
      dplyr::mutate(defect_spline = 100/norm*mu.spline) %>%
      dplyr::mutate(defect_spline_sub = 100-defect_spline)
    write.table(data_all, "Data/150917_growth_outliersrm_raw_defect.csv", sep=",", row.names=F, quote=F)
    
    }
  
  
}

## Growth Curves (aconitase KO)
{
  data <- read.csv("Data/150521_growth_aco_replicates.csv", stringsAsFactors = T)
  data <- data[which(data$Concentration == "0mM"),c(2,9:12)]
  data <- cbind(data, ID=paste(data$Strain, data$Genotype, data$Plasmid))
  summary <- cbind(aggregate(data$mu.model, by=list(data$Strain, data$Genotype, data$Plasmid, data$ID), FUN=mean), aggregate(data$mu.model, by=list(data$Strain, data$Genotype, data$Plasmid, data$ID), FUN=sd)[5])
  names(summary) <- c("Strain", "Genotype", "Plasmid", "ID", "mean", "sd")
  summary <- cbind(summary, pvalue=NA, signif=NA)
  for(i in unique(data$ID)){
    summary$pvalue[which(summary$ID == i)] <- t.test(data$mu.model[which(data$ID == i)], data$mu.model[which(data$ID == "BY4741 wt pRS313")])$p.value
  }
  summary$signif[which(summary$pvalue < 0.01)] <- "*"
  summary$Genotype <- as.character(summary$Genotype)
  summary$Genotype[1:2] <- c("aco1", "aco1")
  
  write.table(summary, "Data/growth_rates_aco.csv", sep=",")
  write.table(data, "Data/growth_rates_aco_raw.csv", sep=",")
}

## ICP-MS Analysis
{
  data <- read.xlsx("Data/2015-10-20 Fe Jakob.xls", sheetIndex=1, header=F)
  
  sample <- data[is.na(data[,3]),2]
  
  for(i in 1:length(sample)){
    t <- data[((i*22)-21):(i*22),]
    t <- t[-1,]
    t <- t[which(!t[,3] %in% c("63", "68", "60")),]
    if(i == 1){
      output <- t[c(-1,-2),c(2,3,4)]
      names(output) <- c("Analyte", "Mass",as.character(sample[i]))  
    } else {
      output <- cbind(output,t[c(-1,-2),c(4)])
      names(output)[i+2] <- as.character(sample[i])    
    }
  }
  
  for(i in 1:8){
    names(output)[grepl(paste("JV-5",i,sep=""), names(output))] <- paste(names(output)[grepl(paste("JV-5",i,sep=""), names(output))], seq(1:length(names(output)[grepl(paste("JV-5",i,sep=""), names(output))])))
  }
  
  
  output <- melt(output, id=c("Analyte", "Mass"))
  output$value <- as.numeric(as.character(output$value))
  
  output <- cbind(output, SampleType="Blank")
  output$SampleType <- as.character(output$SampleType)
  output$SampleType[grepl("JV",output$variable)] <- "Sample"
  output$SampleType[grepl("Standard",output$variable)] <- "Standard"
  output$SampleType[grepl("BIR",output$variable)] <- "Standard"
  output$SampleType[grepl("neat",output$variable)] <- "Sample undiluted"
  
  SampleData <- output[which(output$SampleType == "Sample"),]
  SampleData <- cbind(SampleData, Strain=NA)
  SampleData$Strain[grepl("JV-51", SampleData$variable)] <- "WT 313"
  SampleData$Strain[grepl("JV-52", SampleData$variable)] <- "Rho 313"
  SampleData$Strain[grepl("JV-53", SampleData$variable)] <- "Rho T911A"
  SampleData$Strain[grepl("JV-54", SampleData$variable)] <- "LB"
  SampleData$Strain[grepl("JV-55", SampleData$variable)] <- "WT 313"
  SampleData$Strain[grepl("JV-56", SampleData$variable)] <- "Rho 313"
  SampleData$Strain[grepl("JV-57", SampleData$variable)] <- "Rho T911A"
  SampleData$Strain[grepl("JV-58", SampleData$variable)] <- "LB"
  SampleData <- SampleData[which(!SampleData$Analyte %in% c("In", "Re", "Rh")),]
  SampleData <- cbind(SampleData, replicate=NA)
  SampleData$replicate[grepl(" 1", SampleData$variable)] <- 1
  SampleData$replicate[grepl(" 2", SampleData$variable)] <- 2
  SampleData$replicate[grepl(" 3", SampleData$variable)] <- 3
  SampleData$replicate[grepl(" 4", SampleData$variable)] <- 4
  SampleData$replicate[grepl(" 5", SampleData$variable)] <- 5
  SampleData$replicate[grepl(" 6", SampleData$variable)] <- 6
  SampleData <- cbind(SampleData, extraction=NA)
  SampleData$extraction[grepl("JV-51", SampleData$variable)] <- 1
  SampleData$extraction[grepl("JV-52", SampleData$variable)] <- 1
  SampleData$extraction[grepl("JV-53", SampleData$variable)] <- 1
  SampleData$extraction[grepl("JV-54", SampleData$variable)] <- 1
  SampleData$extraction[grepl("JV-55", SampleData$variable)] <- 2
  SampleData$extraction[grepl("JV-56", SampleData$variable)] <- 2
  SampleData$extraction[grepl("JV-57", SampleData$variable)] <- 2
  SampleData$extraction[grepl("JV-58", SampleData$variable)] <- 2
  
  SumData <- aggregate(SampleData$value, by=list(SampleData$Analyte, SampleData$Mass, SampleData$Strain, SampleData$replicate), FUN=sum)
  names(SumData) <- c("Analyte", "Mass","Strain", "replicate", "sum")
  
  #subtract blank measurements
  t <- aggregate(SumData$sum, by=list(SumData$Analyte, SumData$Strain), FUN=mean)
  names(t) <- c("Analyte", "Strain", "blank")
  t <- t[which(t$Strain == "LB"),]
  SumData <- merge(SumData[which(!SumData$Strain == "LB"),], t[,-2], by="Analyte", all=T)
  SumData <- cbind(SumData, sum_blanked=SumData$sum-SumData$blank)
  
  BCA <- read.xlsx("Data/141029_ICP.xlsx", sheetName="BCA")
  SumData <- merge(SumData, BCA, by="Strain", all=T)
  SumData <- cbind(SumData, ng.Âµg_protein = SumData$sum/SumData$BCA)
  SumData$Mass <- as.numeric(as.character(SumData$Mass))
  SumData <- cbind(SumData, nmol.Âµg_protein = SumData$ng.Âµg_protein/SumData$Mass)
  
  ttest <- data.frame()
  for(i in unique(SumData$Analyte)){
    for(j in unique(SumData$Strain)){
      t <- t.test(SumData$nmol.Âµg_protein[which(SumData$Strain == j & SumData$Analyte == i)], SumData$nmol.Âµg_protein[which(SumData$Strain == "WT 313" & SumData$Analyte == i)])$p.value
      ttest <- rbind(ttest, data.frame(Strain=j, Analyte=i, pvalue=t))
    }}
  
  SumData[which(SumData$Strain == "WT 313" & SumData$Analyte=="Fe"),]
  SumData <- SumData[which(SumData$replicate %in% c(1:4)),]
  
  MeanData <- cbind(aggregate(SumData$nmol.Âµg_protein, by=list(SumData$Strain, SumData$Analyte), FUN=mean), aggregate(SumData$nmol.Âµg_protein, by=list(SumData$Strain, SumData$Analyte), FUN=sd)[3])
  names(MeanData) <- c("Strain", "Analyte", "mean.nmol.Âµg_protein", "sd.nmol.Âµg_protein")
  MeanData <- merge(MeanData, ttest, by=c("Strain", "Analyte"), all.x=T)
  MeanData <- cbind(MeanData, signif=NA)
  MeanData$signif[which(MeanData$pvalue<0.01)] <- "*"
  MeanData$Strain <- factor(MeanData$Strain)
  MeanData$Strain <- factor(MeanData$Strain, levels=levels(MeanData$Strain)[c(3,1,2)])
  
  write.table(MeanData, "Data/151019_ICP_Analysis_MeanData.csv", sep=",")
  write.table(SumData, "Data/151019_ICP_Analysis_SumData.csv", sep=",")
  
}

## 55Fe Incorporation
{
  data <- read.xlsx("Data/160211_ISC_labelling.xlsx", sheetName="results", startRow = 2, colIndex = 1:12)
  
  data <- melt(data, id=c("sample", "replicate", "dilution", "blank"))
  data <- cbind(data, value_blanked=data$value-data$blank)
  data <- cbind(data, value_dil = data$value_blanked*data$dilution)
  
  summary <- cbind(aggregate(data$value_dil, by=list(data$sample, data$variable), FUN=mean, na.rm=T), aggregate(data$value_dil, by=list(data$sample, data$variable), FUN=sd, na.rm=T)[3])
  names(summary) <- c("sample", "variable", "mean", "sd")
  #attach p values
  output <- data.frame()
  for(j in unique(data$sample)[c(1,4,5,8)]){
    for(i in unique(data$variable)){
      output <- rbind(output,
                      data.frame(sample=j, 
                                 variable = i, 
                                 pvalue_aco=t.test(data$value_dil[which(data$variable == i & data$sample == j)], data$value_dil[which(data$variable == "BY4741.pUG35_Aco1.pRS313" & data$sample == j)])$p.value,
                                 pvalue_leu=t.test(data$value_dil[which(data$variable == i & data$sample == j)], data$value_dil[which(data$variable == "BY4741.p416GPD_StrepII.Leu1.pRS313.pLM" & data$sample == j)])$p.value))
    }}
  summary <- merge(summary, output, by=c("sample", "variable"), all.x=T)
  
  names(summary)[1:2] <- c("fraction", "sample")
  names(data)[c(1,5)] <- c("fraction", "sample")
  
  levels(summary$sample) <- c("pUG35_Aco1 wt", "pUG35_Aco1 rho0", "pUG35_Aco1 rho0 ATP3-6", "p416GPD_StrepII-Leu1 wt", "p416GPD_StrepII-Leu1 rho0", "p416GPD_StrepII-Leu1 rho0 ATP3-6", "pUG35_Leu1", "pUG35")
  summary$fraction <- factor(summary$fraction)
  summary$fraction <- factor(summary$fraction, levels=levels(summary$fraction)[c(3,4,5,2,6,7,8,1)])
  
  levels(data$sample) <- c("pUG35_Aco1 wt", "pUG35_Aco1 rho0", "pUG35_Aco1 rho0 ATP3-6", "p416GPD_StrepII-Leu1 wt", "p416GPD_StrepII-Leu1 rho0", "p416GPD_StrepII-Leu1 rho0 ATP3-6", "pUG35_Leu1", "pUG35")
  data$fraction <- factor(data$fraction)
  data$fraction <- factor(data$fraction, levels=levels(data$fraction)[c(3,4,5,2,6,7,8,1)])
  
  
  #correct for protein concentration: 10mg lysate/sample
  summary <- cbind(summary, mean_cor=summary$mean/7.5, sd_cor=summary$sd/7.5)
  data <- cbind(data, cor=data$value_dil/7.5)
  
  #correct for lysate iron
  t <- summary[which(summary$fraction == "lysate II"),c(2,7)]
  summary <- merge(summary, t, by="sample", all.x=T)
  names(summary)[c(7,9)] <- c("mean_cor", "ref")
  summary <- cbind(summary, mean_norm=summary$mean/summary$ref*10, sd_norm=summary$sd_cor/summary$ref*10)
  
  data <- merge(data, t, by = "sample", all.x=T)
  names(data)[10] <- "ref"
  data <- cbind(data, norm=data$cor/data$ref*10)
  
  
  summary <- cbind(summary, category=NA)
  summary$category[which(summary$fraction %in% levels(summary$fraction)[1:3])] <- "medium"
  summary$category[which(!summary$fraction %in% levels(summary$fraction)[1:3])] <- "lysate"
  data <- cbind(data, category=NA)
  data$category[which(data$fraction %in% levels(data$fraction)[1:3])] <- "medium"
  data$category[which(!data$fraction %in% levels(data$fraction)[1:3])] <- "lysate"
  
  
  write.table(summary, "Data/160211_ISC_labelling_summary.tsv", sep="\t")
  write.table(data, "Data/160211_ISC_labelling_data.tsv", sep="\t")
  
}

## Enzyme activities
{
  ### Aconitase: mtAco1
  {
  data <- read.csv("Data/150729_BY_mtAco_output.csv", stringsAsFactors = T)
  data <- cbind(data, ID=paste(data$Strain, data$Genotype, data$Replicate))
  data <- merge(data[which(data$Substrate == "+citrate"),c(4,9)], data[which(data$Substrate == "-citrate"),c(4,9)], by=c("ID"))
  data <- cbind(data, blanked=data$slope.x-data$slope.y)
  
  data <- cbind(data, ID2=data$ID)
  levels(data$ID2) <-   gsub('.{2}$', '', levels(data$ID2))
  summary <- cbind(aggregate(data$blanked, by=list(data$ID2), FUN=mean), aggregate(data$blanked, by=list(data$ID2), FUN=sd)[2], pvalue=NA)
  names(summary)[1:3] <- c("ID", "mean", "sd")
  summary <- cbind(summary, t(as.data.frame(strsplit(as.character(summary$ID), " ")))[,-3])
  names(summary)[5:6] <- c("Strain", "Genotype")
  data <- cbind(data, t(as.data.frame(strsplit(as.character(data$ID), " ")))[,-3])
  names(data)[6:7] <- c("Strain", "Genotype")
  
  for(j in unique(summary$Strain)){
    for(i in unique(summary$Genotype)){
      summary$pvalue[which(summary$Strain == j & summary$Genotype == i)] <- t.test(data$blanked[which(data$Strain == j & data$Genotype == i)], data$blanked[which(data$Strain == j & data$Genotype == "control")])$p.value
    }}
  summary <- cbind(summary, signif=NA)
  summary$signif[which(summary$pvalue < 0.01)] <- "*"
  
  #calculate enzyme units
  #coefficient: 6220 M-1cm-1 or 6.22E-03
  #path length: 0.631cm
  e <- 6.22E-03 * 0.631
  #protein concentration: 20Âµg/mL in 200Âµl
  p <- 20/1000*200 #Âµg
  #slope unit: A340/s
  
  #final unit: mu/mg protein
  summary <- cbind(summary, unit_mean=(summary$mean*60)/e/p*1000, unit_sd=(summary$sd*60)/e/p*1000)
  
  write.table(summary, "Data/Enzyme_assays_p416Aco1_summary_v01_2021.csv", sep=",")
  }
  
  ### Aconitase: Iron
  {
    data <- read.csv("Data/150414_Aco_Iron_output.csv", stringsAsFactors = T)
    data <- cbind(data, ID=paste(gsub('.{11}$', '', as.character(data$StrainID)),data$Replicate))
    data <- merge(data[which(data$Substrate == "+citrate"),c(4,5,8,10,12,13)], data[which(data$Substrate == "-citrate"),c(4,5,12)], by=c("Strain", "Replicate"))
    data <- cbind(data, blanked=data$slope.x-data$slope.y)
    
    data <- cbind(data, ID2=data$ID)
    data$ID2 <- factor(data$ID2)
    levels(data$ID2) <-   gsub('.{2}$', '', levels(data$ID2))
    summary <- cbind(aggregate(data$blanked, by=list(data$ID2), FUN=mean), aggregate(data$blanked, by=list(data$ID2), FUN=sd)[2], pvalue=NA)
    names(summary)[1:3] <- c("ID", "mean", "sd")
    summary <- cbind(summary, t(as.data.frame(strsplit(as.character(summary$ID), " ")))[,c(4,6)])
    names(summary)[5:6] <- c("Deferoxamine", "Iron")
    
    for(i in unique(summary$ID)){
      summary$pvalue[which(summary$ID == i)] <- t.test(data$blanked[which(data$ID2 == i)], data$blanked[which(data$ID2 == "YSBN11 wt pRS313 w/o supplement 0")])$p.value
    }
    summary <- cbind(summary, signif=NA)
    summary$signif[which(summary$pvalue < 0.05)] <- "*"
    
    #calculate enzyme units
    #coefficient: 6220 M-1cm-1 or 6.22E-03
    #path length: 0.631cm
    e <- 6.22E-03 * 0.631
    #protein concentration: 20Âµg/mL in 200Âµl
    p <- 20/1000*200 #Âµg
    #slope unit: A340/s
    
    #final unit: mu/mg protein
    summary <- cbind(summary, unit_mean=(summary$mean*60)/e/p*1000, unit_sd=(summary$sd*60)/e/p*1000)
    
    data <- data %>%
      dplyr::select(ID2, concentration, Iron, blanked) %>%
      dplyr::rename(ID = ID2, Deferoxamine = concentration) %>%
      dplyr::mutate(activity = blanked*60/e/p*1000)
    
    
    write.table(summary, "Data/Enzyme_assays_Iron_summary_v03.csv", sep=",")
    write.table(data, "Data/Enzyme_assays_Iron_data.csv", sep=",")
    
  }
}

## Figure 1----

### Fig. 1C
{
  data <- read.csv("Data/SNP_enrichment.tsv", sep="\t")
  
  Figures[["1C"]] <- ggplot(data, aes(x = generation, y = value_cor, group=SNP))+
    geom_area(stat="identity", position="stack", aes(fill=SNP))+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position=c(0,1), 
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_fill_brewer(breaks=levels(data$SNP)[c(3,1,2)], labels=c(lb9, lb7, lb8), palette=palette)+
    ylab(bquote('allele frequency [%]'))
  
  ### Fig. 1E
  data <- read.csv("Data/Growcurves_YSBN11_130526_fit_values_summary.csv", stringsAsFactors = T)
  data$ATP3 <- factor(data$ATP3, levels=levels(data$ATP3)[c(5,1,2,3,4)])
  data$clone <- factor(data$clone, levels=levels(data$clone)[c(6,5,3,2,4,1)])
  dodge <- position_dodge(width=0.9)
  data$genotype <- factor(data$genotype, levels=levels(data$genotype)[c(3,2,1)])
  ttest <- read.csv("Data/Growcurves_YSBN11_130526_fit_values_ttest.csv")
  data <- merge(data, ttest[which(ttest$Reference == "init"),2:3], by="clone")
  data <- cbind(data, signif=NA)
  data$signif[which(data$pvalue < 0.05)] <- "*"
  
  
  ggplot(data, aes(x = clone, y = mean))+
    geom_bar(stat="identity", aes(fill=clone, group=clone), colour="black")+
    #scale_fill_discrete()+
    #scale_fill_grey(breaks=levels(data$genotype)[c(1,2)], labels=c(lb6,lb5))+
    scale_fill_manual(values=c("black", "red", brewer.pal(n = 8, name = palette)[1:4]), breaks=levels(data$clone), labels=c(lb6, lb5, lb27, lb28, lb29, lb30))+
    theme_classic()+
    #guides(fill=guide_legend(override.aes = list(colour=NULL)))+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,0.15))+ 
    scale_x_discrete(breaks=levels(data$clone), labels=c(lb6, lb5, lb27, lb28, lb29, lb30))+
    ylab(bquote('growth rate [OD min'^-1~']'))+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
    xlab("")+
    #geom_text(aes(label=round(pvalue,4), y = 0.02), position=dodge, angle=90, hjust=0, size=4)+
    geom_text(aes(label=signif, y=mean+sd+0.005), position=dodge, angle=90, hjust=0, vjust=0.8, size=10)
}

### Fig. 1E
{
  data <- read.csv("Data/Growcurves_YSBN11_130526_fit_values_summary.csv", stringsAsFactors = T)
  data$ATP3 <- factor(data$ATP3, levels=levels(data$ATP3)[c(5,1,2,3,4)])
  data$clone <- factor(data$clone, levels=levels(data$clone)[c(6,5,3,2,4,1)])
  dodge <- position_dodge(width=0.9)
  data$genotype <- factor(data$genotype, levels=levels(data$genotype)[c(3,2,1)])
  ttest <- read.csv("Data/Growcurves_YSBN11_130526_fit_values_ttest.csv")
  data <- merge(data, ttest[which(ttest$Reference == "init"),2:3], by="clone")
  data <- cbind(data, signif=NA)
  data$signif[which(data$pvalue < 0.05)] <- "*"
  
  idata <- read.delim("Data/Growcurves_YSBN11_130526_fit_values.csv", stringsAsFactors = T)
  
  cdata <- idata %>%
    dplyr::group_by(clone) %>%
    dplyr::summarise(count = length(mu.spline))
  
  Figures[["1E"]] <- ggplot(data, aes(x = clone, y = mean))+
    geom_bar(stat="identity", aes(fill=clone, group=clone), colour="black")+
    geom_text(aes(label=signif(pvalue,2), y=mean+sd+0.005), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
    geom_jitter(data = idata, aes(y = mu.spline), colour = "grey", width = 0.2)+
    scale_fill_manual(values=c("black", "red", brewer.pal(n = 8, name = palette)[1:4]), breaks=levels(data$clone), labels=c(lb6, lb5, lb27, lb28, lb29, lb30))+
    theme_classic()+
    #guides(fill=guide_legend(override.aes = list(colour=NULL)))+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,0.15))+ 
    scale_x_discrete(breaks=levels(data$clone), labels=c(lb6, lb5, lb27, lb28, lb29, lb30))+
    ylab(bquote('growth rate [OD min'^-1~']'))+
    xlab("")

}

### Fig. 1F
{
  data <- read.csv("Data/Growcurves_YSBN11_130914_fit_values_summary.csv")
  data$Strain <- factor(data$Strain)
  data$Strain <- factor(data$Strain, levels=levels(data$Strain)[c(6,1,2,5,3,4)])
  data$plasmid <- factor(data$plasmid)
  data$plasmid <- factor(data$plasmid, levels=c("pRS313", "pRS313_ATP3", "pRS313_ATP3-6", "pRS313_ATP3-7", "pRS313_ATP3-67"))
  dodge <- position_dodge(width=0.9)
  data$genotype <- factor(data$genotype)
  data$genotype <- factor(data$genotype, levels=levels(data$genotype)[2:1])
  ttest <- read.csv("Data/Growcurves_YSBN11_130914_fit_values_ttest.csv")
  data <- merge(data, ttest[which(ttest$Reference == "rho0_pRS313"),2:3], by="Strain")
  data <- cbind(data, signif=NA)
  data$signif[which(data$pvalue < 0.01)] <- "*"
  
  idata <- read.delim("Data/Growcurves_YSBN11_130914_fit_values.csv", stringsAsFactors = T)
  
  
  Figures[["1F"]] <- ggplot(data, aes(x = Strain, y = mean))+
    geom_bar(stat="identity", aes(fill=Strain, group=Strain), colour="black")+
    geom_text(aes(label=signif(pvalue, 2), y=mean+sd+0.005), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
    geom_jitter(data = idata, aes(y = mu.model), colour = "grey", width = 0.2)+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10,lb11,lb12,lb13,lb14,lb15), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(5,1,2)], "#73A5D4"))+
    theme_classic()+
    #guides(fill=guide_legend(override.aes = list(colour=NULL)))+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_y_continuous(expand = c(0,0), limits = c(0,0.25))+ 
    scale_x_discrete(breaks=levels(data$Strain), labels=c(lb10, lb11, lb12, lb13, lb14, lb15))+
    ylab(bquote('growth rate [OD min'^-1~']'))+
    xlab("")

}

### Fig. 1G
{
  data <- read.csv("Data/Growcurves_YSBN11_130914_summary_scarse.csv")
  fit <- read.csv("Data/Growcurves_YSBN11_130914_fit_summary.csv")
  
  Figures[["1G"]] <- ggplot(data, aes(x = time, y = mean, group=Strain))+
    geom_point(aes(colour=Strain))+
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, fill=Strain), alpha=0.2)+
    geom_line(data=fit, aes(colour=Strain), size=1)+
    xlim(0,25)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_linetype_discrete(guide=F)+
    scale_colour_manual(breaks=levels(data$Strain)[c(6,1,2,5,3,4)], labels=c(lb10,lb11,lb12,lb13,lb14,lb15), values=c("red",brewer.pal(n = 8, name = palette)[c(5,2)], "#73A5D4" , brewer.pal(n = 8, name = palette)[1], "black"))+
    scale_fill_manual(breaks=levels(data$Strain)[c(6,1,2,5,3,4)], labels=c(lb10,lb11,lb12,lb13,lb14,lb15), values=c("red",brewer.pal(n = 8, name = palette)[c(5,2)], "#73A5D4" , brewer.pal(n = 8, name = palette)[1], "black"))+
    ylab(bquote('biomass [OD'[600]~']'))+
    xlab("time [hours]")
}

## Figure 2----

### Fig. 2B
{
  data <- read.csv("Data/mtPotential_YSBN11_130805_summary.csv", stringsAsFactors = T)
  data <- data[which(data$Reference == "rho0 pRS313"),]
  data <- data[which(data$Strain %in% levels(data$Strain)[c(6,1,3,5)]),]
  data$Strain <- factor(data$Strain, levels=levels(data$Strain)[c(6,1,3,5)])
  data$plasmid <- factor(data$plasmid, levels=c("pRS313", "pRS313_ATP3-6", "pRS313_ATP3-7"))
  dodge <- position_dodge(width=0.9)
  data$genotype <- factor(data$genotype, levels=levels(data$genotype)[2:1])
  data <- cbind(data, signif=NA)
  data$signif[which(data$pvalue < 0.01)] <- "*"
  
  idata <- read.delim("Data/mtPotential_YSBN11_130805_raw.csv", sep = ",") %>%
    dplyr::filter(Strain %in% data$Strain)
  
  Figures[["2B"]] <- ggplot(data, aes(x = Strain, y = mean))+
    geom_bar(stat="identity", aes(fill=Strain, group=Strain), colour="black")+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
    geom_text(aes(label=signif(pvalue, 2), y=mean+se+0.01), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_text(aes(label = paste("n =", n)), y=0.05, angle=90, hjust=0, vjust=0.5, size=6, colour = "white")+
  geom_jitter(data = idata, aes(y = Pearsons.Correlation), colour = "grey", width = 0.2)+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    #coord_cartesian(ylim=c(0.4,0.86))+
    scale_x_discrete(breaks=levels(data$Strain), labels=c(lb10, lb11, lb13, lb14))+
    scale_y_continuous(expand = c(0,0), limits = c(0, 0.9))+
    ylab(bquote('Pearson Correlation Coefficient'))+
    xlab("")
}

### Fig. 2C
{
  data <- data.table::fread("Data/20210817_mtPotential_FACS.xls", stringsAsFactors = T) %>%
    dplyr::filter(StrainID %in% c("Rho+ pRS313", "Rho0 pRS313", "Rho0 ATP3 T911A", "Rho0 ATP3 G919C")) %>%
    dplyr::filter(Medium == "SM") %>%
    dplyr::filter(Replicate == 2) %>%
    dplyr::mutate(StrainID = factor(StrainID)) %>%
    dplyr::mutate(StrainID = factor(StrainID, levels = c("Rho+ pRS313", "Rho0 pRS313", "Rho0 ATP3 T911A", "Rho0 ATP3 G919C")))
  
  sdata <- data %>% 
    dplyr::group_by(StrainID, Medium, Replicate) %>%
    dplyr::summarise(median = median(value))
  
  Figures[["2C"]] <- ggplot(data, aes(x = value, y = ..scaled..))+
    #stat_density(data = data, aes(x=value, y=(..count..)/cnt*1000, color=StrainID), position="dodge", geom="line", size = 2)+
    #geom_density(aes(colour = StrainID, y = (..count..)/cnt*1000), rel_min_height = 0.01)+
    geom_density(aes(colour = StrainID, fill = StrainID), alpha = 0.2)+
    scale_x_log10(limits = c(100,1000), name = "DiOC6 Fluorescence")+
    scale_y_continuous(expand = c(0,0), name = "Relative Number of Cells") +
    facet_wrap(~ Replicate)+
    scale_colour_manual(labels=c(lb10b,lb11b,lb13b, lb14b), values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1,2)]))+
    scale_fill_manual(labels=c(lb10b,lb11b,lb13b, lb14b), values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1,2)]))+
    ggtitle("SM")+
    theme_classic()+
    theme(legend.position="bottom", 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.title=element_blank(),
          axis.line = element_line(),
          strip.background = element_blank(),
          legend.background=element_blank())
}

### Fig. 2D
{
  data <- read.xlsx("Data/ATP Assay.xlsx", sheetIndex=1, startRow = 3, endRow = 28)
  
  data <- data[,c(1,6,7,8)]
  data <- data[!is.na(data$strain),]
  names(data) <- c("StrainID", "Protein", "wo_OM", "w_OM")
  data$wo_OM <- as.numeric(as.character(data$wo_OM))
  data$w_OM <- as.numeric(as.character(data$w_OM))
  data$Protein <- as.numeric(as.character(data$Protein))
  
  data <- cbind(data, wo_OM_norm=data$wo_OM/data$Protein*1000, w_OM_norm=data$w_OM/data$Protein*1000)
  
  data <- cbind(data, Strain = t(as.data.frame(strsplit(as.character(data$StrainID), "_")))[,1])
  data <- cbind(data, wo_OM_activity = -data$wo_OM_norm*1000000/(6.22*1000000), w_OM_activity = -data$w_OM_norm*1000000/(6.22*1000000))
  data <- cbind(data, spec_activity=data$wo_OM_activity-data$w_OM_activity)
  data$Strain <- factor(data$Strain)
  data$Strain <- factor(data$Strain, levels = levels(data$Strain)[c(4,1,2,3)])
  
  
  melt <- melt(data[,c(7:9)], id="Strain")
  melt$Strain <- factor(melt$Strain, levels = levels(melt$Strain)[c(4,1,2,3)])
  levels(melt$variable) <- c("without oligomycin", "with 0.2 mg/mL oligomycin")
  
  summary <- cbind(aggregate(data$spec_activity, by=list(data$Strain), FUN=mean), aggregate(data$spec_activity, by=list(data$Strain), FUN=sd)[2], aggregate(data$spec_activity, by=list(data$Strain), FUN=se)[2])
  names(summary) <- c("Strain", "mean", "sd", "se")
  #ttests
  output <- data.frame()
  for(i in unique(data$Strain)){
    t <- t.test(data$spec_activity[which(data$Strain == i)] , data$spec_activity[which(data$Strain == "WT")])$p.value
    output <- rbind(output, data.frame(Strain=i, pvalue=t))
  }
  levels(output$Strain) <- c("wt", "rho0 G919C", "rho0 T911A", "rho0")
  
  sdata <- data %>% 
    dplyr::group_by(Strain) %>% 
    dplyr::summarize(mean = mean(spec_activity), sd = sd(spec_activity), se = se(spec_activity))
  
  signif <- data.frame()
  for(j in unique(data$Strain)){
    for(i in unique(data$Strain)){
      if(!j == i){
        signif <- rbind(signif, data.frame(Reference = j, Query = i, p.value = t.test(spec_activity ~ Strain, data = data[data$Strain %in% c(j, i),])$p.value))
        
      }
    }}
  
  ssignif <- signif %>%
    dplyr::filter(Reference == "WT") %>%
    dplyr::select(Query, p.value) %>%
    dplyr::rename(Strain = Query)
  sdata <- merge(sdata, ssignif, by = "Strain", all=T)
  
  sdata$signif <- NA
  sdata$signif[sdata$p.value  < 0.05] <- "*"
  sdata$p.value[is.na(sdata$p.value)] <- 1
  data$spec_activity[log2(data$spec_activity) < 0] <- 1
  
  Figures[["2D"]] <- ggplot(sdata, aes(x = Strain, y = log2(mean)))+
    #geom_bar(aes(fill = Strain), stat = "identity")
    geom_bar(stat="identity", aes(fill=Strain, group=Strain), colour="black")+
    geom_errorbar(aes(ymin=log2(mean-se), ymax=log2(mean+se)), width=0.2)+
    geom_text(aes(label=signif(p.value, 2), y=log2(mean+se+0.01)), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = data, aes(y = log2(spec_activity)), colour = "grey", width = 0.2)+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb6,lb5,lb27,lb28), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_x_discrete(breaks=levels(data$Strain), labels=c(lb6,lb5,lb27,lb28))+
    scale_y_continuous(expand = c(0,0), limits = c(0, 6))+
    ylab(bquote('Specific ATPase activity [?mol min '^-1~'mg'^-1~']'))+
    xlab("")
}

### Fig. 2E-H
{
  data <- read.csv("Data/fermentations_data.csv", stringsAsFactors = T)
  max <- read.csv("Data/fermentations_diauxic.csv", stringsAsFactors = T)
  gases <- read.csv("Data/fermentations_data_gases.csv", stringsAsFactors = T)
  data <- rbind(data, gases)
  data$data <- factor(data$data)
  data$compound <- factor(data$compound)
  levels(data$compound)[4] <- "biomass"
  data$variable <- factor(data$variable, levels=levels(data$variable)[c(1,2,4,3)])
  max$variable <- factor(max$variable, levels=levels(max$variable)[c(1,2,4,3)])
  
  # faceting according to compounds
  {
    i <- "biomass"
    Figures[["2E"]] <- ggplot(data[which(data$compound == i),], aes(x = timepoint, y = value, group=variable))+
      geom_line(aes(colour=variable), size=1.5, linetype=1)+
      geom_point(aes(colour=variable), size=3)+
      geom_vline(data=max, aes(xintercept=timepoint, colour = variable), linetype=2, size=0.5)+
      theme_classic()+
      theme(legend.justification=c(1,1), 
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            legend.position="none", 
            legend.direction="vertical",
            legend.title=element_blank(),
            #legend.key=element_rect(colour = "black"),
            legend.background=element_blank(),
            strip.background = element_blank())+
      scale_x_continuous("time [hrs]")+
      scale_y_continuous(paste0(i, " (", unique(data$unit[data$compound == i]), ")"), limits = c(0,max(data$value[data$compound == i])))+
      scale_colour_manual(labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
      annotate("text",x = 0, y = max(data$value[data$compound == i]), label=i, parse=F, hjust=0)
    
    i <- "glucose"
    Figures[["2F"]] <- ggplot(data[which(data$compound == i),], aes(x = timepoint, y = value, group=variable))+
      geom_line(aes(colour=variable), size=1.5, linetype=1)+
      geom_point(aes(colour=variable), size=3)+
      geom_vline(data=max, aes(xintercept=timepoint, colour = variable), linetype=2, size=0.5)+
      theme_classic()+
      theme(legend.justification=c(1,1), 
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            legend.position="none", 
            legend.direction="vertical",
            legend.title=element_blank(),
            #legend.key=element_rect(colour = "black"),
            legend.background=element_blank(),
            strip.background = element_blank())+
      scale_x_continuous("time [hrs]")+
      scale_y_continuous(paste0(i, " (", unique(data$unit[data$compound == i]), ")"), limits = c(0,max(data$value[data$compound == i])))+
      scale_colour_manual(labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
      annotate("text",x = 0, y = max(data$value[data$compound == i]), label=i, parse=F, hjust=0)
    
    i <- "ethanol"
    Figures[["2G"]] <- ggplot(data[which(data$compound == i),], aes(x = timepoint, y = value, group=variable))+
      geom_line(aes(colour=variable), size=1.5, linetype=1)+
      geom_point(aes(colour=variable), size=3)+
      geom_vline(data=max, aes(xintercept=timepoint, colour = variable), linetype=2, size=0.5)+
      theme_classic()+
      theme(legend.justification=c(1,1), 
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            legend.position="none", 
            legend.direction="vertical",
            legend.title=element_blank(),
            #legend.key=element_rect(colour = "black"),
            legend.background=element_blank(),
            strip.background = element_blank())+
      scale_x_continuous("time [hrs]")+
      scale_y_continuous(paste0(i, " (", unique(data$unit[data$compound == i]), ")"), limits = c(0,max(data$value[data$compound == i])))+
      scale_colour_manual(labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
      annotate("text",x = 0, y = max(data$value[data$compound == i]), label=i, parse=F, hjust=0)
    
    i <- "glycerol"
    Figures[["2H"]] <- ggplot(data[which(data$compound == i),], aes(x = timepoint, y = value, group=variable))+
      geom_line(aes(colour=variable), size=1.5, linetype=1)+
      geom_point(aes(colour=variable), size=3)+
      geom_vline(data=max, aes(xintercept=timepoint, colour = variable), linetype=2, size=0.5)+
      theme_classic()+
      theme(legend.justification=c(1,1), 
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            legend.position="none", 
            legend.direction="vertical",
            legend.title=element_blank(),
            #legend.key=element_rect(colour = "black"),
            legend.background=element_blank(),
            strip.background = element_blank())+
      scale_x_continuous("time [hrs]")+
      scale_y_continuous(paste0(i, " (", unique(data$unit[data$compound == i]), ")"), limits = c(0,max(data$value[data$compound == i])))+
      scale_colour_manual(labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
      annotate("text",x = 0, y = max(data$value[data$compound == i]), label=i, parse=F, hjust=0)
    
    # #pdf("data/fermentations_data.pdf", width = 9, height = 6)
    # for(i in  unique(data$compound)){
    #   
    #   p <- ggplot(data[which(data$compound == i),], aes(x = timepoint, y = value, group=variable))+
    #     geom_line(aes(colour=variable), size=1.5, linetype=1)+
    #     geom_point(aes(colour=variable), size=3)+
    #     geom_vline(data=max, aes(xintercept=timepoint, colour = variable), linetype=2, size=0.5)+
    #     theme_classic()+
    #     theme(legend.justification=c(1,1), 
    #           legend.position="none", 
    #           legend.direction="vertical",
    #           legend.title=element_blank(),
    #           #legend.key=element_rect(colour = "black"),
    #           legend.background=element_blank(),
    #           strip.background = element_blank())+
    #     scale_x_continuous("time [hrs]")+
    #     scale_y_continuous(paste0(i, " (", unique(data$unit[data$compound == i]), ")"), limits = c(0,max(data$value[data$compound == i])))+
    #     scale_colour_manual(labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
    #     annotate("text",x = 0, y = max(data$value[data$compound == i]), label=i, parse=F, hjust=0)
    #   print(p)
    # }
    #dev.off()
  }
}

## Figure 3----

### Fig. 3A
{
  t <- list.dirs("Data/")
  t <- t[grepl("String", t)]
  
  output <- data.frame()
  
  for(i in t[2:4]){
    u <-list.files(i, full.names = T)
    u <- u[grepl("enrichment", u)]
    exp <- sub("-db_", "", str_extract(i, "-db_.+"))
    for(j in u){
      tdf <- read.delim(j)
      group <- unlist(strsplit(j, "/", fixed=T))
      group <- unlist(strsplit(group[length(group)], ".", fixed=T))
      output <- rbind(output, cbind(tdf, Group = group[2], Experiment = exp))
    }
  }
  
  t <- output[output$Experiment == "WT-Rho0" & output$Group == "RCTM",]
  t$X.term.ID <- factor(t$X.term.ID, levels = unique(t$X.term.ID))
  t$term.description <- factor(t$term.description, levels = rev(unique(t$term.description)))
  t$Label <- paste(t$term.description, "\n(", t$X.term.ID, ")", sep="")
  t$Label <- factor(t$Label, levels = rev(unique(t$Label)))
  t$Enrichment <- 1/t$background.gene.count*t$observed.gene.count
  t$pValue <- as.character(round(t$false.discovery.rate, 3))
  t$pValue[t$pValue == "0"] <- formatC(t$false.discovery.rate[t$pValue == "0"], format = "e", digits = 1)
  
  
  Figures[["3A"]] <- ggplot(t, aes(x = Label, y = -log10(false.discovery.rate)))+
    geom_bar(stat = "identity", aes(fill = Enrichment), colour = "black")+
    geom_text(aes(label = pValue), y = 0.1, hjust = 0)+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_gradient(low = "white", high = "#FF0000")+
    labs(x = "Reactome Pathway", y = "-log10 enrichment p-value")+
    coord_flip()+
    theme_minimal()+
    theme(legend.position = c(1,0),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.ticks.x = element_line(),
          panel.grid = element_blank(),
          #axis.line.y = element_line(),
          axis.line.x = element_line(),
          legend.justification = c(1,0))
}

### Fig. 3B
{
  fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
  outfile <- try(readLines("Data/Metabolites_map_template_v2.svg"))
  poi <- fcsummary %>%
    filter(Strain %in% c("rho+_313", "rho0_313", "rho0_ATP3_T911A"))
  poi$Strain <- factor(poi$Strain)
  poi$Strain <- factor(poi$Strain, levels = levels(poi$Strain))
  
  # make code
  cde <-
    data.frame(
      Compound = unique(odata$Compound),
      Experiment = "proteomics",
      Key = capitalize_str(tolower(paste0(
        sub(";.+", "", unique(odata$Compound)), "p"
      ))),
      style = "st2"
    )
  
  cde$Key[cde$Compound == "YGR012W"] <- "Mcy1p"
  code <- read.delim("Data/Metabolites_map_code.xls")
  code <- rbind(code, cde)
  
  maxval <- 3
  colrange <- seq(-maxval, maxval, by= 0.01)
  coloring_ <- data.frame(value=round(colrange,2), color=bluered(length(colrange)))
  coloring <- data.frame(value=round(colrange, 2))
  coloring <- merge(coloring, coloring_, by="value", all.x=T)
  coloring$color[is.na(coloring$color) & coloring$value < 0] <- "#0000FF"
  coloring$color[is.na(coloring$color) & coloring$value > 0] <- "#FF0000"
  
  ColorKey <- ggplot(coloring, aes(x = value, y = 1))+
    geom_bar(aes(fill = factor(value)), stat = "identity", width = 0.1)+
    xlab("log2 fold change")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.ticks.x = element_line(colour = "black",size = 1, linetype = "solid"))+
    scale_fill_manual(values = as.character(coloring$color))+
    scale_y_continuous(expand = c(0,0), limits = c(0,1))
  
  for(i in unique(poi$Compound)){
    
    spoi <- gsub(",", "", gsub(".", "", i, fixed=T), fixed=T)
    key <- code$Key[code$Compound == i]
    styletem <- code$style[code$Compound == i]
    
    
    if(i %in% poi$Compound[!is.na(poi$mean)] & length(styletem) > 0){
      # get data
      t <- poi[poi$Compound == i & !is.na(poi$mean),] %>%
        dplyr::group_by(Compound, Strain) %>%
        dplyr::summarise(value = mean(mean), .groups = "drop") %>%
        arrange(rev(Strain))
      
      t$value[log2(t$value) < min(coloring$value)] <- 2^min(coloring$value)
      t$value[log2(t$value) > max(coloring$value)] <- 2^max(coloring$value)
      
      # construct styles
      style <- list()
      style[[1]] <- paste("	.st", spoi,"evo{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho0_ATP3_T911A"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      style[[2]] <- paste("	.st", spoi,"rho{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho0_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      style[[3]] <- paste("	.st", spoi,"wt{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho+_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      
      names(style) <- c(paste0("st", spoi,"evo"), paste0("st", spoi,"rho"), paste0("st", spoi,"wt"))
      
      # add styles to style table
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[1]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[2]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[3]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      
      # apply style to element
      
      for(j in unique(which(grepl(paste0(key, "<"), outfile)))){
        rng <- which(grepl(styletem, outfile, fixed=T))
        rng <- rng[rng > j][1:3]
        
        for(k in 1:3){
          outfile[rng[k]] <- sub(styletem, names(style)[k], outfile[rng[k]])
        }
      }
    } else {
      print(i)
    }
  }
  outfile <- outfile[!is.na(outfile)]
  
  writeLines(text=outfile, con="Plots/Figure_3A.svg", sep="\n")
  
}

### Fig. 3C
{
  data <- read.csv("Data/150917_growth_outliersrm.csv", stringsAsFactors = T)
  data_raw <- read.csv("Data/150917_growth_outliersrm_raw_defect.csv", stringsAsFactors = T)
  data$StrainID <- factor(data$StrainID, levels=levels(data$StrainID)[c(3,1,2)])
  data$Mix <- factor(data$Mix, levels=levels(data$Mix)[c(1,6,3,11,5,7,12,8,10,4,9,2)])
  data$Mix <- factor(data$Mix, levels=rev(levels(data$Mix)))
  data_raw$StrainID <- factor(data_raw$StrainID)
  data_raw$StrainID <- factor(data_raw$StrainID, levels = levels(data$StrainID))
  
  Figures[["3C"]] <- ggplot(data, aes(x = Mix, y = mu.spline_mean, group=StrainID, fill=StrainID))+
    geom_bar(stat="identity", position=dodge, colour="black")+
    geom_errorbar(aes(ymin=mu.spline_mean-mu.spline_sd, ymax=mu.spline_mean+mu.spline_sd), width=0.2, position=dodge)+
    geom_point(data = data_raw, aes(y = mu.spline), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]),labels=c(lb10b,lb11b,lb13b))+
    theme_classic()+
    #geom_hline(yintercept = data$mu.spline_mean[which(data$Mix %in% c("+QERL", "-AA") & data$StrainID == "YSBN11 rho0 313")], linetype=2, colour=rho)+
    #geom_hline(yintercept = data$defect_spline_mean_sub[which(data$Mix %in% c("+QERL", "-AA") & data$StrainID == "YSBN11 rho0 ATP3 T911A")], linetype=2, colour=brewer.pal(n = 8, name = palette)[c(1)])+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(1,1), 
          legend.direction="vertical",
          legend.background = element_blank(),
          legend.title=element_blank())+
    #coord_cartesian(ylim=c(0,0.210))+
    coord_flip()+
    scale_y_continuous(expand = c(0,0))+ 
    xlab("amino acid supplementation")+
    ylab(bquote('growth rate [OD min'^-1~']'))
}

### Fig. 3D
{
  data <- read.csv("Data/150917_growth_outliersrm.csv", stringsAsFactors = T)
  data_raw <- read.csv("Data/150917_growth_outliersrm_raw_defect.csv", stringsAsFactors = T) %>%
    dplyr::filter(Mix %in% c("-AA", "+QERL"))
  data$StrainID <- factor(data$StrainID, levels=levels(data$StrainID)[c(3,1,2)])
  data$Mix <- factor(data$Mix, levels=levels(data$Mix)[c(1,6,3,11,5,7,12,8,10,4,9,2)])
  data_raw$StrainID <- factor(data_raw$StrainID)
  data_raw$StrainID <- factor(data_raw$StrainID, levels = levels(data_raw$StrainID)[c(3,1,2)])
  
  
  Figures[["3D"]] <- ggplot(data[which(data$Mix %in% c("-AA", "+QERL")),], aes(x = Mix, y = defect_spline_mean_sub, group=StrainID, fill = StrainID))+
    geom_bar(aes(fill=StrainID), stat="identity", position=dodge, colour="black")+
    geom_errorbar(aes(ymin=defect_spline_mean_sub-defect_spline_sd, ymax=defect_spline_mean_sub+defect_spline_sd), width=0.2, position=dodge)+
    geom_point(data = data_raw, aes(y = defect_spline_sub), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]),labels=c(lb10b,lb11b,lb13b))+
    theme_classic()+
    geom_hline(yintercept = data$defect_spline_mean_sub[which(data$Mix %in% c("+QERL", "-AA") & data$StrainID == "YSBN11 rho0 313")], linetype=2, colour=rho)+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(1,1), 
          legend.direction="vertical",
          legend.background = element_blank(),
          legend.title=element_blank())+
    #coord_cartesian(ylim=c(0,0.210))+
    scale_y_continuous(expand = c(0,0))+ 
    xlab("amino acid supplementation")+
    ylab(bquote('growth rate defect [%]'))
}

## Figure 4----

### Fig. 4A
{
  # svg
  {
    odata <- read.delim("Data/20210816_Proteomics_General.tsv") 
    fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
    
    # make code
    cde <-
      data.frame(
        Compound = unique(odata$Compound),
        Experiment = "proteomics",
        Key = capitalize_str(tolower(paste0(
          sub(";.+", "", unique(odata$Compound)), "p"
        ))),
        style = "st16"
      )
    cde$Key[cde$Compound == "YGR012W"] <- "Mcy1p"
    
    outfile <- try(readLines("Data/TCA_map_template_prot-01.svg"))
    poi <- fcsummary %>%
      filter(Strain %in% c("rho+_313", "rho0_313", "rho0_ATP3_T911A"))
    poi$Strain <- factor(poi$Strain)
    poi$Strain <- factor(poi$Strain, levels = levels(poi$Strain))
    
    code <- read.delim("Data/TCA_map_code.xls")
    code <- rbind(code, cde)
    
    maxval <- 3
    colrange <- seq(-maxval, maxval, by= 0.01)
    coloring_ <- data.frame(value=round(colrange,2), color=bluered(length(colrange)))
    coloring <- data.frame(value=round(colrange, 2))
    coloring <- merge(coloring, coloring_, by="value", all.x=T)
    coloring$color[is.na(coloring$color) & coloring$value < 0] <- "#0000FF"
    coloring$color[is.na(coloring$color) & coloring$value > 0] <- "#FF0000"
    
    ColorKey <- ggplot(coloring, aes(x = value, y = 1))+
      geom_bar(aes(fill = factor(value)), stat = "identity", width = 0.1)+
      xlab("log2 fold change")+
      theme_minimal()+
      theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 12),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 12),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.ticks.x = element_line(colour = "black",size = 1, linetype = "solid"))+
      scale_fill_manual(values = as.character(coloring$color))+
      scale_y_continuous(expand = c(0,0), limits = c(0,1))
    
    pdf("Data/TCA_map_key.pdf")
    print(ColorKey)
    dev.off()
    
    
    for(i in unique(poi$Compound)){
      
      spoi <- gsub(",", "", gsub(".", "", i, fixed=T), fixed=T)
      key <- code$Key[code$Compound == i]
      styletem <- code$style[code$Compound == i]
      
      
      if(i %in% poi$Compound[!is.na(poi$mean)] & length(styletem) > 0){
        # get data
        t <- poi[poi$Compound == i & !is.na(poi$mean),] %>%
          dplyr::group_by(Compound, Strain) %>%
          dplyr::summarise(value = mean(mean), .groups = "drop") %>%
          arrange(rev(Strain))
        
        t$value[log2(t$value) < min(coloring$value)] <- 2^min(coloring$value)
        t$value[log2(t$value) > max(coloring$value)] <- 2^max(coloring$value)
        
        # construct styles
        style <- list()
        style[[1]] <- paste("	.st", spoi,"evo{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho0_ATP3_T911A"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
        style[[2]] <- paste("	.st", spoi,"rho{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho0_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
        style[[3]] <- paste("	.st", spoi,"wt{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho+_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
        
        names(style) <- c(paste0("st", spoi,"evo"), paste0("st", spoi,"rho"), paste0("st", spoi,"wt"))
        
        # add styles to style table
        outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[1]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
        outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[2]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
        outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[3]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
        
        # apply style to element
        
        for(j in unique(which(grepl(paste0(key, "<"), outfile)))){
          rng <- which(grepl(styletem, outfile, fixed=T))
          rng <- rng[rng > j][1:3]
          
          for(k in 1:3){
            outfile[rng[k]] <- sub(styletem, names(style)[k], outfile[rng[k]])
          }
        }
      } else {
        print(i)
      }
    }
    outfile <- outfile[!is.na(outfile)]
    
    writeLines(text=outfile, con="Plots/Figure_4A.svg", sep="\n")
    
  }
}

### Fig. 4B
{
  # heat maps
  {
  fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
  
  for(i in unique(fcsummary$Experiment)){
    exp <- i
    
    mdata <- fcsummary %>%
      filter(Strain %in% c("rho+_313", "rho0_313", "rho0_ATP3_T911A") & Experiment == !!exp) %>%
      select(Compound, Strain, mean) %>%
      tidyr::pivot_wider(names_from = "Strain", values_from = "mean") %>%
      tibble::column_to_rownames(var = "Compound") %>%
      as.matrix(.) %>%
      log2(.)
    
    heatmap.2(mdata, main = exp,
              cexCol=0.5, 
              Colv=F,
              Rowv = F,
              trace="none", 
              breaks=mybreakslog, 
              col=mycollog, 
              cexRow=0.5, 
              na.color="grey",
              keysize=1)
  }
  }
  
  # svg
  {
    odata <- read.delim("Data/20210816_Proteomics_General.tsv") 
    fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
    # make code
    cde <-
      data.frame(
        Compound = unique(odata$Compound),
        Experiment = "proteomics",
        Key = capitalize_str(tolower(paste0(
          sub(";.+", "", unique(odata$Compound)), "p"
        ))),
        style = "st16"
      )
    cde$Key[cde$Compound == "YGR012W"] <- "Mcy1p"
    
    
    outfile <- try(readLines("Data/TCA_map_template_met-01.svg"))
    poi <- fcsummary %>%
      filter(Strain %in% c("rho+_313", "rho0_313", "rho0_ATP3_T911A"))
    poi$Strain <- factor(poi$Strain)
    poi$Strain <- factor(poi$Strain, levels = levels(poi$Strain))
    
    
    # code <- fcdata %>%
    #   distinct(Experiment, Compound) %>%
    #   mutate(Key = NA) %>%
    #   write.table(., "C:/Dropbox (Personal)/AG Ralser_JV/Publications/rho0_project/Figures/data/Metabolites_map_code.xls", sep="\t", row.names=F)
    
    code <- read.delim("Data/TCA_map_code.xls")
    code <- rbind(code, cde)
    
    maxval <- 3
    colrange <- seq(-maxval, maxval, by= 0.01)
    coloring_ <- data.frame(value=round(colrange,2), color=bluered(length(colrange)))
    coloring <- data.frame(value=round(colrange, 2))
    coloring <- merge(coloring, coloring_, by="value", all.x=T)
    coloring$color[is.na(coloring$color) & coloring$value < 0] <- "#0000FF"
    coloring$color[is.na(coloring$color) & coloring$value > 0] <- "#FF0000"
    
    ColorKey <- ggplot(coloring, aes(x = value, y = 1))+
      geom_bar(aes(fill = factor(value)), stat = "identity", width = 0.1)+
      xlab("log2 fold change")+
      theme_minimal()+
      theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 12),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 12),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.ticks.x = element_line(colour = "black",size = 1, linetype = "solid"))+
      scale_fill_manual(values = as.character(coloring$color))+
      scale_y_continuous(expand = c(0,0), limits = c(0,1))
    
    for(i in unique(poi$Compound[poi$Experiment %in% c("TCA", "PPP")])){
      
      spoi <- gsub(",", "", gsub(".", "", i, fixed=T), fixed=T)
      key <- code$Key[code$Compound == i]
      styletem <- code$style[code$Compound == i]
      
      
      if(i %in% poi$Compound[!is.na(poi$mean)] & length(styletem) > 0){
        # get data
        t <- poi[poi$Compound == i & !is.na(poi$mean),] %>%
          dplyr::group_by(Compound, Strain) %>%
          dplyr::summarise(value = mean(mean), .groups = "drop") %>%
          arrange(rev(Strain))
        
        t$value[log2(t$value) < min(coloring$value)] <- 2^min(coloring$value)
        t$value[log2(t$value) > max(coloring$value)] <- 2^max(coloring$value)
        
        # construct styles
        style <- list()
        style[[1]] <- paste("	.st", spoi,"evo{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho0_ATP3_T911A"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
        style[[2]] <- paste("	.st", spoi,"rho{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho0_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
        style[[3]] <- paste("	.st", spoi,"wt{fill:", coloring$color[coloring$value == round(log2(t$value[t$Strain == "rho+_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
        
        names(style) <- c(paste0("st", spoi,"evo"), paste0("st", spoi,"rho"), paste0("st", spoi,"wt"))
        
        # add styles to style table
        outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[1]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
        outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[2]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
        outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[3]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
        
        # apply style to element
        
        for(j in unique(which(grepl(paste0(key, "<"), outfile)))){
          if(j > 240){
            rng <- which(grepl(styletem, outfile, fixed=T))
            rng <- rng[rng > j][1:3]
            
            for(k in 1:3){
              outfile[rng[k]] <- sub(styletem, names(style)[k], outfile[rng[k]])
            }
          }}
      } else {
        print(i)
      }
    }
    
    outfile <- outfile[!is.na(outfile)]
    
    writeLines(text=outfile, con="Plots/Figure_4B.svg", sep="\n")
  }
}

### Fig. 4C
{
  fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
  signif <- read.delim("Data/20210816_CombinedData_Signif.tsv")
  poi <- c("ACO1", "IDH1", "IDP1")
  poidata <- fcsummary[fcsummary$Compound %in% poi,]
  poisdata <- fcdata %>%
    filter(Strain %in% c("rho+_313", "rho0_313", "rho0_ATP3_T911A")) %>%
    filter(Compound %in% poi)
  
  poidata <- merge(poidata, signif, by=c("Compound", "Strain"), all.x=T)
  poidata$label <- round(poidata$p.value, 5)
  
  
  Figures[["4C"]] <- ggplot(poidata, aes(x = Compound, y = mean*100, fill=Strain))+
    geom_bar(stat="identity", position=dodge,aes(fill=Strain, group=Strain), colour="black")+
    geom_errorbar(aes(ymin=mean*100-sd*100, ymax=mean*100+sd*100), width=0.2, position=dodge)+
    geom_point(data = poisdata, aes(y = foldchange*100), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.15))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]), labels=c(lb10b,lb11b,lb13b))+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.background = element_blank(),
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,550))+ 
    xlab("")+
    geom_text(aes(label=signif(p.value, 2), y=mean*100+sd*100+5), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    ylab("protein abundance [%]")+
    xlab("")
}

### Fig. 4D
{
  #ACO
  data <- read.csv("Data/150207_Aco_fit.csv")
  data <- data[which(data$Cofactor == "Aconitase"),]
  data <- cbind(data, ID=paste(data$Genotype, data$Plasmid, data$Replicate))
  data <- merge(data[which(data$Isocitrate == "+citrate"),c(4,8,11)], data[which(data$Isocitrate == "-citrate"),c(4,8,11)], by=c("ID", "Cofactor"))
  data <- cbind(data, blanked=data$slope.x-data$slope.y)
  data <- cbind(data, ID2=data$ID)
  data$ID2 <-   gsub('.{2}$', '', data$ID2)
  data$ID2 <- sub("ATP3", "ATP3T911A", data$ID2)
  
  data_aco <- data
  
  #IDH/NADP
  data <- read.csv("Data/150130_NAD_NADP_fit.csv")
  data <- cbind(data, ID=paste(data$Genotype, data$Plasmid, data$Replicate))
  data <- merge(data[which(data$Isocitrate == "+isocitrate"),c(4,8,11)], data[which(data$Isocitrate == "-isocitrate"),c(4,8,11)], by=c("ID", "Cofactor"))
  data <- cbind(data, blanked=data$slope.x-data$slope.y)
  
  data <- cbind(data, ID2=data$ID)
  data$ID2 <-   gsub('.{2}$', '', data$ID2)
  data$ID2 <- sub("ATP3", "ATP3T911A", data$ID2)
  data_idh <- data
  
  #MDH
  data <- read.csv("Data/150305_Aco_MDH_fit.csv")
  data <- data[which(data$Cofactor == "MDH"),]
  data <- cbind(data, ID=paste(data$Genotype, data$Plasmid, data$Replicate))
  data <- merge(data[which(data$Isocitrate == "+malate"),c(4,8,11)], data[which(data$Isocitrate == "-malate"),c(4,8,11)], by=c("ID", "Cofactor"))
  data <- cbind(data, blanked=data$slope.x-data$slope.y)
  
  data <- cbind(data, ID2=data$ID)
  data$ID2 <-   gsub('.{2}$', '', data$ID2)
  data_mdh <- data
  
  edata <- rbind(data_aco, data_idh, data_mdh)
  ndata <- edata %>% filter(ID2 == "wt pRS313") %>% group_by(Cofactor) %>% dplyr::summarise(norm = mean(blanked))
  
  edata <- edata %>%
    left_join(ndata, by = "Cofactor") %>%
    filter(!Cofactor == "MDH") %>%
    mutate(norm_value=100/norm*blanked) %>%
    mutate(Enzyme = factor(Cofactor, levels = c("Aconitase", "NAD", "NADP"))) %>%
    mutate(Strain = factor(ID2, levels = c("wt pRS313", "rho0 pRS313", "rho0 pRS313_ATP3T911A")))
  
  sedata <- edata %>%
    group_by(Strain, Enzyme) %>%
    dplyr::summarise(norm_mean = mean(norm_value), norm_sd = sd(norm_value))
  
  signif <- data.frame()
  for(i in unique(sedata$Enzyme)){
    t <- edata[edata$Enzyme == i,]
    for(j in unique(t$Strain)){
      signif <- rbind(signif, data.frame(Enzyme = i, Strain = j, p.value = t.test(t$norm_value[t$Strain == j], t$norm_value[t$Strain == "wt pRS313"], var.equal = T)$p.value))
      
    }
  }
  signif$signif <- F
  signif$signif[signif$p.value < 0.01] <- T
  
  
  sedata <- merge(sedata, signif, by=c("Enzyme", "Strain"), all.x=T)
  sedata$label <- round(sedata$p.value, 5)
  
  Figures[["4D"]] <- ggplot(sedata, aes(x = Enzyme, y = norm_mean, fill=Strain))+
    geom_bar(stat="identity", position=dodge,aes(fill=Strain, group=Strain), colour="black")+
    geom_errorbar(aes(ymin=norm_mean-norm_sd, ymax=norm_mean+norm_sd), width=0.2, position=dodge)+
    geom_point(data = edata, aes(y = norm_value), colour = "grey", position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]), breaks=levels(sedata$Strain), labels=c(lb10b,lb11b,lb13b))+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="vertical",
          legend.background = element_blank(),
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(-25,240))+ 
    scale_x_discrete(breaks=levels(sedata$Enzyme), labels=c("Aconitase", lb32, lb33))+
    xlab("")+
    geom_text(aes(label=signif(p.value, 2), y=norm_mean+norm_sd+5), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    ylab("enzyme activity [%]")+
    xlab("")
}

### Fig. 4F
{
  file <- "Data/labeling_data.xlsx" 
  sheet <- "slimmed"
  data <- readxl::read_xlsx(file,sheet = sheet,col_types = "text")
  
  data2 <- data %>% 
    mutate(labeling_percent = as.numeric(labeling_percent)) %>%
    mutate(time = as.numeric(time)) %>%
    mutate(abundance_nmol = as.numeric(abundance_nmol))
  
  temp <- data2 %>% filter(tracer == "glucose")
  
  temp <- temp %>% filter(metabolite %in% c("pyruvate", "citrate", "cis-aconitate", "glutamate", "gaba", "glutamine","succinate")) # "a-ketoglutarate", "fumarate", "malate"
  
  #remove outlier by sample (metabolites clearly affected in one sample)
  temp <- temp %>% filter(!(tracer == "glucose" & time == 30 & genotype == "rho0 313" & replicate == 2))
  
  temp2 <- temp %>% group_by(tracer, metabolite, time, genotype) %>%
    dplyr::mutate(mean_labeling_percent = mean(labeling_percent, na.rm = T)) %>%
    dplyr::mutate(sd_labeling_percent = sd(labeling_percent, na.rm = T)) %>%
    
    dplyr::mutate(mean_abundance_nmol = mean(abundance_nmol, na.rm = T)) %>%
    dplyr::mutate(sd_abundance_nmol = sd(abundance_nmol, na.rm = T)) %>%
    ungroup() 
  
  temp2$genotype <- factor(temp2$genotype, levels = c("wt 313", "rho0 313", "rho0 T911A"))
  
  # plotting single timepoints
  time_filter <- 30
  temp3 <- temp2 %>% filter(time == time_filter)
  
  temp3$genotype <- factor(temp3$genotype, levels = c("wt 313", "rho0 313", "rho0 T911A"))
  temp3$metabolite <- factor(temp3$metabolite, levels = metabolites)
  
  # labelling - facet wrap
  Figures[["4F"]] <- temp3 %>% 
    ggplot(aes(x=genotype, y=mean_labeling_percent, fill = genotype)) +
    geom_point(aes(x=genotype, y=labeling_percent), inherit.aes = T, colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.15)) + 
    geom_bar(stat="identity", position=position_dodge(width = 0.9), alpha = 1, colour = "black") + 
    geom_errorbar(aes(ymin=mean_labeling_percent - sd_labeling_percent, ymax=mean_labeling_percent + sd_labeling_percent), width=0.2) + 
    facet_wrap(~ metabolite, scales = "free_y", ncol = length(unique(temp2$metabolite))) +
    scale_colour_manual(values = c("wt 313" = wt, "rho0 313" = rho, "rho0 T911A" = brewer.pal(n = 8, name = palette)[c(1)])) +
    scale_fill_manual(values = c("wt 313" = wt, "rho0 313" = rho, "rho0 T911A" = brewer.pal(n = 8, name = palette)[c(1)])) +
    scale_y_continuous(expand = c(0,0))+ 
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="vertical",
          legend.background = element_blank(),
          strip.background = element_blank(),
          legend.title=element_blank())+
    ylab("metabolite labelling (%)")+
    xlab("")
}

### Fig. 4G
{
  file <- "Data/labeling_data.xlsx" 
  sheet <- "slimmed"
  data <- readxl::read_xlsx(file,sheet = sheet,col_types = "text")
  
  data2 <- data %>% 
    mutate(labeling_percent = as.numeric(labeling_percent)) %>%
    mutate(time = as.numeric(time)) %>%
    mutate(abundance_nmol = as.numeric(abundance_nmol))
  
  temp <- data2 %>% filter(tracer == "glucose")
  
  temp <- temp %>% filter(metabolite %in% c("pyruvate", "citrate", "cis-aconitate", "glutamate", "gaba", "glutamine","succinate")) # "a-ketoglutarate", "fumarate", "malate"
  
  #remove outlier by sample (metabolites clearly affected in one sample)
  temp <- temp %>% filter(!(tracer == "glucose" & time == 30 & genotype == "rho0 313" & replicate == 2))
  
  # remove outlier by metabolite
  # temp <- temp %>% filter(outlier != 1)
  
  temp2 <- temp %>% group_by(tracer, metabolite, time, genotype) %>%
    dplyr::mutate(mean_labeling_percent = mean(labeling_percent, na.rm = T)) %>%
    dplyr::mutate(sd_labeling_percent = sd(labeling_percent, na.rm = T)) %>%
    
    dplyr::mutate(mean_abundance_nmol = mean(abundance_nmol, na.rm = T)) %>%
    dplyr::mutate(sd_abundance_nmol = sd(abundance_nmol, na.rm = T)) %>%
    ungroup() 
  
  temp2$genotype <- factor(temp2$genotype, levels = c("wt 313", "rho0 313", "rho0 T911A"))
  
  # plotting single timepoints
  time_filter <- 30
  temp3 <- temp2 %>% filter(time == time_filter)
  
  temp3$genotype <- factor(temp3$genotype, levels = c("wt 313", "rho0 313", "rho0 T911A"))
  temp3$metabolite <- factor(temp3$metabolite, levels = metabolites)
  
  # abundance - facet wrap
  Figures[["4G"]] <- temp3 %>% 
    ggplot(aes(x=genotype, y=mean_abundance_nmol, fill = genotype)) +
    geom_point(aes(x=genotype, y=abundance_nmol), inherit.aes = T, colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.15)) + 
    geom_bar(stat="identity", position=position_dodge(width = 0.9), alpha = 1, colour = "black") + 
    geom_errorbar(aes(ymin=mean_abundance_nmol - sd_abundance_nmol, ymax=mean_abundance_nmol + sd_abundance_nmol), width=0.2) + 
    facet_wrap(~ metabolite, scales = "free_y", ncol = length(unique(temp2$metabolite))) +
    scale_colour_manual(values = c("wt 313" = wt, "rho0 313" = rho, "rho0 T911A" = brewer.pal(n = 8, name = palette)[c(1)])) +
    scale_fill_manual(values = c("wt 313" = wt, "rho0 313" = rho, "rho0 T911A" = brewer.pal(n = 8, name = palette)[c(1)])) +
    scale_y_continuous(expand = c(0,0))+ 
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="vertical",
          legend.background = element_blank(),
          strip.background = element_blank(),
          legend.title=element_blank())+
    ylab("metabolite abundance (nmol)")+
    xlab("")
}

## Figure 5----

### Fig. 5A-B
{
  # data
  fcsummary <- read.delim("Data/20210816_IronData_FoldChange_Summarized.tsv")
  fcsummary$Compound <- paste0(str_to_title(tolower(fcsummary$Compound)), "p")
  pdata <- read.delim("Data/20210816_IronData_Raw_Data.tsv")
  
  #PCA analysis
  odata <- pdata %>%
    dplyr::filter(Experiment == "Fe") %>%
    dplyr::mutate(experiment = "proteomics") %>%
    dplyr::select(R.FileName, PG.Genes, PG.Quantity, experiment, ID) %>%
    dplyr::rename(Sample = R.FileName, Compound = PG.Genes, value = PG.Quantity)
  
  odata <- odata %>%
    dplyr::filter(!Compound %in% "RSM26")
  
  odata$Compound <- paste0(str_to_title(tolower(odata$Compound)), "p")
  
  matrix <- cast(odata, Sample ~ Compound, value="value")
  
  matrix <- cbind(t(as.data.frame(strsplit(as.character(matrix$Sample), "_"))), matrix)
  matrix <- cbind(Strain=paste(matrix[,3], matrix[,4], matrix[,5], sep="_"), matrix[,c(-1,-2,-3, -4, -5)])
  names(matrix)[2] <- c("Replicate")
  
  for(i in 4:ncol(matrix)){
    if(T %in% is.na(matrix[,i])){
      t <- aggregate(matrix[,i], by=list(matrix$Strain), FUN=mean, na.rm=T)
      for(j in unique(matrix$Strain)){
        matrix[which(matrix$Strain == j),i][is.na(matrix[which(matrix$Strain == j),i])] <- t$x[which(t$Group.1 == j)]
      }
    }
  }
  matrix <- na.omit(matrix)
  
  matrix.pca <- prcomp(matrix[,-1:-3], scale. = T)
  loadings <- data.frame(matrix.pca$rotation)
  
  GO <- data.frame()
  for(i in list.files(path = "Data/", pattern = "_annotations.txt", full.names = T)){
    t <- read.csv(i, sep="\t", skip=8, stringsAsFactors = T)
    GO <- rbind(GO, t)
  }
  
  GO <- GO %>%
    dplyr::distinct(Gene, Gene.Ontology.Term) %>%
    dplyr::rename(Compound = Gene)
  levels(GO$Gene.Ontology.Term)
  
  #PC1
  loadings <- loadings[order(-abs(loadings$PC1)),]
  goi <- row.names(loadings)[1:50]
  
  
  subdata <- fcsummary %>%
    dplyr::filter(Compound %in% goi) %>%
    dplyr::mutate(sgdName = factor(Compound)) %>%
    dplyr::mutate(StrainID = factor(Strain, levels = c("WT Control + Iron", "Rho0 Control + Iron", "Rho0 ATP3-T911A + Iron", "WT Control - Iron", "Rho0 Control - Iron", "Rho0 ATP3-T911A - Iron")))
  
  ISC <- read.delim("Data/ironsulfur_cluster_binding_annotations.txt", sep="\t", skip=8, stringsAsFactors = T)
  for(i in unique(GO$Gene.Ontology.Term)){
    print(i)
    print(unique(subdata$Compound[subdata$Compound %in% paste0(str_to_title(tolower(GO$Compound[GO$Gene.Ontology.Term == i])), "p")]))
  }
  
 
  arrdata <- subdata %>%
    dplyr::group_by(sgdName) %>%
    dplyr::summarise(meanval = mean(log2(mean))) %>%
    dplyr::arrange(-meanval)
  
  subdata$sgdName <- factor(subdata$sgdName, levels = unique(arrdata$sgdName))
  
  Figures[["5A"]] <- ggplot(data = subdata, aes(x = sgdName, y = StrainID))+
    geom_raster(aes(fill = log2(mean)))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))

  #PC2
  loadings <- loadings[order(-abs(loadings$PC2)),]
  goi <- row.names(loadings)[1:50]
  
  
  subdata <- fcsummary %>%
    dplyr::filter(Compound %in% goi) %>%
    dplyr::mutate(sgdName = factor(Compound)) %>%
    dplyr::mutate(StrainID = factor(Strain, levels = c("WT Control + Iron", "Rho0 Control + Iron", "Rho0 ATP3-T911A + Iron", "WT Control - Iron", "Rho0 Control - Iron", "Rho0 ATP3-T911A - Iron")))
  
  
  arrdata <- subdata %>%
    dplyr::group_by(sgdName) %>%
    dplyr::summarise(meanval = mean(log2(mean))) %>%
    dplyr::arrange(-meanval)
  
  subdata$sgdName <- factor(subdata$sgdName, levels = unique(arrdata$sgdName))
  
  Figures[["5B"]] <- ggplot(data = subdata, aes(x = sgdName, y = StrainID))+
    geom_raster(aes(fill = log2(mean)))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))
}

### Fig. 5C
{
  data <- read.csv("Data/151019_ICP_Analysis_MeanData.csv", stringsAsFactors = T)
  data$Strain <- factor(data$Strain, levels=levels(data$Strain)[c(3,1,2)])
  names(data) <- sub("Â", "", names(data))
  idata <- read.csv("Data/151019_ICP_Analysis_SumData.csv", stringsAsFactors = T)
  idata$Strain <- factor(idata$Strain, levels=levels(idata$Strain)[c(3,1,2)])
  names(idata) <- sub("Â", "", names(idata))
  
  Figures[["5C"]] <- ggplot(data[which(data$Analyte == "Fe"),], aes(x = Strain, y = mean.nmol.µg_protein))+
    geom_bar(stat="identity", position=dodge, aes(fill=Strain), colour="black")+
    geom_errorbar(aes(ymin=mean.nmol.µg_protein-sd.nmol.µg_protein, ymax=mean.nmol.µg_protein+sd.nmol.µg_protein), width=0.2, position=dodge)+
    geom_text(aes(label=signif(pvalue, 2), y=mean.nmol.µg_protein+sd.nmol.µg_protein+0.01), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = idata[which(idata$Analyte == "Fe"),], aes(y = nmol.µg_protein), colour = "grey", width = 0.2)+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,0.23))+ 
    scale_x_discrete(labels=c(lb10, lb11, lb13))+
    xlab("")+
    ylab("cellular Iron [nmol/µg protein]")
}

### Fig. 5D
{
  data2 <- read.delim("Data/210301_Iron_Depletion_Processed.tsv")
  sdata <- data2 %>%
    dplyr::group_by(condition, concentration) %>%
    dplyr::summarise(mean = mean(mu.spline.norm), sd = sd(mu.spline.norm))
  
  Figures[["5D"]] <- ggplot(data2, aes(x = concentration, y = mu.spline.norm)) +
    geom_point(aes(colour = condition)) +
    geom_line(data = sdata, aes(y = mean, colour = condition))+
    geom_ribbon(data = sdata, aes(ymin = mean-sd, ymax = mean+sd, y = mean, fill = condition), alpha = 0.2 )+
    geom_hline(yintercept = 1, linetype = 2, colour = "grey")+
    scale_colour_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]),labels=c(lb10b,lb11b,lb13b))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]),labels=c(lb10b,lb11b,lb13b))+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(1,1), 
          legend.direction="vertical",
          legend.background = element_blank(),
          legend.title=element_blank())+
    #coord_cartesian(ylim=c(0,0.210))+
    scale_y_continuous(expand = c(0,0), limits = c(0.2, 1.05))+ 
    xlab("DIP concentration [mM]")+
    ylab(bquote('normalized growth rate'))
}

### Fig. 5E
{
  data <- read.csv("Data/160211_ISC_labelling_summary.tsv", sep="\t", stringsAsFactors = T) %>%
    dplyr::filter(sample %in% c("pUG35_Aco1 wt",
                                "pUG35_Aco1 rho0",
                                "pUG35_Aco1 rho0 ATP3-6",
                                "pUG35")) %>%
    dplyr::filter(fraction == "eluate") %>%
    dplyr::mutate(sample = factor(sample, levels = c("pUG35_Aco1 wt",
                                                     "pUG35_Aco1 rho0",
                                                     "pUG35_Aco1 rho0 ATP3-6",
                                                     "pUG35"))) %>%
    dplyr::mutate(signif_aco = ifelse(pvalue_aco < 0.01, "*", NA))
  idata <- read.csv("Data/160211_ISC_labelling_data.tsv", sep="\t", stringsAsFactors = T) %>%
    dplyr::filter(sample %in% c("pUG35_Aco1 wt",
                                "pUG35_Aco1 rho0",
                                "pUG35_Aco1 rho0 ATP3-6",
                                "pUG35")) %>%
    dplyr::filter(fraction == "eluate") %>%
    dplyr::mutate(sample = factor(sample, levels = c("pUG35_Aco1 wt",
                                                     "pUG35_Aco1 rho0",
                                                     "pUG35_Aco1 rho0 ATP3-6",
                                                     "pUG35")))

  
  Figures[["5E"]] <- ggplot(data, aes(x = sample, y = mean_cor))+
    geom_bar(stat="identity", position=dodge, aes(fill=sample), colour="black")+
    geom_errorbar(aes(ymin=mean_cor-sd_cor, ymax=mean_cor+sd_cor), width=0.2, position=dodge)+
    geom_text(aes(label=signif(pvalue_aco, 2), y=mean_cor+sd_cor+50), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = idata, aes(y = cor), colour = "grey", width = 0.2)+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)], "lightgrey"))+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,8000))+ 
    scale_x_discrete(labels=c(lb10, lb11, lb13, "GFP"))+
    xlab("")+
    ylab("Aco1 iron incorporation [cpm/mg protein]")
}

### Fig. 5F
{
  edata <- read.csv("Data/Enzyme_assays_p416Aco1_summary_v01_2021.csv", stringsAsFactors = T)
  edata$Strain <- factor(edata$Strain, levels=levels(edata$Strain)[c(3,2,1)])
  edata$Genotype <- factor(edata$Genotype, levels=levels(edata$Genotype)[c(2,1,3)])
  edata <- edata %>%
    dplyr::filter(!Strain == "aco1")
    
  
  sdata <- edata %>%
    dplyr::group_by(Strain, Genotype) %>%
    dplyr::summarise(unit_sd = sd(unit_mean, na.rm=T), unit_mean = mean(unit_mean), pvalue = mean(pvalue)) %>%
    dplyr::filter(!Strain == "aco1") %>%
    dplyr::mutate(signif = ifelse(pvalue < 0.01, "*", NA))
  
  Figures[["5F"]] <- ggplot(sdata, aes(x = Strain, y = unit_mean, group=Genotype, fill = Genotype))+
    geom_bar(stat="identity", position=dodge,aes(fill=Genotype, group=Genotype), colour="black")+
    geom_errorbar(aes(ymin=unit_mean-unit_sd, ymax=unit_mean+unit_sd), width=0.2, position=dodge)+
    geom_text(aes(label=signif(pvalue, 2), y=unit_mean+unit_sd+0.5), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_point(data = edata, colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    scale_fill_grey(breaks=levels(sdata$Genotype), labels=c("p416GPD","p416GPD_Aco1","p416GPD_mtAco1"))+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_x_discrete(breaks=levels(sdata$Strain), labels=c("wild type", lb5, lb34))+
    xlab("")+
    ylab("aconitase activity [mU/mg protein]")
}

### Fig. 5G
{
  edata <- read.csv("Data/Enzyme_assays_Iron_summary_v01.csv", stringsAsFactors = T)
  edata$Deferoxamine <- factor(edata$Deferoxamine, levels=levels(edata$Deferoxamine)[c(5,1,3,2,4)])
  edata$Iron <- factor(edata$Iron)
  
  Figures[["5G"]] <- ggplot(edata[which(edata$Deferoxamine == "w/o"),], aes(x = Iron, y = unit_mean, group=Deferoxamine))+
    geom_line(colour="black")+
    geom_ribbon(aes(ymin=unit_mean-unit_sd, ymax=unit_mean+unit_sd), alpha=0.2, fill="grey")+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="vertical",
          legend.title=element_blank())+
    geom_text(aes(label=signif(pvalue, 2), y=unit_mean+unit_sd), position=dodge, angle=90, hjust=0, vjust=0.8, size=10)+
    ylab("aconitase activity [mU/mg protein]")+
    xlab("Iron [II] chloride [ÂµM]")
}

### Fig. 5H
{
  edata <- read.csv("Data/Enzyme_assays_Iron_summary_v03.csv", stringsAsFactors = T)
  edata$Deferoxamine <- factor(edata$Deferoxamine, levels=levels(edata$Deferoxamine)[c(5,1,3,2,4)])
  edata$Iron <- factor(edata$Iron)
  t <- levels(edata$ID)[c(17,21,22)]
  edata$ID <- factor(edata$ID, levels=levels(edata$ID)[order(edata$Deferoxamine)])
  
  
  idata <- read.csv("Data/Enzyme_assays_Iron_data.csv", stringsAsFactors = T) %>%
    dplyr::filter(ID %in% t)
  
  Figures[["5H"]] <- ggplot(edata[which(edata$ID %in% t),], aes(x = ID, y = unit_mean))+
    geom_bar(stat="identity", position=dodge,aes(fill=Deferoxamine, group=Deferoxamine), colour="black")+
    geom_errorbar(aes(ymin=unit_mean-unit_sd, ymax=unit_mean+unit_sd), width=0.2, position=dodge)+
    geom_text(aes(label=signif(pvalue, 2), y=unit_mean+unit_sd+0.5), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = idata, aes(y = activity), colour = "grey", width = 0.2)+
    scale_fill_grey(breaks=levels(edata$Deferoxamine)[c(1,5)], labels=c("control","20mM DFO-B"))+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(1,1), 
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,16))+ 
    scale_x_discrete(labels=c("control", lb35, lb35))+
    xlab("")+
    ylab("aconitase activity [mU/mg protein]")
}

### Fig. 5I
{
  data <- read.csv("Data/growth_rates_aco.csv", stringsAsFactors = T)
  idata <- read.csv("Data/growth_rates_aco_raw.csv", stringsAsFactors = T)
  data$Genotype <- factor(data$Genotype, levels=levels(data$Genotype)[c(3,2,1)])
  data$ID <- factor(data$ID, levels=levels(data$ID)[c(5,6,3,4,1,2)])
  idata$Genotype <- factor(idata$Genotype, levels=levels(idata$Genotype)[c(3,2,1)])
  idata$ID <- factor(idata$ID, levels=levels(idata$ID)[c(5,6,3,4,1,2)])
  
  Figures[["5I"]] <- ggplot(data, aes(x = ID, y = mean))+
    geom_bar(stat="identity", position=dodge,aes(fill=ID), colour="black")+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, position=dodge)+
    geom_text(aes(label=signif(pvalue, 2), y=mean+sd+0.005), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = idata, aes(y = mu.model), colour = "grey", width = 0.2)+
    geom_jitter(data = idata[idata$ID == "BY4741 wt pRS313_ATP3T911A" & idata$mu.model < 0.16,], aes(y = mu.model), colour = "black", width = 0.2)+
    scale_fill_manual(values=c(wt, "grey", rho, brewer.pal(n = 8, name = palette)[c(1)], aco, "darkorange2"))+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,0.22))+ 
    scale_x_discrete(labels=c(lb10, lb10c, lb11, lb13, lb36, lb37))+
    xlab("")+
    ylab(bquote('growth rate [OD min'^-1~']'))
}

### Fig. 5J
{
  data <- data.table::fread("Data/20210816_AcoData_FoldChange_Data.tsv")
  matrix <- cast(data, Sample ~ sgdName, value="foldchange")
  matrix$Sample <- factor(matrix$Sample)
  levels(matrix$Sample) <- sub("_1$", " 1", levels(matrix$Sample))
  levels(matrix$Sample) <- sub("_2$", " 2", levels(matrix$Sample))
  levels(matrix$Sample) <- sub("_3$", " 3", levels(matrix$Sample))
  matrix <- cbind(t(as.data.frame(strsplit(as.character(matrix$Sample), " "))), matrix)
  names(matrix)[1:2] <- c("Strain", "Replicate")
  
  for(i in 4:ncol(matrix)){
    if(T %in% is.na(matrix[,i])){
      t <- aggregate(matrix[,i], by=list(matrix$Strain), FUN=mean, na.rm=T)
      for(j in unique(matrix$Strain)){
        matrix[which(matrix$Strain == j),i][is.na(matrix[which(matrix$Strain == j),i])] <- t$x[which(t$Group.1 == j)]
      }
    }
  }
  matrix <- na.omit(matrix)
  
  matrix.pca <- prcomp(matrix[,-1:-3], scale. = T)
  
  matrix.class <- factor(na.omit(matrix$Strain))
  matrix.class <- factor(matrix.class, levels=levels(matrix.class)[c(5,6,3,4,1,2)])
  
  
  Figures[["5J"]] <- ggbiplot(matrix.pca,
           choices = c(2,3), 
           groups=matrix.class,
           var.scale=0,
           ellipse = T,
           var.alpha=0,
           select = NA)+
    theme_classic() +
    theme(legend.justification=c(0,0), 
          legend.position="bottom", 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.direction="vertical",
          legend.background = element_blank(),
          legend.title=element_blank())+
    #scale_color_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]), labels=c(lb10b,lb11b,lb13b))+
    geom_hline(yintercept=0, alpha=0.5)+
    geom_vline(xintercept=0, alpha=0.5)+
    xlim(-2.5,1.5)+
    ylim(-1.5,1.5)
}

## Supplementary Figures----

### Suppl. Fig. 5
{
  cdata <- read.delim("Data/20130502_atp3 growth curves_combined_edited_tech triplicates.csv", sep=",")
  time <- as.numeric(unlist(cdata[1,]))
  cdata <- as.data.frame(t(cdata))
  names(cdata) <- c("Time", as.character(unlist(cdata[1,]))[-1])
  cdata <- cdata[c(-1:-3),]
  cdata <- reshape::melt(cdata, id = "Time")
  names(cdata)[2] <- "Strain"
  cdata$value <- as.numeric(as.character(cdata$value))
  cdata$StrainID <- sub(" c.+", "", as.character(cdata$Strain))
  cdata$Time <- as.numeric(as.character(cdata$Time))
  cdata <- cdata[!cdata$Strain == "WT pRS313_ATP3_mut c1",]
  
  sdata <- cdata %>% dplyr::group_by(StrainID, Time) %>% dplyr::summarise(mean = mean(value), sd = sd(value))
  sdata <- as.data.frame(sdata)
  
  sdata <- sdata[grepl("WT", sdata$StrainID),]
  sdata$StrainID <- factor(sdata$StrainID)
  sdata$StrainID <- factor(sdata$StrainID, levels = levels(sdata$StrainID)[c(1,2,5,3,4)])
  
  Figures[["S5"]] <- ggplot(sdata[grepl("WT", sdata$StrainID),], aes(x = Time, y = mean, group = StrainID))+
    geom_line(aes(colour = StrainID))+
    geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd, fill = StrainID), alpha = 0.2)+
    scale_fill_manual(breaks=levels(sdata$StrainID), labels=c(lb10,lb12,lb13,lb14,lb15), values=c("black", "red", brewer.pal(n = 8, name = palette)[1:4]))+
    scale_colour_manual(breaks=levels(sdata$StrainID), labels=c(lb10,lb12,lb13,lb14,lb15), values=c("black", "red", brewer.pal(n = 8, name = palette)[1:4]))+
    theme_classic()+
    #guides(fill=guide_legend(override.aes = list(colour=NULL)))+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="bottom", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_y_continuous(expand = c(0,0))+ 
    ylab(bquote('growth rate [OD min'^-1~']'))
}

### Suppl. Fig. 7A
{
  data <- read.csv("Data/mtPotential_YSBN1_130815_selected_summary.csv", stringsAsFactors = T)
  data <- cbind(data, Genotype="rho0 evolved", signif=NA)
  data$Genotype <- as.character(data$Genotype)
  data$Genotype[grepl("wt", data$genotype)] <- "wild type"
  data$Genotype[grepl("Rho0 Chr16", data$genotype)] <- "rho0"
  data$Genotype <- factor(data$Genotype)
  data$Genotype <- factor(data$Genotype, levels=levels(data$Genotype)[c(3,1,2)])
  data$signif[which(data$pValue < 0.01)] <- "*"
  levels(data$genotype) <- c("ATP3-9", "ATP3-7", "ATP3-8", "ATP3-6", "wild type ", "wild type")
  data$genotype <- factor(data$genotype, levels=levels(data$genotype)[c(6,5,4,2,3,1)])
  data$clone <- factor(data$clone, levels=levels(data$clone)[c(6,5,3,2,4,1)])
  
  idata <- read.csv("Data/mtPotential_YSBN1_130815_selected_data.csv", stringsAsFactors = T) %>%
    dplyr::inner_join(., data[,c("clone", "genotype")], by = "clone") %>%
    dplyr::mutate(Genotype = ifelse(grepl("wt", genotype), "wild type", NA)) %>%
    dplyr::mutate(Genotype = ifelse(grepl("Rho0 Chr16", genotype), "rho0", Genotype)) %>%
    dplyr::mutate(Genotype = ifelse(is.na(Genotype), "rho0 evolved", Genotype))
  
  idata$Genotype <- factor(idata$Genotype)
  idata$Genotype <- factor(idata$Genotype, levels=levels(idata$Genotype)[c(3,1,2)])
  levels(idata$genotype) <- c("ATP3-9", "ATP3-7", "ATP3-8", "ATP3-6", "wild type ", "wild type")
  idata$genotype <- factor(idata$genotype, levels=levels(idata$genotype)[c(6,5,4,2,3,1)])
  idata$clone <- factor(idata$clone, levels=levels(idata$clone)[c(6,5,3,2,4,1)])
  
  Figures[["S7A"]] <- ggplot(data, aes(x = clone, y = GPC_mean))+
    geom_bar(stat="identity", aes(fill=clone, group=clone), colour="black")+
    scale_fill_manual(values=c("black", "red", brewer.pal(n = 8, name = palette)[1:4]), breaks=levels(data$clone), labels=c(lb6, lb5, lb27, lb28, lb29, lb30))+
    geom_text(aes(label=signif(pValue, 2), y=GPC_mean+GPC_se+0.01), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_text(aes(label = paste("n =", n)), y = 0.4, angle = 90, hjust = 0, vjust = 0.5, fontface = "bold", size = 6)+
    geom_text(data = data[data$clone %in% c("wt", "init"),], aes(label = paste("n =", n)), y = 0.4, angle = 90, hjust = 0, vjust = 0.5, fontface = "bold", colour = "white", size = 6)+
    #geom_jitter(data = idata, aes(y = GPC))+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    coord_cartesian(ylim=c(0.4,0.8))+
    scale_x_discrete(breaks=levels(data$clone), labels=c(lb6, lb5, lb27, lb28, lb29, lb30))+
    ylab(bquote('Pearson Correlation Coefficient'))+
    geom_errorbar(aes(ymin=GPC_mean-GPC_se, ymax=GPC_mean+GPC_se), width=0.2)+
    xlab("")
}

### Suppl. Fig. 7B
{
  data_psi <- read.csv("Data/mtPotential_YSBN11_130805_summary.csv") %>%
    dplyr::filter(Reference == "wt pRS313") %>%
    mutate(genotype = sub("wt", "wild type", genotype)) %>%
    select(genotype, plasmid, mean, se) %>%
    dplyr::rename(mean_psi = mean, sd_psi = se)
  
  data_grwth <- read.csv("Data/Growcurves_YSBN11_130914_fit_values_summary.csv") %>%
    select(genotype, plasmid, mean, sd) %>%
    dplyr::rename(mean_growth = mean, sd_growth = sd) %>%
    inner_join(data_psi, by = c("genotype", "plasmid")) %>%
    filter(!plasmid == "pRS313_ATP3-67") %>%
    mutate(StrainID = factor(paste(genotype, plasmid), levels = c("wild type pRS313", 
                                                                  "rho0 pRS313", 
                                                                  "rho0 pRS313_ATP3",
                                                                  "rho0 pRS313_ATP3-6",
                                                                  "rho0 pRS313_ATP3-7") ))
  
  
  t <- summary(lm(data = data_grwth, mean_growth ~ mean_psi))
  t$r.squared # R2 = 0.93
  model <- as.data.frame(t$coefficients)
  
  Figures[["S7B"]] <- ggplot(data_grwth, aes(x = mean_psi, y = mean_growth))+
    geom_point(aes(colour = StrainID))+
    geom_errorbar(aes(ymin = mean_growth-sd_growth, ymax = mean_growth+sd_growth, colour = StrainID))+
    geom_errorbarh(aes(xmin = mean_psi-sd_psi, xmax = mean_psi+sd_psi, colour = StrainID))+
    geom_abline(intercept = model$Estimate[1], slope = model$Estimate[2], linetype = 2, size = 1)+
    scale_colour_manual(values = c(wt, rho, brewer.pal(n = 8, name = palette)[c(5,1,2)]))+
    labs(y = bquote('growth rate [OD min'^-1~']'), x = "Mitochondrial Membrane Potential (PCC)")+
    theme_classic()+
    theme(legend.position = c(1,0),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.justification = c(1,0))
  
}

### Suppl. Fig. 8
{
  data <- read.csv("Data/151108_fscore.csv", stringsAsFactors = T)
  data$ConditionID <- factor(data$ConditionID, levels=levels(data$ConditionID)[c(4,1,3,2)])
  data <- cbind(data, signif=NA)
  data$signif[which(data$ttest_wt < 0.01)] <- "*"
  
  idata <- read.csv("Data/151108_fscore_data.csv") %>%
    dplyr::mutate(ConditionID = factor(ConditionID, levels = levels(data$ConditionID)))
  
  sdata <- idata %>% group_by(ConditionID) %>% dplyr::summarise(sum = length(f))
  
  Figures[["S8"]] <- ggplot(data, aes(x = ConditionID, y = f_mean))+
    geom_bar(stat="identity", aes(fill=ConditionID, group=ConditionID), colour="black")+
    geom_errorbar(aes(ymin=f_mean-f_sem, ymax=f_mean+f_sem), width=0.2)+
    geom_text(aes(label=signif(ttest_rho0, 2), y=f_mean+f_sem+5), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    #geom_jitter(data = idata, aes(y = f), colour = "grey", width = 0.2)+
    geom_text(data = sdata, aes(label = paste("n =", sum)), y = 5, angle = 90, fontface = "bold", size = 6)+
    geom_text(data = sdata[sdata$ConditionID == "wt_313",], aes(label = paste("n =", sum)), y = 5, angle = 90, fontface = "bold", colour = "white", size = 6)+
    scale_fill_manual(values=c("black", "red", brewer.pal(n = 8, name = palette)[1:2]), breaks=levels(data$ConditionID), labels=c(lb10,lb11,lb13,lb14))+
    theme_classic()+
    #guides(fill=guide_legend(override.aes = list(colour=NULL)))+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,80))+ 
    scale_x_discrete(breaks=levels(data$ConditionID), labels=c(lb10,lb11,lb13,lb14))+
    ylab(bquote('fragmented mitochondria [%]'))+
    xlab("")
}

### Suppl. Fig. 9
{
  sdata <- read.csv("Data/time_course_viability_OD_EC_summary.tsv", sep="\t", stringsAsFactors = T) %>%
    dplyr::filter(name == "EC") %>%
    dplyr::filter(round(timestamp) == 3)
  
  idata <- read.csv("Data/time_course_viability_OD_EC.tsv", sep="\t", stringsAsFactors = T) %>%
    dplyr::filter(name == "EC") %>%
    dplyr::filter(round(timestamp_) == 3)
  
  Figures[["S9"]] <- ggplot(sdata, aes(x = Strain, y = value_mean))+
    geom_bar(aes(fill = Strain), colour = "black",stat = "identity", position = position_dodge(width = 0.9))+
    geom_errorbar(aes(ymin = value_mean-value_sd, ymax = value_mean+value_sd), width = 0.2)+
    geom_jitter(data = idata, aes(y = value), colour = "grey", width = 0.2)+
    scale_fill_manual(labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_y_continuous(expand = c(0,0), breaks=seq(0.1,1.0, by=0.1), limits = c(0,1), name = "energy charge")+
    scale_x_discrete(breaks=levels(sdata$Strain), labels=c(lb10,lb11,lb13))+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())
}

### Suppl. Fig. 10
{
  data <- read.csv("Data/time_course_viability_OD_EC_summary.tsv", sep="\t", stringsAsFactors = T)
  data$name <- factor(data$name, levels=levels(data$name)[c(4,5,7,6,3,1,2)])
  labels <- list('biomass [OD'[600]*']', 
                 'medium glucose [mM]',
                 'cell viability [cfu/OD]',
                 'energy charge',
                 'ATP [mM]',
                 'ADP [mM]',
                 'AMP [mM]')
  
  p <- list()
  p[[1]] <- ggplot(data[data$name == levels(data$name)[1] & !is.na(data$value_mean),], aes(x = time_sec, y = value_mean, group=Strain))+
    geom_line(aes(colour=Strain))+
    geom_ribbon(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd, fill=Strain), alpha=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_x_continuous(expand = c(0,0), limits=c(0,80))+
    ylab(bquote('biomass [OD'[600]*']'))+
    xlab("culture time [h]")+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_colour_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))
  
  p[[2]] <- ggplot(data[data$name == levels(data$name)[2] & !is.na(data$value_mean),], aes(x = time_sec, y = value_mean, group=Strain))+
    geom_line(aes(colour=Strain))+
    geom_ribbon(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd, fill=Strain), alpha=0.2)+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(1,1), 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_x_continuous(expand = c(0,0), limits=c(0,80))+
    ylab(bquote('medium glucose [mM]'))+
    xlab("culture time [h]")+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_colour_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))
  
  p[[3]] <- ggplot(data[data$name == levels(data$name)[3] & !is.na(data$value_mean),], aes(x = time_sec, y = value_mean, group=Strain))+
    geom_line(aes(colour=Strain))+
    geom_ribbon(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd, fill=Strain), alpha=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_x_continuous(expand = c(0,0), limits=c(0,80))+
    scale_y_log10()+
    annotation_logticks(sides = "l")+
    ylab(bquote('cell viability [cfu/OD]'))+
    xlab("culture time [h]")+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_colour_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))
  
  p[[4]] <- ggplot(data[data$name == levels(data$name)[4] & !is.na(data$value_mean),], aes(x = time_sec, y = value_mean, group=Strain))+
    geom_line(aes(colour=Strain))+
    geom_ribbon(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd, fill=Strain), alpha=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_x_continuous(expand = c(0,0), limits=c(10,60))+
    scale_y_continuous(expand = c(0,0), breaks=seq(0.1,0.9, by=0.1))+
    ylab(bquote('energy charge'))+
    xlab("culture time [h]")+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_colour_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))
  
  p[[5]] <- ggplot(data[data$name == levels(data$name)[5] & !is.na(data$value_mean),], aes(x = time_sec, y = value_mean, group=Strain))+
    geom_line(aes(colour=Strain))+
    geom_ribbon(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd, fill=Strain), alpha=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position="none", 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_x_continuous(expand = c(0,0), limits=c(0,80))+
    ylab(bquote('ATP [mM]'))+
    xlab("culture time [h]")+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_colour_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))
  
  p[[6]] <- ggplot(data[data$name == levels(data$name)[6] & !is.na(data$value_mean),], aes(x = time_sec, y = value_mean, group=Strain))+
    geom_line(aes(colour=Strain))+
    geom_ribbon(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd, fill=Strain), alpha=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position="none", 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_x_continuous(expand = c(0,0), limits=c(0,80))+
    ylab(bquote('ADP [mM]'))+
    xlab("culture time [h]")+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_colour_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))
  
  p[[7]] <- ggplot(data[data$name == levels(data$name)[7] & !is.na(data$value_mean),], aes(x = time_sec, y = value_mean, group=Strain))+
    geom_line(aes(colour=Strain))+
    geom_ribbon(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd, fill=Strain), alpha=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_x_continuous(expand = c(0,0), limits=c(0,80))+
    ylab(bquote('AMP [mM]'))+
    xlab("culture time [h]")+
    scale_fill_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_colour_manual(breaks=levels(data$Strain), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho,brewer.pal(n = 8, name = palette)[c(1)]))
  
}

### Suppl. Fig. 11
{
  data <- read.csv("Data/fermentations_data.csv", stringsAsFactors = T)
  max <- read.csv("Data/fermentations_diauxic.csv", stringsAsFactors = T)
  gases <- read.csv("Data/fermentations_data_gases.csv", stringsAsFactors = T)
  data <- rbind(data, gases)
  data$data <- factor(data$data)
  data$compound <- factor(data$compound)
  levels(data$compound)[4] <- "biomass"
  data$variable <- factor(data$variable, levels=levels(data$variable)[c(1,2,4,3)])
  max$variable <- factor(max$variable, levels=levels(max$variable)[c(1,2,4,3)])
  
  # faceting according to compounds
  {
    #pdf("data/fermentations_data.pdf", width = 9, height = 6)
    #for(i in  unique(data$compound[grepl("O2", data$compound)])){
    i <- "nCO2" 
    Figures[["S11A"]] <- ggplot(data[which(data$compound == i),], aes(x = timepoint, y = value, group=variable))+
      geom_line(aes(colour=variable), size=1.5, linetype=1)+
      geom_point(aes(colour=variable), size=3)+
      geom_vline(data=max, aes(xintercept=timepoint, colour = variable), linetype=2, size=0.5)+
      theme_classic()+
      theme(legend.justification=c(1,1), 
            legend.position="none", 
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            legend.direction="vertical",
            legend.title=element_blank(),
            #legend.key=element_rect(colour = "black"),
            legend.background=element_blank(),
            strip.background = element_blank())+
      scale_x_continuous("time [hrs]")+
      scale_y_continuous(paste0(i, " (", unique(data$unit[data$compound == i]), ")"), limits = c(0,max(data$value[data$compound == i])))+
      scale_colour_manual(labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
      annotate("text",x = 0, y = max(data$value[data$compound == i]), label=i, parse=F, hjust=0)
    
    i <- "nO2" 
    Figures[["S11B"]] <- ggplot(data[which(data$compound == i),], aes(x = timepoint, y = value, group=variable))+
      geom_line(aes(colour=variable), size=1.5, linetype=1)+
      geom_point(aes(colour=variable), size=3)+
      geom_vline(data=max, aes(xintercept=timepoint, colour = variable), linetype=2, size=0.5)+
      theme_classic()+
      theme(legend.justification=c(1,1), 
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            legend.position="none", 
            legend.direction="vertical",
            legend.title=element_blank(),
            #legend.key=element_rect(colour = "black"),
            legend.background=element_blank(),
            strip.background = element_blank())+
      scale_x_continuous("time [hrs]")+
      scale_y_continuous(paste0(i, " (", unique(data$unit[data$compound == i]), ")"), limits = c(0,max(data$value[data$compound == i])))+
      scale_colour_manual(labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
      annotate("text",x = 0, y = max(data$value[data$compound == i]), label=i, parse=F, hjust=0)
    
    i <- "CO2-rate" 
    Figures[["S11C"]] <- ggplot(data[which(data$compound == i),], aes(x = timepoint, y = value, group=variable))+
      geom_line(aes(colour=variable), size=1.5, linetype=1)+
      geom_point(aes(colour=variable), size=3)+
      geom_vline(data=max, aes(xintercept=timepoint, colour = variable), linetype=2, size=0.5)+
      theme_classic()+
      theme(legend.justification=c(1,1), 
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            legend.position="none", 
            legend.direction="vertical",
            legend.title=element_blank(),
            #legend.key=element_rect(colour = "black"),
            legend.background=element_blank(),
            strip.background = element_blank())+
      scale_x_continuous("time [hrs]")+
      scale_y_continuous(paste0(i, " (", unique(data$unit[data$compound == i]), ")"), limits = c(0,max(data$value[data$compound == i])))+
      scale_colour_manual(labels=c(lb10,lb11,lb13,lb14), values=c("black", "red",brewer.pal(n = 8, name = palette)[c(1,2)]))+
      annotate("text",x = 0, y = max(data$value[data$compound == i]), label=i, parse=F, hjust=0)
    #   print(p)
    # }
    #dev.off()
  }
}

### Suppl. Fig. 12A
{
  fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
  signif <- read.delim("Data/20210816_CombinedData_Signif.tsv")
  poi <- c("GDH1", "GLT1", "MEU1", "PUT2", "GLN1")
  sdata <- fcsummary[which(fcsummary$Compound %in% poi),]
  ssignif <- signif[signif$Compound %in% poi,]
  ssignif$p.value <- signif(ssignif$p.value, 2)
  ssignif$plabel <- NA
  ssignif$plabel[ssignif$p.value < 1E-5] <- "p < 1E-5"
  ssignif$plabel[!ssignif$p.value < 1E-5] <- paste0("p = ", ssignif$p.value[!ssignif$p.value < 1E-5])
  poisdata <- fcdata[fcdata$Compound %in% poi,]
  sdata$sgdName <- paste0(str_to_title(tolower(sdata$Compound)), "p")
  
  Figures[["S12A"]] <- ggplot(sdata, aes(x = Strain, y = log2(mean), group=Compound, fill = Strain))+
    geom_bar(stat="identity", aes(fill=Strain), colour="black", position=dodge)+
    geom_errorbar(aes(ymin=log2(mean-sd), ymax=log2(mean+sd)), width=0.2, position=dodge)+
    geom_point(data = poisdata, aes(y = log2(foldchange)), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    geom_text(data = ssignif, aes(label = signif(p.value, 2), y = 2.5), angle = 90, hjust = 0.5, position = dodge)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          legend.title=element_blank())+
    ylab("log2 fold change")+
    ylim(min(log2(sdata$mean-sdata$sd))-0.8, 3)+
    geom_hline(yintercept=0, colour="grey")+
    xlab("")+
    geom_text(aes(label=sgdName, y=min(log2(sdata$mean-sdata$sd))-0.05), position=dodge, angle=90, hjust=1, vjust=0.5, size=6)+
    #scale_x_discrete(breaks=levels(Pdata$Sample), labels=c(lb10, lb11, lb13))+
    scale_fill_manual(breaks=unique(sdata$Strain), labels=c(lb10,lb11,lb13), values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))
}

### Suppl. Fig. 12B
{
  fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
  signif <- read.delim("Data/20210816_CombinedData_Signif.tsv")
  poi <- c("ILV2", "ILV5", "ILV3", "LEU9", "LEU4", "LEU1", "LEU2", "BAT1", "BAT2")
  sdata <- fcsummary[which(fcsummary$Compound %in% poi),]
  ssignif <- signif[signif$Compound %in% poi,]
  ssignif$p.value <- signif(ssignif$p.value, 2)
  ssignif$plabel <- NA
  ssignif$plabel[ssignif$p.value < 1E-5] <- "p < 1E-5"
  ssignif$plabel[!ssignif$p.value < 1E-5] <- paste0("p = ", ssignif$p.value[!ssignif$p.value < 1E-5])
  
  
  poisdata <- fcdata[fcdata$Compound %in% poi,]
  sdata$sgdName <- paste0(str_to_title(tolower(sdata$Compound)), "p")
  
  Figures[["S12B"]] <- ggplot(sdata, aes(x = Strain, y = log2(mean), group=Compound, fill = Strain))+
    geom_bar(stat="identity", aes(fill=Strain), colour="black", position=dodge)+
    geom_errorbar(aes(ymin=log2(mean-sd), ymax=log2(mean+sd)), width=0.2, position=dodge)+
    geom_point(data = poisdata, aes(y = log2(foldchange)), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    geom_text(data = ssignif, aes(label = signif(p.value, 2), y = 2.5), angle = 90, hjust = 0.5, position = dodge)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position="none", 
          legend.direction="vertical",
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          legend.title=element_blank())+
    ylab("log2 fold change")+
    ylim(min(log2(sdata$mean-sdata$sd))-0.8, 3)+
    geom_hline(yintercept=0, colour="grey")+
    xlab("")+
    geom_text(aes(label=sgdName, y=min(log2(sdata$mean-sdata$sd))-0.05), position=dodge, angle=90, hjust=1, vjust=0.5, size=6)+
    #scale_x_discrete(breaks=levels(Pdata$Sample), labels=c(lb10, lb11, lb13))+
    scale_fill_manual(breaks=unique(sdata$Strain), labels=c(lb10,lb11,lb13), values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))
}

### Suppl. Fig. 12C
{
  fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
  signif <- read.delim("Data/20210816_CombinedData_Signif.tsv")
  poi <- c("PUT1", "PUT2", "ARG5,6", "ARG7", "ARG8", "ORT1", "CPA1", "CPA2", "ARG1", "ARG3", "ARG4", "CAN1", "CAR1", "PRO3")
  sdata <- fcsummary[which(fcsummary$Compound %in% poi),]
  ssignif <- signif[signif$Compound %in% poi,]
  ssignif$p.value <- signif(ssignif$p.value, 2)
  ssignif$plabel <- NA
  ssignif$plabel[ssignif$p.value < 1E-5] <- "p < 1E-5"
  ssignif$plabel[!ssignif$p.value < 1E-5] <- paste0("p = ", ssignif$p.value[!ssignif$p.value < 1E-5])
  poisdata <- fcdata[fcdata$Compound %in% poi,]
  sdata$sgdName <- paste0(str_to_title(tolower(sdata$Compound)), "p")
  
  Figures[["S12C"]] <- ggplot(sdata, aes(x = Strain, y = log2(mean), group=Compound, fill = Strain))+
    geom_bar(stat="identity", aes(fill=Strain), colour="black", position=dodge)+
    geom_errorbar(aes(ymin=log2(mean-sd), ymax=log2(mean+sd)), width=0.2, position=dodge)+
    geom_point(data = poisdata, aes(y = log2(foldchange)), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    geom_text(data = ssignif, aes(label = signif(p.value, 2), y = 2.5), angle = 90, hjust = 0.5, position = dodge)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position="none", 
          legend.direction="vertical",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          legend.title=element_blank())+
    ylab("log2 fold change")+
    ylim(min(log2(sdata$mean-sdata$sd))-0.8, 3)+
    geom_hline(yintercept=0, colour="grey")+
    xlab("")+
    geom_text(aes(label=sgdName, y=min(log2(sdata$mean-sdata$sd))-0.05), position=dodge, angle=90, hjust=1, vjust=0.5, size=6)+
    #scale_x_discrete(breaks=levels(Pdata$Sample), labels=c(lb10, lb11, lb13))+
    scale_fill_manual(breaks=unique(sdata$Strain), labels=c(lb10,lb11,lb13), values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))
}

### Suppl. Fig. 13A
{
  data <- read.csv("Data/150917_growth_outliersrm.csv", stringsAsFactors = T)
  data_raw <- read.csv("Data/150917_growth_outliersrm_raw_defect.csv", stringsAsFactors = T)
  data$StrainID <- factor(data$StrainID, levels=levels(data$StrainID)[c(3,1,2)])
  data$Mix <- factor(data$Mix, levels=levels(data$Mix)[c(1,6,3,11,5,7,12,8,10,4,9,2)])  
  data <- cbind(data, signif = NA)
  data$signif[which(data$pvalue < 0.01)] <- "*"
  
  data <- data %>%
    dplyr::filter(Mix %in% c("-AA", "+L"))
  
  data_raw <- data_raw %>%
    dplyr::filter(Mix %in% c("-AA", "+L"))
  
  Figures[["S13A"]] <- ggplot(data, aes(x = StrainID, y = mu.spline_mean, group=Mix, fill = StrainID))+
    geom_bar(aes(fill=StrainID, alpha=Mix), stat="identity", position=dodge, colour="black")+
    geom_hline(yintercept = data$mu.spline_mean[data$Mix == "-AA" & data$StrainID == "YSBN11 wt 313"], linetype = 2)+
    geom_text(aes(label=signif(pvalue, 2), y=mu.spline_mean+mu.spline_sd+0.005), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_point(data = data_raw, aes(y = mu.spline), position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2), colour = "grey")+
    geom_errorbar(data = data, aes(ymin=mu.spline_mean-mu.spline_sd, ymax=mu.spline_mean+mu.spline_sd), width=0.2, position=dodge)+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]),labels=c(lb10b,lb11b,lb13b), guide=F)+
    theme_classic()+
    scale_alpha_discrete(range = c(0.5,1), guide=F)+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="horizontal",
          legend.background = element_blank(),
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,0.18))+ 
    xlab("amino acid supplementation")+
    ylab(bquote('growth rate [OD min'^-1~']'))
}

### Suppl. Fig. 13B
{
  data <- read.csv("Data/150917_growth_outliersrm.csv", stringsAsFactors = T)
  data_raw <- read.csv("Data/150917_growth_outliersrm_raw_defect.csv", stringsAsFactors = T)
  data$StrainID <- factor(data$StrainID, levels=levels(data$StrainID)[c(3,1,2)])
  data$Mix <- factor(data$Mix, levels=levels(data$Mix)[c(1,6,3,11,5,7,12,8,10,4,9,2)])  
  data <- cbind(data, signif = NA)
  data$signif[which(data$pvalue < 0.01)] <- "*"
  
  
  data <- data %>%
    dplyr::filter(Mix %in% c("-AA", "+R"))
  
  data_raw <- data_raw %>%
    dplyr::filter(Mix %in% c("-AA", "+R"))
  
  Figures[["S13B"]] <- ggplot(data, aes(x = StrainID, y = mu.spline_mean, group=Mix, fill = StrainID))+
    geom_bar(aes(fill=StrainID, alpha=Mix), stat="identity", position=dodge, colour="black")+
    geom_hline(yintercept = data$mu.spline_mean[data$Mix == "-AA" & data$StrainID == "YSBN11 wt 313"], linetype = 2)+
    geom_text(aes(label=signif(pvalue, 2), y=mu.spline_mean+mu.spline_sd+0.005), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_point(data = data_raw, aes(y = mu.spline), position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2), colour = "grey")+
    geom_errorbar(data = data, aes(ymin=mu.spline_mean-mu.spline_sd, ymax=mu.spline_mean+mu.spline_sd), width=0.2, position=dodge)+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]),labels=c(lb10b,lb11b,lb13b), guide=F)+
    theme_classic()+
    scale_alpha_discrete(range = c(0.5,1), guide=F)+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="horizontal",
          legend.background = element_blank(),
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,0.18))+ 
    xlab("amino acid supplementation")+
    ylab(bquote('growth rate [OD min'^-1~']'))
}

### Suppl. Fig. 13C
{
  data <- read.csv("Data/150917_growth_outliersrm.csv", stringsAsFactors = T)
  data_raw <- read.csv("Data/150917_growth_outliersrm_raw_defect.csv", stringsAsFactors = T)
  data$StrainID <- factor(data$StrainID, levels=levels(data$StrainID)[c(3,1,2)])
  data$Mix <- factor(data$Mix, levels=levels(data$Mix)[c(1,6,3,11,5,7,12,8,10,4,9,2)])  
  data <- cbind(data, signif = NA)
  data$signif[which(data$pvalue < 0.01)] <- "*"
  
  data <- data %>%
    dplyr::filter(Mix %in% c("-AA", "+Q", "+E"))
  
  data_raw <- data_raw %>%
    dplyr::filter(Mix %in% c("-AA", "+Q", "+E"))
  
  
  Figures[["S13C"]] <- ggplot(data, aes(x = StrainID, y = mu.spline_mean, group=Mix, fill = StrainID))+
    geom_bar(aes(fill=StrainID, alpha=Mix), stat="identity", position=dodge, colour="black")+
    geom_hline(yintercept = data$mu.spline_mean[data$Mix == "-AA" & data$StrainID == "YSBN11 wt 313"], linetype = 2)+
    geom_text(aes(label=signif(pvalue, 2), y=mu.spline_mean+mu.spline_sd+0.005), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_point(data = data_raw, aes(y = mu.spline), position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2), colour = "grey")+
    geom_errorbar(data = data, aes(ymin=mu.spline_mean-mu.spline_sd, ymax=mu.spline_mean+mu.spline_sd), width=0.2, position=dodge)+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]),labels=c(lb10b,lb11b,lb13b), guide=F)+
    theme_classic()+
    scale_alpha_discrete(range = c(0.5,1), guide=F)+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="horizontal",
          legend.background = element_blank(),
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,0.18))+ 
    xlab("amino acid supplementation")+
    ylab(bquote('growth rate [OD min'^-1~']'))
}

### Suppl. Fig. 14
{
  data <- data.table::fread("Data/20210817_mtPotential_FACS.xls", stringsAsFactors = T) %>%
    dplyr::filter(StrainID %in% c("Rho+ pRS313", "Rho0 pRS313", "Rho0 ATP3 T911A", "Rho0 ATP3 G919C")) %>%
    dplyr::filter(Replicate == 2) %>%
    dplyr::mutate(StrainID = factor(StrainID)) %>%
    dplyr::mutate(StrainID = factor(StrainID, levels = c("Rho+ pRS313", "Rho0 pRS313", "Rho0 ATP3 T911A", "Rho0 ATP3 G919C"))) %>%
    dplyr::mutate(Sample = paste(StrainID, Medium))

  Figures[["S14"]] <- ggplot(data, aes(x = value, y = ..scaled..))+
    geom_density(aes(colour = Sample, fill = Sample), alpha = 0.2)+
    scale_x_log10(limits = c(100,1000), name = "DiOC6 Fluorescence")+
    scale_y_continuous(expand = c(0,0), name = "Relative Number of Cells") +
    facet_wrap(~ StrainID)+
    theme_classic()+
    theme(legend.position="bottom", 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.title=element_blank(),
          axis.line = element_line(),
          strip.background = element_blank(),
          legend.background=element_blank())
}

### Suppl. Fig. 15B
{
  data <- read.csv("Data/140407_replicates_summary.csv", stringsAsFactors = T)
  data$condition <- factor(data$condition, levels=levels(data$condition)[c(1,5,6,8,7,2,3,4)])
  data <- data %>%
    dplyr::filter(condition %in% levels(data$condition)[c(1:3,5:6)] & compound == levels(data$compound)[1])
  
  idata <- read.csv("Data/140407_replicates_data.csv", stringsAsFactors = T) %>%
    dplyr::filter(strainID %in% levels(data$condition)[c(1:3,5:6)] & variable == levels(data$compound)[1]) %>%
    dplyr::select(strainID, variable, value) %>%
    dplyr::rename(condition = strainID, compound = variable)
  
  
  
  Figures[["S15B"]] <- ggplot(data, aes(x = condition, y = mean, group=condition))+
    geom_bar(aes(fill=condition), stat="identity", colour = "black")+
    geom_text(aes(label=signif(pvalue, 2), y=mean+sd+0.2), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = idata, aes(y = value), width = 0.2, colour = "grey")+
    geom_errorbar(data = data, aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,40))+
    ylab(bquote('GSH [mM]'))+
    xlab("")+
    scale_x_discrete(labels=c(lb10,lb11, lb12, lb14, lb38))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(5)], brewer.pal(n = 8, name = palette)[c(2)], wt))
}

### Suppl. Fig. 15C
{
  data <- read.csv("Data/140407_replicates_summary.csv", stringsAsFactors = T)
  data$condition <- factor(data$condition, levels=levels(data$condition)[c(1,5,6,8,7,2,3,4)])
  data <- data %>%
    dplyr::filter(condition %in% levels(data$condition)[c(1:3,5:6)] & compound == levels(data$compound)[3])
  
  idata <- read.csv("Data/140407_replicates_data.csv", stringsAsFactors = T) %>%
    dplyr::filter(strainID %in% levels(data$condition)[c(1:3,5:6)] & variable == levels(data$compound)[3]) %>%
    dplyr::select(strainID, variable, value) %>%
    dplyr::rename(condition = strainID, compound = variable)
  
  

  Figures[["S15C"]] <- ggplot(data, aes(x = condition, y = mean, group=condition))+
    geom_bar(aes(fill=condition), stat="identity", colour = "black")+
    geom_text(aes(label=signif(pvalue, 2), y=mean+sd+0.2), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = idata, aes(y = value), width = 0.2, colour = "grey")+
    geom_errorbar(data = data, aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,5))+
    ylab(bquote('GSSG [mM]'))+
    xlab("")+
    scale_x_discrete(labels=c(lb10,lb11, lb12, lb14, lb38))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(5)], brewer.pal(n = 8, name = palette)[c(2)], wt))
}

### Suppl. Fig. 15D
{
  data <- read.csv("Data/140407_replicates_summary.csv", stringsAsFactors = T)
  data$condition <- factor(data$condition, levels=levels(data$condition)[c(1,5,6,8,7,2,3,4)])
  
  data <- data %>%
    dplyr::filter(condition %in% levels(data$condition)[c(1:3,5:6)] & compound == levels(data$compound)[2])
  
  idata <- read.csv("Data/140407_replicates_data.csv", stringsAsFactors = T) %>%
    dplyr::filter(strainID %in% levels(data$condition)[c(1:3,5:6)] & variable == levels(data$compound)[2]) %>%
    dplyr::select(strainID, variable, value) %>%
    dplyr::rename(condition = strainID, compound = variable)
  
  Figures[["S15D"]] <- ggplot(data, aes(x = condition, y = mean, group=condition))+
    geom_bar(aes(fill=condition), stat="identity", colour = "black")+
    geom_text(aes(label=signif(pvalue, 2), y=mean+sd+0.2), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = idata, aes(y = value), width = 0.2, colour = "grey")+
    geom_errorbar(data = data, aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,22))+
    ylab(bquote('GSH/GSSG'))+
    xlab("")+
    scale_x_discrete(labels=c(lb10,lb11, lb12, lb14, lb38))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(5)], brewer.pal(n = 8, name = palette)[c(2)], wt))
}

### Suppl. Fig. 15E
{
  data <- read.xlsx("Data/151110_DHE_Summary.xlsx", sheetName="summary")
  data$Strain <- factor(data$Strain)
  data$Strain <- factor(data$Strain, levels=levels(data$Strain)[c(3,1,2,4)])

  idata <- readxl::read_excel("Data/151110_DHE_Summary.xlsx", sheet="percentages") %>%
    dplyr::select(-1)
  
  ncells <- idata %>%
    dplyr::group_by(Strain) %>%
    dplyr::summarise(cells = sum(Cells), n = n())
  
  idata$Strain <- factor(idata$Strain)
  idata$Strain <- factor(idata$Strain, levels=levels(idata$Strain)[c(3,1,2,4)])
  levels(idata$Strain)[4] <- "wt+tertBOOH"
  
  Figures[["S15E"]] <- ggplot(data, aes(x = Strain, y = mean))+
    geom_bar(stat="identity", aes(fill=Strain, group=Strain), colour="black")+
    geom_text(aes(label=signif(pvalue, 2), y=mean+sd+1), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_jitter(data = idata, aes(y = DHE_percent), colour = "grey", width = 0.2)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,45))+ 
    scale_x_discrete(breaks=levels(data$Strain), labels=c(lb10, lb11, lb13, lb31))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)], wt), breaks=levels(data$Strain), labels=c(lb10, lb11, lb13, lb31))+
    ylab("number of DHE-positive cells [%]")+
    xlab("")
}

### Suppl. Fig. 16
{
  data <- read.csv("Data/141119_Rapamycin_data_summary.csv", stringsAsFactors = T)
  fit <- read.csv("Data/141119_Rapamycin_fit_summary.csv", stringsAsFactors = T)
  EC50 <- read.csv("Data/141119_Rapamycin_summary.csv", stringsAsFactors = T)
  data$StrainID <- factor(data$StrainID, levels=levels(data$StrainID)[c(3,1,2)])
  fit$StrainID <- factor(fit$StrainID, levels=levels(fit$StrainID)[c(3,1,2)])
  names(EC50)[1] <- "StrainID"
  EC50$StrainID <- factor(EC50$StrainID, levels=levels(EC50$StrainID)[c(3,1,2)])
  levels(EC50$StrainID) <- paste("YSBN11", levels(EC50$StrainID))
  
  Figures[["S16"]] <- ggplot(fit, aes(x = concentration, group=StrainID))+
    #geom_point(data=data, aes(y = mean, colour=StrainID))+
    geom_line(aes(y = mean, colour=StrainID))+
    geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd, fill=StrainID), alpha=0.2)+
    geom_vline(data=EC50[which(EC50$Parameter == "ED50:(Intercept)"),], aes(xintercept=mean, colour=StrainID), alpha=0.5, linetype=2)+
    scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
    ylab(bquote('growh rate [OD min'^-1~']'))+
    xlab("Rapamycin [µg/mL]")+
    annotation_logticks(sides="b")+
    theme_classic()+
    theme(legend.justification=c(0,0), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,0), 
          legend.direction="vertical",
          legend.title=element_blank(),
          #legend.key=element_rect(colour = "black"),
          legend.background=element_blank())+
    scale_colour_manual(breaks=levels(data$StrainID), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))+
    scale_fill_manual(breaks=levels(data$StrainID), labels=c(lb10b,lb11b,lb13b), values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))
  
  
}

### Suppl. Fig. 17B
{
  data <- read.csv("Data/140227_PTS1_results.csv", stringsAsFactors = T)
  names(data)[1] <- "Strain"
  data <- cbind(data, genotype=c("rho0", "rho0","wild type"))
  data <- cbind(data, plasmid=c("pRS313", "pRS313_ATP3-6", "pRS313"))
  data$Strain <- factor(data$Strain, levels=levels(data$Strain)[c(3,1,2)])
  data$genotype <- factor(data$genotype, levels=levels(data$genotype)[2:1])
  data <- cbind(data, signif=NA)
  data$signif[which(data$p.value < 0.01)] <- "*"
  
  idata <- read.csv("Data/140227_PTS1_results_raw.csv", stringsAsFactors = T) %>%
    dplyr::rename(Strain = Metadata_genotype)
  
  Figures[["S17B"]] <- ggplot(data, aes(x = Strain, y = mean))+
    geom_bar(stat="identity", aes(fill=Strain, group=Strain), colour="black")+
    geom_text(aes(label=signif(p.value, 2), y=mean+sem+1), position=dodge, angle=90, hjust=0, vjust=0.5, size=6)+
    geom_text(aes(label = paste("n =", n)), angle = 90, fontface = "bold", y = 1, hjust=0, vjust=0.5, size=6)+
    geom_text(data = data[data$Strain == "WT_313",], aes(label = paste("n =", n)), angle = 90, fontface = "bold", y = 1, hjust=0, vjust=0.5, size=6, colour = "white")+
    #geom_jitter(data = idata, aes(y = Children_PXS_Count))+
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0.2)+
    theme_classic()+
    theme(legend.justification=c(1,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="none", 
          legend.direction="vertical",
          legend.title=element_blank())+
    scale_y_continuous(expand = c(0,0), limits=c(0,22))+ 
    scale_x_discrete(breaks=levels(data$Strain), labels=c(lb10, lb11, lb13))+
    scale_fill_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]), breaks=levels(data$Strain), labels=c(lb10, lb11, lb13))+
    xlab("")+
    ylab("number of peroxisomes per cell")+
    xlab("")
}

### Suppl. Fig. 17C
{
  fcsummary <- read.delim("Data/20210816_CombinedData_FoldChange_Summarized.tsv")
  signif <- read.delim("Data/20210816_CombinedData_Signif.tsv")
  poi <- c("ACO1", "CIT2", "DLD3", "IDH1", "IDH2", "PYC1", "RTG2")
  sdata <- fcsummary[which(fcsummary$Compound %in% poi),]
  ssignif <- signif[signif$Compound %in% poi,]
  ssignif$p.value <- signif(ssignif$p.value, 2)
  ssignif$plabel <- NA
  ssignif$plabel[ssignif$p.value < 1E-5] <- "p < 1E-5"
  ssignif$plabel[!ssignif$p.value < 1E-5] <- paste0("p = ", ssignif$p.value[!ssignif$p.value < 1E-5])
  sdata$Compound <- factor(sdata$Compound)
  sdata$Compound <- factor(sdata$Compound, levels = poi)
  
  
  poisdata <- fcdata[fcdata$Compound %in% poi,]
  sdata$sgdName <- paste0(str_to_title(tolower(sdata$Compound)), "p")
  
  Figures[["S17C"]] <- ggplot(sdata, aes(x = Compound, y = log2(mean), group=Strain, fill = Strain))+
    geom_bar(stat="identity", aes(fill=Strain), colour="black", position=dodge)+
    geom_errorbar(aes(ymin=log2(mean-sd), ymax=log2(mean+sd)), width=0.2, position=dodge)+
    geom_point(data = poisdata, aes(y = log2(foldchange)), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    geom_text(data = ssignif, aes(label = signif(p.value, 2), y = 2.0), angle = 90, hjust = 0.5, position = dodge)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position="none", 
          legend.direction="vertical",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          legend.title=element_blank())+
    ylab("log2 fold change")+
    ylim(-0.8, 2.5)+
    geom_hline(yintercept=0, colour="grey")+
    xlab("")+
    geom_text(aes(label=sgdName, y=min(log2(sdata$mean-sdata$sd))-0.05), position=dodge, angle=90, hjust=1, vjust=0.5, size=6)+
    #scale_x_discrete(breaks=levels(Pdata$Sample), labels=c(lb10, lb11, lb13))+
    scale_fill_manual(breaks=unique(sdata$Strain), labels=c(lb10,lb11,lb13), values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))
}

### Suppl. Fig. 18
{
  zsummary <- read.delim("Data/20210816_CombinedData_ZScores_Summarized.tsv")
  
  outfile <- try(readLines("Data/Metabolites_map_template_v2.svg"))
  poi <- zsummary %>%
    filter(Strain %in% c("rho+_313", "rho0_313", "rho0_ATP3_T911A"))
  poi$Strain <- factor(poi$Strain)
  poi$Strain <- factor(poi$Strain, levels = levels(poi$Strain))
  
  # make code
  cde <-
    data.frame(
      Compound = unique(odata$Compound),
      Experiment = "proteomics",
      Key = capitalize_str(tolower(paste0(
        sub(";.+", "", unique(odata$Compound)), "p"
      ))),
      style = "st2"
    )
  
  cde$Key[cde$Compound == "YGR012W"] <- "Mcy1p"
  
  code <- read.delim("Data/Metabolites_map_code.xls")
  code <- rbind(code, cde)
  
  maxval <- 2
  colrange <- seq(-maxval, maxval, by= 0.01)
  coloring_ <- data.frame(value=round(colrange,2), color=bluered(length(colrange)))
  coloring <- data.frame(value=round(colrange, 2))
  coloring <- merge(coloring, coloring_, by="value", all.x=T)
  coloring$color[is.na(coloring$color) & coloring$value < 0] <- "#0000FF"
  coloring$color[is.na(coloring$color) & coloring$value > 0] <- "#FF0000"
  
  ColorKey <- ggplot(coloring, aes(x = value, y = 1))+
    geom_bar(aes(fill = factor(value)), stat = "identity", width = 0.1)+
    xlab("z-score")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.ticks.x = element_line(colour = "black",size = 1, linetype = "solid"))+
    scale_fill_manual(values = as.character(coloring$color))+
    scale_y_continuous(expand = c(0,0), limits = c(0,1))
  
  pdf("Data/Figure_S18_key.pdf")
  print(ColorKey)
  dev.off()
  
  for(i in unique(poi$Compound)){
    
    spoi <- gsub(",", "", gsub(".", "", i, fixed=T), fixed=T)
    key <- code$Key[code$Compound == i]
    styletem <- code$style[code$Compound == i]
    
    
    if(i %in% poi$Compound[!is.na(poi$mean)] & 
       length(styletem) > 0 & 
       length(unique(which(grepl(paste0(key, "<"), outfile)))) > 0){
      # get data
      t <- poi[poi$Compound == i & !is.na(poi$mean),] %>%
        dplyr::group_by(Compound, Strain) %>%
        dplyr::summarise(value = mean(mean), .groups = "drop") %>%
        arrange(rev(Strain))
      
      t$value[(t$value) < min(coloring$value)] <- min(coloring$value)
      t$value[(t$value) > max(coloring$value)] <- max(coloring$value)
      
      # construct styles
      style <- list()
      style[[1]] <- paste("	.st", spoi,"evo{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho0_ATP3_T911A"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      style[[2]] <- paste("	.st", spoi,"rho{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho0_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      style[[3]] <- paste("	.st", spoi,"wt{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho+_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      
      names(style) <- c(paste0("st", spoi,"evo"), paste0("st", spoi,"rho"), paste0("st", spoi,"wt"))
      
      # add styles to style table
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[1]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[2]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[3]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      
      # apply style to element
      
      for(j in unique(which(grepl(paste0(key, "<"), outfile)))){
        rng <- which(grepl(styletem, outfile, fixed=T))
        rng <- rng[rng > j][1:3]
        
        for(k in 1:3){
          outfile[rng[k]] <- sub(styletem, names(style)[k], outfile[rng[k]])
        }
      }
    } else {
      print(i)
    }
  }
  outfile <- outfile[!is.na(outfile)]
  
  writeLines(text=outfile, con="Plots/Figure_S18.svg", sep="\n")
  
}

### Suppl. Fig. 19A
{
  zsummary <- read.delim("Data/20210816_CombinedData_ZScores_Summarized.tsv")
  
  # make code
  cde <-
    data.frame(
      Compound = unique(zsummary$Compound),
      Experiment = "proteomics",
      Key = capitalize_str(tolower(paste0(
        sub(";.+", "", unique(zsummary$Compound)), "p"
      ))),
      style = "st16"
    )
  cde$Key[cde$Compound == "YGR012W"] <- "Mcy1p"
  
  outfile <- try(readLines("Data/TCA_map_template_prot-01.svg"))
  poi <- zsummary %>%
    filter(Strain %in% c("rho+_313", "rho0_313", "rho0_ATP3_T911A"))
  poi$Strain <- factor(poi$Strain)
  poi$Strain <- factor(poi$Strain, levels = levels(poi$Strain))
  
  
  # code <- fcdata %>%
  #   distinct(Experiment, Compound) %>%
  #   mutate(Key = NA) %>%
  #   write.table(., "C:/Dropbox (Personal)/AG Ralser_JV/Publications/rho0_project/Figures/data/Metabolites_map_code.xls", sep="\t", row.names=F)
  
  code <- read.delim("Data/TCA_map_code.xls")
  code <- rbind(code, cde)
  
  maxval <- 2
  colrange <- seq(-maxval, maxval, by= 0.01)
  coloring_ <- data.frame(value=round(colrange,2), color=bluered(length(colrange)))
  coloring <- data.frame(value=round(colrange, 2))
  coloring <- merge(coloring, coloring_, by="value", all.x=T)
  coloring$color[is.na(coloring$color) & coloring$value < 0] <- "#0000FF"
  coloring$color[is.na(coloring$color) & coloring$value > 0] <- "#FF0000"
  
  ColorKey <- ggplot(coloring, aes(x = value, y = 1))+
    geom_bar(aes(fill = factor(value)), stat = "identity", width = 0.1)+
    xlab("z-score")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.ticks.x = element_line(colour = "black",size = 1, linetype = "solid"))+
    scale_fill_manual(values = as.character(coloring$color))+
    scale_y_continuous(expand = c(0,0), limits = c(0,1))
  
  pdf("Plots/Figure_S19A_key.pdf")
  print(ColorKey)
  dev.off()
  
  
  for(i in unique(poi$Compound)){
    
    spoi <- gsub(",", "", gsub(".", "", i, fixed=T), fixed=T)
    key <- code$Key[code$Compound == i]
    styletem <- code$style[code$Compound == i]
    
    
    if(i %in% poi$Compound[!is.na(poi$mean)] & length(styletem) > 0){
      # get data
      t <- poi[poi$Compound == i & !is.na(poi$mean),] %>%
        dplyr::group_by(Compound, Strain) %>%
        dplyr::summarise(value = mean(mean), .groups = "drop") %>%
        arrange(rev(Strain))
      
      t$value[(t$value) < min(coloring$value)] <- min(coloring$value)
      t$value[(t$value) > max(coloring$value)] <- max(coloring$value)
      
      # construct styles
      style <- list()
      style[[1]] <- paste("	.st", spoi,"evo{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho0_ATP3_T911A"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      style[[2]] <- paste("	.st", spoi,"rho{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho0_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      style[[3]] <- paste("	.st", spoi,"wt{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho+_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      
      names(style) <- c(paste0("st", spoi,"evo"), paste0("st", spoi,"rho"), paste0("st", spoi,"wt"))
      
      # add styles to style table
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[1]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[2]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[3]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      
      # apply style to element
      
      for(j in unique(which(grepl(paste0(key, "<"), outfile)))){
        rng <- which(grepl(styletem, outfile, fixed=T))
        rng <- rng[rng > j][1:3]
        
        for(k in 1:3){
          outfile[rng[k]] <- sub(styletem, names(style)[k], outfile[rng[k]])
        }
      }
    } else {
      print(i)
    }
  }
  outfile <- outfile[!is.na(outfile)]
  
  writeLines(text=outfile, con="Plots/Figure_S19A.svg", sep="\n")
  
}

### Suppl. Fig. 19B
{
  zsummary <- read.delim("Data/20210816_CombinedData_ZScores_Summarized.tsv")
  
  # make code
  cde <-
    data.frame(
      Compound = unique(zsummary$Compound),
      Experiment = "proteomics",
      Key = capitalize_str(tolower(paste0(
        sub(";.+", "", unique(zsummary$Compound)), "p"
      ))),
      style = "st16"
    )
  cde$Key[cde$Compound == "YGR012W"] <- "Mcy1p"
  
  outfile <- try(readLines("Data/TCA_map_template_met-01.svg"))
  poi <- zsummary %>%
    filter(Strain %in% c("rho+_313", "rho0_313", "rho0_ATP3_T911A"))
  poi$Strain <- factor(poi$Strain)
  poi$Strain <- factor(poi$Strain, levels = levels(poi$Strain))
  
  
  # code <- fcdata %>%
  #   distinct(Experiment, Compound) %>%
  #   mutate(Key = NA) %>%
  #   write.table(., "C:/Dropbox (Personal)/AG Ralser_JV/Publications/rho0_project/Figures/data/Metabolites_map_code.xls", sep="\t", row.names=F)
  
  code <- read.delim("Data/TCA_map_code.xls")
  code <- rbind(code, cde)
  
  maxval <- 2
  colrange <- seq(-maxval, maxval, by= 0.01)
  coloring_ <- data.frame(value=round(colrange,2), color=bluered(length(colrange)))
  coloring <- data.frame(value=round(colrange, 2))
  coloring <- merge(coloring, coloring_, by="value", all.x=T)
  coloring$color[is.na(coloring$color) & coloring$value < 0] <- "#0000FF"
  coloring$color[is.na(coloring$color) & coloring$value > 0] <- "#FF0000"
  
  
  for(i in unique(poi$Compound[poi$Experiment %in% c("TCA", "PPP")])){
    
    spoi <- gsub(",", "", gsub(".", "", i, fixed=T), fixed=T)
    key <- code$Key[code$Compound == i]
    styletem <- code$style[code$Compound == i]
    
    
    if(i %in% poi$Compound[!is.na(poi$mean)] & length(styletem) > 0){
      # get data
      t <- poi[poi$Compound == i & !is.na(poi$mean),] %>%
        dplyr::group_by(Compound, Strain) %>%
        dplyr::summarise(value = mean(mean), .groups = "drop") %>%
        arrange(rev(Strain))
      
      t$value[(t$value) < min(coloring$value)] <- min(coloring$value)
      t$value[(t$value) > max(coloring$value)] <- max(coloring$value)
      
      # construct styles
      style <- list()
      style[[1]] <- paste("	.st", spoi,"evo{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho0_ATP3_T911A"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      style[[2]] <- paste("	.st", spoi,"rho{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho0_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      style[[3]] <- paste("	.st", spoi,"wt{fill:", coloring$color[coloring$value == round((t$value[t$Strain == "rho+_313"]), 2)], ";stroke:#000000;stroke-miterlimit:10;}", sep="")
      
      names(style) <- c(paste0("st", spoi,"evo"), paste0("st", spoi,"rho"), paste0("st", spoi,"wt"))
      
      # add styles to style table
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[1]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[2]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      outfile <- c(outfile[1:which(grepl("<style type=", outfile))], style[[3]], outfile[which(grepl("<style type=", outfile))+1:length(outfile)])
      
      # apply style to element
      
      for(j in unique(which(grepl(paste0(key, "<"), outfile)))){
        if(j > 240){
          rng <- which(grepl(styletem, outfile, fixed=T))
          rng <- rng[rng > j][1:3]
          
          for(k in 1:3){
            outfile[rng[k]] <- sub(styletem, names(style)[k], outfile[rng[k]])
          }
        }}
    } else {
      print(i)
    }
  }
  
  outfile <- outfile[!is.na(outfile)]
  
  writeLines(text=outfile, con="Plots/Figure_S19B.svg", sep="\n")
  
}

### Suppl. Fig. 20A
{
  file <- "Data/labeling_data.xlsx" 
  sheet <- "slimmed"
  data <- readxl::read_xlsx(file,sheet = sheet,col_types = "text")
  
  data2 <- data %>% 
    mutate(labeling_percent = as.numeric(labeling_percent)) %>%
    mutate(time = as.numeric(time)) %>%
    mutate(abundance_nmol = as.numeric(abundance_nmol))
  
  temp <- data2 %>% filter(tracer == "glucose")
  
  temp <- temp %>% filter(metabolite %in% c("pyruvate", "citrate", "cis-aconitate", "glutamate", "gaba", "glutamine","succinate")) # "a-ketoglutarate", "fumarate", "malate"
  
  #remove outlier by sample (metabolites clearly affected in one sample)
  temp <- temp %>% filter(!(tracer == "glucose" & time == 30 & genotype == "rho0 313" & replicate == 2))
  
  temp2 <- temp %>% group_by(tracer, metabolite, time, genotype) %>%
    dplyr::mutate(mean_labeling_percent = mean(labeling_percent, na.rm = T)) %>%
    dplyr::mutate(sd_labeling_percent = sd(labeling_percent, na.rm = T)) %>%
    
    dplyr::mutate(mean_abundance_nmol = mean(abundance_nmol, na.rm = T)) %>%
    dplyr::mutate(sd_abundance_nmol = sd(abundance_nmol, na.rm = T)) %>%
    ungroup() 
  
  temp2$genotype <- factor(temp2$genotype, levels = c("wt 313", "rho0 313", "rho0 T911A"))
  
  # plotting summary
  # labelling - facet wrap
  Figures[["S20A"]] <- temp2 %>% 
    ggplot(aes(x=time, y=labeling_percent, colour = genotype)) +
    geom_point() + 
    geom_line(aes(y=mean_labeling_percent)) +
    geom_errorbar(aes(ymin=mean_labeling_percent - sd_labeling_percent, ymax=mean_labeling_percent + sd_labeling_percent)) + 
    facet_wrap(~ metabolite, scales = "free_y", ncol = length(unique(temp2$metabolite))) +
    scale_colour_manual(values = c("wt 313" = wt, "rho0 313" = rho, "rho0 T911A" = brewer.pal(n = 8, name = palette)[c(1)])) +
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.direction="vertical",
          legend.background = element_blank(),
          strip.background = element_blank(),
          legend.title=element_blank())+
    ylab("metabolite labelling (%)")+
    xlab("")
}

### Suppl. Fig. 20B
{
  file <- "Data/labeling_data.xlsx" 
  sheet <- "slimmed"
  data <- readxl::read_xlsx(file,sheet = sheet,col_types = "text")
  
  data2 <- data %>% 
    mutate(labeling_percent = as.numeric(labeling_percent)) %>%
    mutate(time = as.numeric(time)) %>%
    mutate(abundance_nmol = as.numeric(abundance_nmol))
  
  temp <- data2 %>% filter(tracer == "glucose")
  
  temp <- temp %>% filter(metabolite %in% c("pyruvate", "citrate", "cis-aconitate", "glutamate", "gaba", "glutamine","succinate")) # "a-ketoglutarate", "fumarate", "malate"
  
  #remove outlier by sample (metabolites clearly affected in one sample)
  temp <- temp %>% filter(!(tracer == "glucose" & time == 30 & genotype == "rho0 313" & replicate == 2))
  
  temp2 <- temp %>% group_by(tracer, metabolite, time, genotype) %>%
    dplyr::mutate(mean_labeling_percent = mean(labeling_percent, na.rm = T)) %>%
    dplyr::mutate(sd_labeling_percent = sd(labeling_percent, na.rm = T)) %>%
    
    dplyr::mutate(mean_abundance_nmol = mean(abundance_nmol, na.rm = T)) %>%
    dplyr::mutate(sd_abundance_nmol = sd(abundance_nmol, na.rm = T)) %>%
    ungroup() 
  
  temp2$genotype <- factor(temp2$genotype, levels = c("wt 313", "rho0 313", "rho0 T911A"))
  
  # plotting summary
  # labelling - facet wrap
  Figures[["S20B"]] <- temp2 %>% 
    ggplot(aes(x=time, y=abundance_nmol, colour = genotype)) +
    geom_point() + 
    geom_line(aes(y=mean_abundance_nmol)) +
    geom_errorbar(aes(ymin=mean_abundance_nmol - sd_abundance_nmol, ymax=mean_abundance_nmol + sd_abundance_nmol)) + 
    facet_wrap(~ metabolite, scales = "free_y", ncol = length(unique(temp2$metabolite))) +
    scale_colour_manual(values = c("wt 313" = wt, "rho0 313" = rho, "rho0 T911A" = brewer.pal(n = 8, name = palette)[c(1)])) +
    theme_classic()+
    theme(legend.justification=c(0,1), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position=c(0,1), 
          legend.direction="vertical",
          legend.background = element_blank(),
          strip.background = element_blank(),
          legend.title=element_blank())+
    ylab("metabolite abundance (nmol)")+
    xlab("")
}

### Suppl. Fig. 21
{
  data <- data.table::fread("Data/20210816_IronData_FoldChange_Data.tsv")
  matrix <- cast(data, Sample ~ sgdName, value="foldchange")
  matrix$Sample <- factor(matrix$Sample)
  levels(matrix$Sample) <- sub("Fe_", "Fe ", levels(matrix$Sample))
  matrix <- cbind(t(as.data.frame(strsplit(as.character(matrix$Sample), " "))), matrix)
  names(matrix)[1:2] <- c("Strain", "Replicate")
  
  for(i in 4:ncol(matrix)){
    if(T %in% is.na(matrix[,i])){
      t <- aggregate(matrix[,i], by=list(matrix$Strain), FUN=mean, na.rm=T)
      for(j in unique(matrix$Strain)){
        matrix[which(matrix$Strain == j),i][is.na(matrix[which(matrix$Strain == j),i])] <- t$x[which(t$Group.1 == j)]
      }
    }
  }
  matrix <- na.omit(matrix)
  
  matrix.pca <- prcomp(matrix[,-1:-3], scale. = T)
  
  matrix.class <- factor(na.omit(matrix$Strain))
  matrix.class <- factor(matrix.class, levels=levels(matrix.class)[c(6,2,4,5,1,3)])
  
  #averaged GO variables
  GOI <- read.delim("Data/GOI_Summary.tsv")
  GOI$sgdName <-  paste0(str_to_title(tolower(GOI$sgdName)), "p")

  t <- cbind(sgdName=row.names(matrix.pca$rotation), matrix.pca$rotation)
  GOI <- merge(GOI, t[,], by="sgdName", all.x=T)
  for(i in 3:ncol(GOI)){
    GOI[,i] <- as.numeric(as.character(GOI[,i]))
  }
  GOI <- GOI[!duplicated(GOI[,1:2]) | !duplicated(GOI[,1:2], fromLast=T),]
  select <- which(row.names(matrix.pca$rotation) %in% GOI$sgdName[which(GOI$GO == unique(GOI$GO)[1])])
  goi <- GOI[,1:2]
  
  GOI <- aggregate(GOI[,3:ncol(GOI)], by=list(GOI$GO), FUN=mean, na.rm=T)
  row.names(GOI) <- GOI$Group.1
  GOI <- GOI[,-1]
  
  Figures[["S21"]] <- ggbiplot(matrix.pca, 
                  groups=matrix.class,
                  var.scale=0.5,
                  ellipse = T,
                  var.alpha=0.2,
                  select = select,
                  custom.var.palette = c("black", "darkred"),
                  custom.var=T,
                  custom.var.linetype=2,
                  custom.varname.size = 4,
                  custom.var.size = 1,
                  custom.var.data=GOI)+
    theme_classic() +
    theme(legend.justification=c(0,0), 
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position="bottom", 
          legend.direction="vertical",
          legend.background = element_blank(),
          legend.title=element_blank())+
    #scale_color_manual(values=c(wt, rho, brewer.pal(n = 8, name = palette)[c(1)]), labels=c(lb10b,lb11b,lb13b))+
    geom_hline(yintercept=0, alpha=0.5)+
    geom_vline(xintercept=0, alpha=0.5)+
    xlim(-2,2)+
    ylim(-1.5,1.5)
}

### Suppl. Fig. 22A
{
  fcsummary <- read.delim("Data/20210816_IronData_FoldChange_Summarized.tsv")
  signif <- read.delim("Data/20210816_CombinedData_Signif.tsv")
  
  GO <- data.frame()
  for(i in list.files(path = "Data/", pattern = "_annotations.txt", full.names = T)){
    t <- read.csv(i, sep="\t", skip=8, stringsAsFactors = T)
    GO <- rbind(GO, t)
  }
  
  GO <- GO %>%
    distinct(Gene, Gene.Ontology.Term) %>%
    dplyr::rename(Compound = Gene)
  levels(GO$Gene.Ontology.Term)
  
  GOI <- read.csv("Data/Iron_GOI.tsv", sep="\t")
  GOI <- GOI[which(!GOI$GO %in% c("CIA", "Heme")),]
  GOI$GO <- factor(GOI$GO)
  levels(GOI$GO)[1:6] <- "Iron transport"
  levels(GOI$GO)[2] <- "Iron-sulphur proteins"
  GOI <- GOI[!duplicated(GOI$sgdName),]
  
  poi <- unique(GO$Compound[GO$Gene.Ontology.Term %in% levels(GO$Gene.Ontology.Term)[c(1,2,3,6,7)]])
  sdata <- fcsummary[which(fcsummary$Compound %in% poi),]
  ssignif <- signif[signif$Compound %in% poi,]
  ssignif$p.value <- signif(ssignif$p.value, 2)
  ssignif$plabel <- NA
  ssignif$plabel[ssignif$p.value < 1E-5] <- "p < 1E-5"
  ssignif$plabel[!ssignif$p.value < 1E-5] <- paste0("p = ", ssignif$p.value[!ssignif$p.value < 1E-5])
  poisdata <- fcdata[fcdata$Compound %in% poi,]
  sdata$sgdName <- paste0(str_to_title(tolower(sdata$Compound)), "p")
  
  sdata <- sdata %>%
    filter(!sgdName %in% c("Lip5p"))
  
  poisdata <- poisdata %>%
    filter(!sgdName %in% c("Lip5p"))
  
  sdata <- sdata %>%
    mutate(Strain = factor(Strain, levels = c("WT Control - Iron", "Rho0 Control - Iron", "Rho0 ATP3-T911A - Iron", "WT Control + Iron", "Rho0 Control + Iron", "Rho0 ATP3-T911A + Iron")))
  
  
  Figures[["S22A"]] <- ggplot(sdata, aes(x = sgdName, y = log2(mean), fill = Strain))+
    geom_bar(stat="identity", aes(fill=Strain), colour="black", position=dodge)+
    geom_errorbar(aes(ymin=log2(mean-sd), ymax=log2(mean+sd)), width=0.2, position=dodge)+
    geom_point(data = poisdata, aes(y = log2(foldchange)), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    #geom_text(data = ssignif[ssignif$p.value < 0.01,], aes(label = plabel, y = 2.5), angle = 90, hjust = 0.5, position = dodge)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position="none", 
          legend.direction="vertical",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          legend.title=element_blank())+
    ylab("log2 fold change")+
    ylim(min(log2(sdata$mean-sdata$sd))-0.8, 3)+
    geom_hline(yintercept=0, colour="grey")+
    xlab("")+
    geom_text(aes(label=sgdName, y=min(log2(sdata$mean-sdata$sd))-0.05), position=dodge, angle=90, hjust=1, vjust=0.5, size=6)+
    #scale_x_discrete(breaks=levels(Pdata$Sample), labels=c(lb10, lb11, lb13))+
    scale_fill_manual(breaks=unique(sdata$Strain), values=c(wtnoIron, rhonoIron, evonoIron, wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))
}

### Suppl. Fig. 22B
{
  fcsummary <- read.delim("Data/20210816_IronData_FoldChange_Summarized.tsv")
  signif <- read.delim("Data/20210816_CombinedData_Signif.tsv")
  
  GOI <- read.csv("Data/Iron_GOI.tsv", sep="\t")
  GOI <- GOI[which(!GOI$GO %in% c("CIA", "Heme")),]
  GOI$GO <- factor(GOI$GO)
  levels(GOI$GO)[1:6] <- "Iron transport"
  levels(GOI$GO)[2] <- "Iron-sulphur proteins"
  GOI <- GOI[!duplicated(GOI$sgdName),]
  
  poi <- GOI$sgdName[GOI$GO == "Iron-sulphur proteins"]
  sdata <- fcsummary[which(fcsummary$Compound %in% poi),]
  ssignif <- signif[signif$Compound %in% poi,]
  ssignif$p.value <- signif(ssignif$p.value, 2)
  ssignif$plabel <- NA
  ssignif$plabel[ssignif$p.value < 1E-5] <- "p < 1E-5"
  ssignif$plabel[!ssignif$p.value < 1E-5] <- paste0("p = ", ssignif$p.value[!ssignif$p.value < 1E-5])
  poisdata <- fcdata[fcdata$Compound %in% poi,]
  sdata$sgdName <- paste0(str_to_title(tolower(sdata$Compound)), "p")
  
  sdata <- sdata %>%
    filter(!sgdName %in% c("Lip5p"))
  
  poisdata <- poisdata %>%
    filter(!sgdName %in% c("Lip5p"))

  sdata <- sdata %>%
    mutate(Strain = factor(Strain, levels = c("WT Control - Iron", "Rho0 Control - Iron", "Rho0 ATP3-T911A - Iron", "WT Control + Iron", "Rho0 Control + Iron", "Rho0 ATP3-T911A + Iron")))
  
  Figures[["S22B"]] <- ggplot(sdata, aes(x = sgdName, y = log2(mean), fill = Strain))+
    geom_bar(stat="identity", aes(fill=Strain), colour="black", position=dodge)+
    geom_errorbar(aes(ymin=log2(mean-sd), ymax=log2(mean+sd)), width=0.2, position=dodge)+
    geom_point(data = poisdata, aes(y = log2(foldchange)), colour = "grey", position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05))+
    #geom_text(data = ssignif[ssignif$p.value < 0.01,], aes(label = plabel, y = 2.5), angle = 90, hjust = 0.5, position = dodge)+
    theme_classic()+
    theme(legend.justification=c(0,1), 
          legend.position="none", 
          legend.direction="vertical",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          legend.title=element_blank())+
    ylab("log2 fold change")+
    ylim(min(log2(sdata$mean-sdata$sd))-0.8, 3)+
    geom_hline(yintercept=0, colour="grey")+
    xlab("")+
    geom_text(aes(label=sgdName, y=min(log2(sdata$mean-sdata$sd))-0.05), position=dodge, angle=90, hjust=1, vjust=0.5, size=6)+
    #scale_x_discrete(breaks=levels(Pdata$Sample), labels=c(lb10, lb11, lb13))+
    scale_fill_manual(breaks=unique(sdata$Strain), values=c(wtnoIron, rhonoIron, evonoIron, wt, rho, brewer.pal(n = 8, name = palette)[c(1)]))
}

### Suppl. Fig. 25A
{
  # data
  fcsummary <- read.delim("Data/20210816_AcoData_FoldChange_Summarized.tsv")
  fcsummary$Compound <- paste0(str_to_title(tolower(fcsummary$Compound)), "p")
  pdata <- read.delim("Data/20210816_AcoData_Raw_Data.tsv")
  
  #PCA analysis
  odata <- pdata %>%
    dplyr::select(R.FileName, PG.Genes, PG.Quantity, ID) %>%
    dplyr::rename(Sample = R.FileName, Compound = PG.Genes, value = PG.Quantity)
  
  odata$Compound <- paste0(str_to_title(tolower(odata$Compound)), "p")
  
  matrix <- cast(odata, Sample ~ Compound, value="value")
  
  matrix <- cbind(t(as.data.frame(strsplit(as.character(matrix$Sample), "_"))), matrix)
  matrix <- cbind(Strain=paste(matrix[,3], matrix[,4], matrix[,5], sep="_"), matrix[,c(-1,-2,-3, -4, -5, -7)])
  names(matrix)[2] <- c("Replicate")
  
  for(i in 3:ncol(matrix)){
    if(T %in% is.na(matrix[,i])){
      t <- aggregate(matrix[,i], by=list(matrix$Strain), FUN=mean, na.rm=T)
      for(j in unique(matrix$Strain)){
        matrix[which(matrix$Strain == j),i][is.na(matrix[which(matrix$Strain == j),i])] <- t$x[which(t$Group.1 == j)]
      }
    }
  }
  matrix <- na.omit(matrix)
  
  matrix.pca <- prcomp(matrix[,-1:-2], scale. = T)
  loadings <- data.frame(matrix.pca$rotation)
  
  #PC1
  loadings <- loadings[order(-abs(loadings$PC1)),]
  goi <- row.names(loadings)[1:50]
  
  subdata <- fcsummary %>%
    filter(Compound %in% goi) %>%
    mutate(sgdName = factor(Compound)) %>%
    mutate(StrainID = factor(Strain, levels = c("WT Control",      "WT ATP3-T911A", "Rho0 Control",    "Rho0 ATP3-T911A", "aco1 Control",    "aco1 ATP3-T911A"  )))
  
  
  arrdata <- subdata %>%
    dplyr::group_by(sgdName) %>%
    dplyr::summarise(meanval = mean(log2(mean))) %>%
    dplyr::arrange(-meanval)
  
  subdata$sgdName <- factor(subdata$sgdName, levels = unique(arrdata$sgdName))
  
  Figures[["S25A"]] <- ggplot(data = subdata, aes(x = sgdName, y = StrainID))+
    geom_raster(aes(fill = log2(mean)))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))
}

### Suppl. Fig. 25B
{
  # data
  fcsummary <- read.delim("Data/20210816_AcoData_FoldChange_Summarized.tsv")
  fcsummary$Compound <- paste0(str_to_title(tolower(fcsummary$Compound)), "p")
  pdata <- read.delim("Data/20210816_AcoData_Raw_Data.tsv")
  
  #PCA analysis
  odata <- pdata %>%
    dplyr::select(R.FileName, PG.Genes, PG.Quantity, ID) %>%
    dplyr::rename(Sample = R.FileName, Compound = PG.Genes, value = PG.Quantity)
  
  odata$Compound <- paste0(str_to_title(tolower(odata$Compound)), "p")
  
  matrix <- cast(odata, Sample ~ Compound, value="value")
  
  matrix <- cbind(t(as.data.frame(strsplit(as.character(matrix$Sample), "_"))), matrix)
  matrix <- cbind(Strain=paste(matrix[,3], matrix[,4], matrix[,5], sep="_"), matrix[,c(-1,-2,-3, -4, -5, -7)])
  names(matrix)[2] <- c("Replicate")
  
  for(i in 3:ncol(matrix)){
    if(T %in% is.na(matrix[,i])){
      t <- aggregate(matrix[,i], by=list(matrix$Strain), FUN=mean, na.rm=T)
      for(j in unique(matrix$Strain)){
        matrix[which(matrix$Strain == j),i][is.na(matrix[which(matrix$Strain == j),i])] <- t$x[which(t$Group.1 == j)]
      }
    }
  }
  matrix <- na.omit(matrix)
  
  matrix.pca <- prcomp(matrix[,-1:-2], scale. = T)
  loadings <- data.frame(matrix.pca$rotation)
  
  #PC1
  loadings <- loadings[order(-abs(loadings$PC2)),]
  goi <- row.names(loadings)[1:50]
  
  subdata <- fcsummary %>%
    dplyr::filter(Compound %in% goi) %>%
    dplyr::mutate(sgdName = factor(Compound)) %>%
    dplyr::mutate(StrainID = factor(Strain, levels = c("WT Control",      "WT ATP3-T911A", "Rho0 Control",    "Rho0 ATP3-T911A", "aco1 Control",    "aco1 ATP3-T911A"  )))
  
  
  arrdata <- subdata %>%
    dplyr::group_by(sgdName) %>%
    dplyr::summarise(meanval = mean(log2(mean))) %>%
    dplyr::arrange(-meanval)
  
  subdata$sgdName <- factor(subdata$sgdName, levels = unique(arrdata$sgdName))
  
  Figures[["S25B"]] <- ggplot(data = subdata, aes(x = sgdName, y = StrainID))+
    geom_raster(aes(fill = log2(mean)))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))
}

### Suppl. Fig. 25C
{
  # data
  fcsummary <- read.delim("Data/20210816_AcoData_FoldChange_Summarized.tsv")
  fcsummary$Compound <- paste0(str_to_title(tolower(fcsummary$Compound)), "p")
  pdata <- read.delim("Data/20210816_AcoData_Raw_Data.tsv")
  
  #PCA analysis
  odata <- pdata %>%
    dplyr::select(R.FileName, PG.Genes, PG.Quantity, ID) %>%
    dplyr::rename(Sample = R.FileName, Compound = PG.Genes, value = PG.Quantity)
  
  odata$Compound <- paste0(str_to_title(tolower(odata$Compound)), "p")
  
  matrix <- cast(odata, Sample ~ Compound, value="value")
  
  matrix <- cbind(t(as.data.frame(strsplit(as.character(matrix$Sample), "_"))), matrix)
  matrix <- cbind(Strain=paste(matrix[,3], matrix[,4], matrix[,5], sep="_"), matrix[,c(-1,-2,-3, -4, -5, -7)])
  names(matrix)[2] <- c("Replicate")
  
  for(i in 3:ncol(matrix)){
    if(T %in% is.na(matrix[,i])){
      t <- aggregate(matrix[,i], by=list(matrix$Strain), FUN=mean, na.rm=T)
      for(j in unique(matrix$Strain)){
        matrix[which(matrix$Strain == j),i][is.na(matrix[which(matrix$Strain == j),i])] <- t$x[which(t$Group.1 == j)]
      }
    }
  }
  matrix <- na.omit(matrix)
  
  matrix.pca <- prcomp(matrix[,-1:-2], scale. = T)
  loadings <- data.frame(matrix.pca$rotation)
  
  #PC1
  loadings <- loadings[order(-abs(loadings$PC3)),]
  goi <- row.names(loadings)[1:50]
  
  subdata <- fcsummary %>%
    dplyr::filter(Compound %in% goi) %>%
    dplyr::mutate(sgdName = factor(Compound)) %>%
    dplyr::mutate(StrainID = factor(Strain, levels = c("WT Control",      "WT ATP3-T911A", "Rho0 Control",    "Rho0 ATP3-T911A", "aco1 Control",    "aco1 ATP3-T911A"  )))
  
  
  arrdata <- subdata %>%
    dplyr::group_by(sgdName) %>%
    dplyr::summarise(meanval = mean(log2(mean))) %>%
    dplyr::arrange(-meanval)
  
  subdata$sgdName <- factor(subdata$sgdName, levels = unique(arrdata$sgdName))
  
  Figures[["S25C"]] <- ggplot(data = subdata, aes(x = sgdName, y = StrainID))+
    geom_raster(aes(fill = log2(mean)))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))
}

## Plotting----

for(i in names(Figures)){
  Dimensions[[i]] <- c(6,4)
}

# PCA
Dimensions[["5A"]] <- c(18.6,4.9)
Dimensions[["5B"]] <- c(18.6,4.9)
Dimensions[["5C"]] <- c(2.3, 4.5)
Dimensions[["5J"]] <- c(8.3,8.3)
Dimensions[["2B"]] <- c(4.9,4.9)
Dimensions[["2D"]] <- c(4.9,4.9)
Dimensions[["3D"]] <- c(4.6,5.2)
Dimensions[["4C"]] <- c(5.0,4.6)
Dimensions[["4D"]] <- c(5.0,4.6)
Dimensions[["5F"]] <- c(3.6,5.2)
Dimensions[["5E"]] <- c(2.3, 4.5)
Dimensions[["5I"]] <- c(3.6,5.2)
Dimensions[["5H"]] <- c(2.3, 5.2)
Dimensions[["S13A"]] <- c(4.9, 5.3)
Dimensions[["S13B"]] <- c(4.9, 5.3)
Dimensions[["S13C"]] <- c(6, 5.3)
Dimensions[["S15B"]] <- c(4.9, 4.7)
Dimensions[["S15C"]] <- c(4.9, 4.7)
Dimensions[["S15D"]] <- c(4.9, 4.7)
Dimensions[["S15E"]] <- c(4.9, 4.7)
Dimensions[["S17B"]] <- c(4.4, 4.4)


for(i in names(Figures)){
  pdf(paste0("Plots/Figure_", i, ".pdf"), width = Dimensions[[i]][1], height = Dimensions[[i]][2])
  print(Figures[[i]])
  dev.off()
}

## End----