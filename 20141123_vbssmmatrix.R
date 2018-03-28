###WORKING FUNCTIONS
vbssmmatrix <- function(HLonly=FALSE, treat='Light'){
  cat('This function takes a tab-delimited genelist of a single column with header,\n overlaps the genelist with the timeseries \n and creates a matrix ready to use for VBSSM modelling\n\n')
  cat('If HLonly=TRUE, only the genes that respond exclusively to high light (and not to temperature) are selected\nn')
  cat('If treat="Light", the dataset "Light" from Hlcombo2 is selected. If treat="Control", the dataset "Control" is selected instead\nn')
  cat('It also creates a CATMAprobe_to_TAIR mapping file to use as node attributes with Cytoscape')
  
  #Read the data to make the matrix for VBSSM
  HSF <- lapply(list.files(pattern='mylist.txt'), read.delim)
  names(HSF) <- lapply(HSF, function(a) colnames(a))
  HSF <- lapply(HSF, unique)
  
  #For HLonly genes
  func.extract <- function(a,b) {
    a1<- droplevels(a[a[,1] %in% b[,1],1])
    a1<- droplevels(a1[a1 %in% timeseries[,2]])
    return(a1)
  }
  
  #For all genes
  func.extract2 <- function(a){
    a1 <- droplevels(a[a[,1] %in% timeseries[,2],])
    return(a1)
  }
  
  #summary
  l1<- lapply(HSF, function(b) lapply(list('HLonly'=HL.only.m, 'Temponly'=temp.only.m, 'TempHLinfl'=Temp.HLinfl.m), function(a) summary(a[a[,1] %in% b[,1],] %in% timeseries[,2])))
  #genelists
  capture.output(cat('Overlaps of HSFdatasets with HL.only, temp.only and temp.HL-influenced genes within the HL timeseries. \n FALSE shows the number of no overlaps, TRUE the number of overlaps \n\n'), l1, file='overlaps_of_genelists_with_timeseries.txt')
  
  #If HLonly==TRUE, Retrieve HL.only genelists, and add the respective HSFA factor to the dataset. If FALSE, overlap all the genes in the 'mylist' with the timeseries
  if(HLonly==TRUE){
    l2<- lapply(HSF, function(b) lapply(list('HLonly'=HL.only.m, 'Temponly'=temp.only.m, 'TempHLinfl'=Temp.HLinfl.m), func.extract,b))
    ##Generation of matrices for VBSSM modeling
    l2.vssm <- lapply(l2, function(a) a[[1]])
    names(l2.vssm) <- paste0(names(l2.vssm), '.HLonly')
    #extract the genes from l2.vssm and convert to matrix
  } else l2.vssm <- lapply(HSF, func.extract2)
  
  l2.vssm2 <- lapply(l2.vssm, function(a) as.matrix(HLcombo2[HLcombo2[,2] %in% a,]))
  

  #IMPORTANT! Extract either the 'Light' or the 'Control' treatments. The VBSSM cannot handle data with both treatments simultaneously
  
  func.matrix <- function(a, treat) {
    ifelse(treat=='Light',a <- a[,c(1,2,55:106)], a <- a[, c(1:54)])
    repnames <- do.call(c,lapply(strsplit(colnames(a), '\\.'), function(x) x[3]))
    reptimepoints <- gsub('X', '', do.call(c,lapply(strsplit(colnames(a), '\\.'), function(x) x[1])))
    reptimepoints[1:2] <- NA
    a <- rbind(repnames, reptimepoints, a)
    return(a)
  }
  #VSSM files with the HSFs
  l2.vssm3 <- lapply(l2.vssm2, func.matrix, treat)
  ifelse(HLonly==TRUE, f.name <- 'VBSSM_HLonly_', f.name <- 'VBSSM_Allgenes_')
  ifelse(treat=='Light', f.name <- paste0(f.name, 'Light_'), f.name <- paste0(f.name, 'Control_'))
  mapply(function(a,b) write.table(a, paste0(f.name, b, '.txt'), na='', quote=F, row.names=F, col.names=F, sep='\t'), l2.vssm3, names(HSF))
  
  #map the CATMA probes to TAIR to use in Cytoscape
  write.table(Catma.map, 'CATMAmapping_forCytoscape.txt', sep='\t', row.names=FALSE, quote=F)
}











#######Original script############

##Create files for VBSSM modelling
HL30min <- read.delim('mergewithcontrastmatrix.txt')

#LogFC.x == genes DE in HL treament
#LogFC.y == genes DE in temp (27 C)
#LogFC == genes DE between HL and temp treatments

#LogFC.x  LogFC.y   LogFC
#DE       DE        NA -temp_only genes
#DE       DE        DE -HL_influenced, temp-responsive genes
#DE       NA        DE -HL_only responsive genes

#Map ATG codes to the probes
map <- read.delim('probe_to_name_mapping.txt') #the mapping file was created from TAIR
HL30min <- merge(HL30min[,-3], map, by='ID')

#Clean the merged file
HL30min <- HL30min[HL30min$Name!='no_match',]
HL30minsplit <- HL30min[grep(';', HL30min$Name),]
HL30min <- HL30min[-grep(';', HL30min$Name),]
spl <- do.call(rbind, strsplit(as.character(HL30minsplit$Name), ";", fixed=T))[,-3] #it throws a column number error because some probes bind to 2 or 3 genes. I deleted the 3rd column because they happened to be chloroplastic genes.
HL30minsplit <- cbind(HL30minsplit[,-15], spl)
a <- HL30minsplit[,c(1:15,16)]
colnames(a)[16] <- 'Name'
b <- HL30minsplit[,c(1:15,17)]
colnames(b)[16] <- 'Name'

HL30min <- rbind(HL30min,a,b)
merindex <- HL30min[, c(15, 3, 7, 11)]

#Select X_only genes
temp.only <- HL30min[apply(merindex[,c(2:4)], 1, function(a) all(is.na(a)==c(F,F,T))==T),]
HL.only <- HL30min[apply(merindex[,c(2:4)], 1, function(a) all(is.na(a)==c(F,T,F))==T),]
Temp.HLinfl <- HL30min[apply(merindex[,c(2:4)], 1, function(a) all(is.na(a)==c(F,F,F))==T),]
HL.only <- HL.only[!(HL.only$Name %in% Temp.HLinfl$Name),] #PROBLEM: there are different probes for the same gene, and in some cases one says the gene is HL.only and the other says the gene is Temp.HLinfl. The solution I have opted for is stringent: I put those genes as Temp.HLinfl

lapply(list(temp.only, HL.only, Temp.HLinfl), function(a) summary(duplicated(a$Name))) #check for gene duplicates

#Remove duplicates, retain the probe with highest absolute expression value
func.order <- function(a,b) {
  a <- a[order(a$Name, -a[,b]),]
  a <- a[!duplicated(a$Name),]
}
temp.only <- func.order(temp.only, "logFC.y")
HL.only <- func.order(HL.only, "logFC.x")
Temp.HLinfl <- func.order(Temp.HLinfl, "logFC")

###
#Create X_only files for merging for the 'Selecting genes for VBSSM model with HSF'
temp.only.m <- data.frame('Temp.only'=droplevels(temp.only[,15]))
HL.only.m <- data.frame('HL.only'=droplevels(HL.only[,15]))
Temp.HLinfl.m <- data.frame('HLinfl.T-responsive'=droplevels(Temp.HLinfl[,15]))
#double-check for duplicates
summary(duplicated(do.call(c, lapply(list(temp.only.m, HL.only.m, Temp.HLinfl.m), function(a) levels(a[,1])))))

#Delete unwanted variables
rm(list=c('a', 'b', 'HL30minsplit', 'merindex'))

#Read the data to make the matrix for VBSSM
HSF <- lapply(list.files(pattern='mylist.txt'), read.delim)
names(HSF) <- lapply(HSF, function(a) colnames(a))

#TimeSeries DE genes
timeseries <- read.delim('Clustered Genesco2.txt')[,c(1,2)]

##Check number of overlaps between the HL and temp datasets and the HSF datasets
#not overlapped with the timeseries
lapply(HSF, function(b) lapply(list('HLonly'=HL.only.m, 'Temponly'=temp.only.m, 'TempHLinfl'=Temp.HLinfl.m), function(a) summary(a[,1] %in% b[,1])))
#overlapped with the timeseries

func.extract <- function(a,b) {
  a1<- droplevels(a[a[,1] %in% b[,1],1])
  a1<- droplevels(a1[a1 %in% timeseries[,2]])
  return(a1)
}

#summary
l1<- lapply(HSF, function(b) lapply(list('HLonly'=HL.only.m, 'Temponly'=temp.only.m, 'TempHLinfl'=Temp.HLinfl.m), function(a) summary(a[a[,1] %in% b[,1],] %in% timeseries[,2])))
#genelists
l2<- lapply(HSF, function(b) lapply(list('HLonly'=HL.only.m, 'Temponly'=temp.only.m, 'TempHLinfl'=Temp.HLinfl.m), func.extract,b))
#genelists with mapped genenames and descriptions
l3<- lapply(l2, function(a) lapply(a, function(b) merge(b, TAIRmapping[,-1], by=1, all.x=TRUE, sort=FALSE)))

capture.output(cat('Overlaps of HSFdatasets with HL.only, temp.only and temp.HL-influenced genes within the HL timeseries. \n FALSE shows the number of no overlaps, TRUE the number of overlaps \n\n'), l1)

##Generation of matrices for VBSSM modeling
#Retrieve HL.only genelists, and add the respective HSFA factor to the dataset
l2.vssm <- lapply(l2, function(a) a[[1]])
names(l2.vssm) <- paste0(names(l2.vssm), '.HLonly')
#get the HighLight-Combo2
HLcombo2 <- read.delim('HighLight-Combo2.txt')
#extract the genes from l2.vssm and convert to matrix
l2.vssm2 <- lapply(l2.vssm, function(a) as.matrix(HLcombo2[HLcombo2[,2] %in% a,]))

#IMPORTANT! Extract only the 'Light' treatments. The VBSSM could not handle the file otherwise

func.matrix <- function(a) {
  a <- a[,c(1,2,55:106)]
  repnames <- do.call(c,lapply(strsplit(colnames(a), '\\.'), function(x) x[3]))
  reptimepoints <- gsub('X', '', do.call(c,lapply(strsplit(colnames(a), '\\.'), function(x) x[1])))
  reptimepoints[1:2] <- NA
  a <- rbind(repnames, reptimepoints, a)
  return(a)
}
#VSSM files with the HSFs
l2.vssm3 <- lapply(l2.vssm2, func.matrix)
mapply(function(a,b) write.table(a, paste0('VSSM_HLonly_', b, '.txt'), na='', quote=F, row.names=F, col.names=F, sep='\t'), l2.vssm3, )

#map the CATMA probes to TAIR to use in Cytoscape
Catma.map <- merge(HLcombo2[,1:2], TAIRmapping, by=2, sort=FALSE, all.x=TRUE)[-3]
colnames(Catma.map)[c(1,2)] <- c('TAIR', 'ID')
#Catma.map[,2] <- toupper(Catma.map[,2])
Catma.map <- Catma.map[, c(2,1,3,4)]
write.table(Catma.map, 'CATMAmapping.txt', sep='\t', row.names=FALSE, quote=F)
##############END ORIGINAL SCRIPT#############

