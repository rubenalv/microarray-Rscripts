#v25: added extra page size for genelists with >30 elements
#     changed the cutoff for the pheatmap graphs to 1e-5, deleting the 1e-3
#     exchanged 'letters' by 'numbers' (vector c.numbers) in the overlapping calculations, to allow for more than 26 lists (=lenght(letters))
#     corrected the pheatmaps to both have '-log10(phyper + 1e-100)' values (the timeseries pheatmap was +1e-50)

#v27. 20150225: expression profiles are plotted by ATG or CATMA probe entered directly on the command line (to be able to choose which probes are in the models. Lists of CATMA probes are not supported yet)
#     fixed bug with the expression profile heatmap (y-axis names were ordered alphabetically, and that messed cluster order) 

#SCRIPT USAGE
#First of all, load the libraries required typing "startit()".
#Use <Interrogate.timeseries(..., ATG=TRUE, use.mylists=TRUE, plot.phyper=TRUE, plot.profile=FALSE, v.plotname="geneplot.pdf", stringent=FALSE, p.threshold=1)
##(...) Search terms.
      #Write an AGI code, e.g. Interrogate.timeseries('at4g01570'), and it will return the cluster and expression profile of the gene (no phyper plots for a single gene).
      #Alternatively, write a CATMA probe code (e.g. Interrogate.timeseries('CATMA3A21000'). The script takes either AGI's or CATMA's, but not both, and if you enter both only the AGI's will be read in.
			#LEAVE IT EMPTY if 'mylist.txt' lists with AGI codes are in the working folder (supply tab-delimited files with the string 'mylist' in them. 
				#The files can have one or several columns with data, and the first column has to contain the AGI codes.
				#Columns must have a header, and it is recommended that the AGI column header is the name of the experiment, e.g. 'HL'.
				#Any duplicated AGI codes or ';'-separated codes are removed by the script and not used to compute the phyper)
			#WRITE ONE OR MORE TERMS you want to look for in the GO database. 
				#Terms are quoted and comma-separated: 'auxin', 'cytokinin', 'heat'.
##ATG=TRUE						#If TRUE, search for ATG codes and printout of expression profiles is allowed.
##use.mylists=TRUE	#If TRUE, it reads 'mylist.txt' files from the working directory.
##plot.phyper=TRUE		#If TRUE, a plot with the log10 hypergeometric probabilities and a plot with the proportion of genes in the list supplied with respect to the number of genes in the cluster is produced.
##plot.proportions=FALSE #It requires plot.phyper=TRUE. If TRUE, a plot with the proportion of genes in the list supplied with respect to the number of genes in the cluster is produced.
##plot.profile=FALSE	#If TRUE, a plot with all the expression profiles of genes matching to the timeseries database, and individual plots per cluster are produced.
##v.plotname=".pdf"	  #The name given to the plots. It has to be quoted and the file extension is required, and can be '.pdf', '.png', '.jpg' or '.tif'.
##p.threshold=1				#Maximum threshold probability (not in log10) to be plotted, both for the phyper plots and the expression profile plots. E.g. p.threshold=0.001, no clusters with P>0.001 will be plotted or profiled.
##stringent=TRUE			#If more than one list or search term is provided and stringent=TRUE, only those clusters which are under the threshold are plotted and profiled.
##create.timeseriestxt=F #It creates a file called 'TimeSeries_mylist.txt' in the working folder. This file can then be read by the main function to calculate pairwise overlaps of genes including the timeseries

#################
#Interrogate.GOdatabase(...)                              #Finds a term (or list) in the GOterms of the full database and produces a list of AGI codes. E.g. find all the genes related to 'auxin', 'immun' (note that partial matching is allowed with grep), etc.
#Interrogate.timeseries(...)   														#Main function
#Phyper.this(...)                    										  #Normal phyper
#Expression.profile(...)                                  #Expression profile plotting function

###########PACKAGE REQUIREMENTS#########
startit <- function(){
	func.require <- function(a) {
		if(!require(a, character.only=TRUE) == TRUE) {
			cat("... trying to install package ", a, "...")
			install.packages(a, repos="http://www.stats.bris.ac.uk/R")}
		require(a, character.only=TRUE)
}

    packlist <- list('ggplot2', 'gridExtra', 'reshape2', 'pheatmap', 'scales')
    lapply(packlist, func.require)
cat('########Script version 23 20140530########\n\n\n')
cat('#USE <Interrogate.timeseries(..., ATG=TRUE, use.mylists=TRUE, plot.phyper=TRUE, plot.heatmap=TRUE, heatmap.filename="heatmap.pdf", plot.profile=FALSE, v.plotname="geneplot.pdf", stringent=FALSE, p.threshold=1, phyper.pairwise=FALSE, pairwise.only=FALSE)\n\n')
cat('##(...) SEARCH TERMS.\n')
cat('\t\t#Write an AGI code, e.g. Interrogate.timeseries("at4g01570"), and it will return the cluster and expression profile of the gene (no phyper plots for a single gene).\n')
cat('\t\t#LEAVE IT EMPTY if *mylist.txt* lists with AGI codes are in the working folder (supply tab-delimited files with the string *mylist* in them.\n')
cat('\t\t#The files can have one or several columns with data, and the first column has to contain the AGI codes.\n')
cat('\t\t#Columns must have a header, and it is recommended that the AGI column header is the name of the experiment, e.g. "HL".\n')
cat('\t\t#Any duplicated AGI codes or', ';','-separated codes are removed by the script and not used to compute the phyper).\n')
cat('\t\t#WRITE ONE OR MORE TERMS you want to look for in the GO database.\n')
cat('\t\t#Terms are quoted and comma-separated: "auxin", "cytokinin", "heat".\n\n')
cat('##ATG=TRUE							#If TRUE, search for ATG codes and printout of expression profiles is allowed.\n')	
cat('##use.mylists=TRUE			#If TRUE, the "mylist" files are taken from the working folder and appended to the GO database search terms.\n')
cat('##remove.chloromito=TRUE  #If TRUE, it removes chloroplast (ATCxxxxx) and mitochondrial (ATMxxxxx) identifiers from the genelists.\n\n')
cat('##plot.phyper=TRUE		  #If TRUE, a plot with the log10 hypergeometric probabilities and a plot with the proportion of genes in the list supplied with respect to the number of genes in the cluster is produced.\n\n')
cat('##plot.proportions=FALSE #It requires plot.phyper=TRUE. If TRUE, a plot with the proportion of genes in the list supplied with respect to the number of genes in the cluster is produced.\n\n')
cat('##plot.heatmap=TRUE  	#If TRUE, it plots a heatmap with the log10 hypergeometric probabilities. "Cluster 0" is the full list of genes in the timeseries, not divided in clusters.\n\n')
cat('##plot.heatmap.pergene=FALSE    #If TRUE, it plots a heatmap with the overlaps of the genes on the TimeSeries.\n\n')
cat('##heatmap.filename=".pdf" #The name given to the heatmap. It has to be quoted and the file extension is required, and can be ".pdf", ".png", ".jpg" or ".tif".\n\n')
cat('##plot.profile=FALSE	  #If TRUE, a plot and a heatmap with all the expression profiles of genes matching to the timeseries database, and individual plots per cluster are produced. If ="heatmap", only the heatmap is produced\n\n')
cat('##v.plotname=".pdf"	  #The name given to the plots. It has to be quoted and the file extension is required, and can be ".pdf", ".png", ".jpg" or ".tif".\n\n')
cat('##stringent=FALSE  		#If more than one list or search term is provided and stringent=TRUE, only those clusters which are under the threshold are plotted and profiled.\n\n')
cat('##p.threshold=1				#Maximum threshold probability (not in log10) to be plotted, both for the phyper plots and the expression profile plots. E.g. p.threshold=0.001, no clusters with P>0.001 will be plotted or profiled.\n\n')
cat('##phyper.pairwise=FALSE  #If TRUE, the GOterm lists or mylists are compared pairwise as well.\n\n')	
cat('##pairwise.only=FALSE  #If TRUE, the script runs only the pairwise comparisons between lists and does not do the timeseries. If pairwise.only=TRUE, phyper.pairwise is automatically set to TRUE.\n\n\n')

cat('##On the phyperplot bar graph and the results:\n')
cat('#"n.genes" is the number of genes in the cluster (or full timeseries.\n')
cat('#"prob.T" is the probability (p-value) of the list of genes in the cluster and the test list to be similar.\n')
cat('#"overlaps" is the number of genes that overlap between the cluster (of full timeseries list) and the test list of genes.\n')
cat('#"proportion" is the "proportion (%) of overlapping genes, i.e. 2*number.of.overlapping.genes*100/(genes.in.list.A + genes.in.list.B).\n\n\n')

cat('##FILES CREATED:\n')
	cat('#results_X_in_AGIdatabase.txt... Being X a search term for the GO database (e.g. "heat"), this list contains all the genes that have this GO descriptor and overlap with the timeseries\n')
  cat('#RES_phyper_results_probabilitiesonly.txt... the probatilities of the phyper analysis.\n')
	cat('#RES_phyper_allresults.txt... the full phyper results file (number of genes per cluster, overlaps, probability and proportion of overlapping genes per cluster for each dataset.\n')
	cat('#RES_mergedAGIs.txt... file that merges the timeseries AGI codes with the AGI codes of each of the datasets.\n')
	cat('#RES_mergedAGIs_listsonly.txt... Same as above, but one the genes in the timeseries that match any of the datasets are shown.\n')
	cat('#pairwise_merge.txt... merged list for pairwise comparisons of lists of genes.\n')
	cat('#pairwise_phyper.txt... full phyper results file for pairwise comparisons (ID1 and ID2 are the names of the genelists, n.genes is the total number of genes in each list, overlaps are the genes that two lists have in common, prob the likelihood (as p-value) of those to lists being similar, and proportion the proportion of genes in common as defined above: (2*number of overlaps*100)/sum(n.genes1 + n.genes2)))\n')
  cat('#.\n')
	
}
########################################
#######UPDATE FUNCTION##########
#Use this function to update the GO database from TAIR
update.GO <- function() {
  download.file('ftp://ftp.arabidopsis.org/home/tair/Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt', 'GOslim.txt')
  GOslim <- read.delim("GOslim.txt", header=FALSE)[,c(1,4,5,6,9)]
  colnames(GOslim) <- c('TAIR', 'Relationship', 'Key', 'GO', 'Key2')
  GOslim <<- unique(GOslim) #For some reason they have a good number of duplicated entries
  save.image(paste0(Sys.Date(),'_Timeseries_phyper and expression profiling.RData'))
}
################################



Interrogate.GOdatabase <- function(...) {
  v.index <- lapply(list(...), function(a) unique(c(grep(a, GOslim$Key, ignore.case=TRUE), grep(a, GOslim$Key2, ignore.case=TRUE))))
	names(v.index) <- c(...)
	v.nozero <- do.call(c, lapply(v.index, function(a) length(a)!=0)) #indexing to eliminate the values with no matches.
	v.index <- v.index[v.nozero]
	if(any(v.nozero==FALSE)) cat('Some of your search terms do not have a match in the GO database: ', paste(c(...)[!v.nozero], collapse=', '),'.\n\n')
	result <- lapply(v.index, function(a) a<-droplevels(GOslim[a,]))
  
  func.write <- function(a,b) write.table(a, paste0("results_",b, "_in_AGIdatabase.txt"), sep='\t', row.names=FALSE)
  mapply(func.write, result, names(result))
  result <- lapply(result, function(a) a <- a[!duplicated(a[,1]),])
  for(i in 1:length(names(result))) colnames(result[[i]])[1] <- 'GO.'
  for(i in 1:length(names(result))) colnames(result[[i]]) <- paste0(colnames(result[[i]]),names(result)[i])
  names(result) <- paste0('GO.', names(result))
  return(result)
}

func.checkmyfiles <- function(remove.chloromito) {
	filelist <- list.files(pattern="mylist", ignore.case=TRUE)
	files <- lapply(filelist, read.delim)
  names(files) <- do.call(c, lapply(files, function(a) colnames(a)[1]))
	files <- lapply(files, as.data.frame)
	files <- lapply(files, func.clean,remove.chloromito)
	return(files)
}

func.clean <- function(a,remove.chloromito){
  a[,1] <- toupper(a[,1])
  mylist <- a
  dups <- duplicated(mylist[,1])
  if(!all(dups==FALSE)){
    cat("WARNING! There were --",summary(dups)[[3]], "-- duplicates in file ",colnames(mylist)[1], ", out of --", nrow(mylist), "-- entries.", sep="")
    mylist <- mylist[!duplicated(mylist[,1]),]
    if(class(mylist)=='character') {mylist <- data.frame(mylist)
                                    colnames(mylist) <- colnames(a)} #this is to prevent an error when the genelist has only one column. If that is the case, it is converted to character and the colnames are missed.
    cat("The duplicates were removed, and now", colnames(mylist)[1], "contains", summary(dups)[[2]], "unique genes.")
    cat("\n\n")
  }
  #Duplicate the columns in the lists to create a column called "ID" for merging.
  mylist <- cbind(mylist[1], mylist)
  colnames(mylist)[1] <- "ID"
  #Remove semicolons
  semicolon <- grepl(";", mylist[,1])
  if (!all(semicolon==FALSE)){
    cat("###WARNING! The file", colnames(mylist)[2], "has some entries that cannot be matched:\n",as.character(mylist[semicolon,1]), "\n")
    cat("Please use single TAIR names. Those entries have been DELETED before proceeding with the analysis.")
    cat('\n\n')
    mylist <- mylist[!semicolon,]
  }
  #Remove probe (Agilent) names, starting by "A_"
  probe <- grepl("A_", mylist[,1])
  if (!all(probe==FALSE)){
    cat("###WARNING! The file", colnames(mylist)[2], "has some Agilent probe codes that will not be matched:\n",mylist[probe,1], "\n")
    cat("Please use single TAIR names. Those entries have been DELETED before proceeding with the analysis.")
    cat('\n\n')
    mylist <- mylist[!probe,]
  }
  #Remove chroplastic and mitochondrial genes
  chloromito <- grepl("ATC|ATM", mylist[,1])
  if(remove.chloromito==TRUE) {
    if(any(chloromito)) cat("### Mitochondrial and/or chloroplastic gene ID's have been removed from list", colnames(mylist)[2],"\n", as.character(mylist[chloromito,1]),"\n\n")
    mylist <- mylist[!chloromito,]
  }
  
  #Final nomenclature check
  finalcheck <- grepl("AT.G", mylist[,1])
  char <- nchar(as.character(mylist[,1]))
  if(!all(finalcheck&(char==9)==TRUE)) {
    cat("### WARNING! The file", colnames(mylist)[2], "has some codes that do not match TAIR nomenclature 'ATxG.....' (being the dots 5 numbers):\n",mylist[!(finalcheck&char==9),1], "\n")
    cat("Only TAIR names are valid for the matching. Invalid values are counted in the hypergeometric calculations and will lead to wrong probabilities.\n")
    cat("Would you like to delete them (d) and continue, or to stop (s)?")
    if(interactive()) input <- scan(what=character(), n=1, quiet=TRUE)
    if(input=='d') mylist <- mylist[(finalcheck&char==9),] else stop('Terminated by user', call.=FALSE)
    cat('\n\n')
  }
   return(mylist)
}

Interrogate.timeseries <- function(..., ATG=TRUE, use.mylists=TRUE, remove.chloromito=TRUE, plot.phyper=TRUE, plot.proportions=FALSE, plot.heatmap=TRUE, plot.heatmap.pergene=FALSE, heatmap.filename='heatmap.pdf', plot.profile=FALSE, v.plotname="geneplot.pdf", stringent=FALSE, p.threshold=1, phyper.pairwise=FALSE, pairwise.only=FALSE, create.timeseriestxt=F){ #profile.all=FALSE implies that only the genes that went over a 0.001 threshold in phyper are plotted.

  if(create.timeseriestxt==T) write.table(data.frame('TimeSeries'=tclusters[,2]), 'TimeSeries_mylist.txt', sep='\t', row.names=F, quote=F)
  
  justATG=FALSE
  if(length(c(...))==0) {
  	files <- func.checkmyfiles(remove.chloromito)
  	if(length(files)==0) stop("Sorry, you have not supplied any 'mylist.txt' files", call.=FALSE) 
  	} else {
      if(ATG==TRUE) {
      	ATGgrep <- toupper(grep('\\<AT.G', c(...), ignore.case=TRUE, value=TRUE))
      	CATMAgrep <- grep('\\<CATMA', c(...),ignore.case=TRUE, value=TRUE)
          if(length(ATGgrep)!=0) {
          		ATG.check <- match(ATGgrep, tclusters$ID)
              #if(all(is.na(ATG.check))) stop(paste('Sorry, but no match for:',paste(ATGgrep, collapse=', '),'\n You can get the profile suppliying a _mylist.txt file instead\n'), call.=FALSE)
              if(any(is.na(ATG.check))) {
    						cat('Some of the entered AGI codes do not match the timeseries:', paste(ATGgrep[is.na(ATG.check)], collapse=', '))
    						readline(prompt="Press [enter] to continue")
    					  }
              justATG=TRUE
          	  plot.profile=TRUE}
        	if(length(CATMAgrep)!=0) {
          	  CATMA.check <- match(CATMAgrep, TimeSeries$probe)
          	  if(any(is.na(CATMA.check))) {
          	    cat('Some of the entered CATMA probes do not match the timeseries:', paste(CATMAgrep[is.na(CATMA.check)], collapse=', '))
          	    readline(prompt="Press [enter] to continue")
          	  }
          	  justATG=TRUE
          	  plot.profile=TRUE}
          if(all(ATGgrep==0, CATMAgrep==0)) {if(use.mylists == TRUE) files.1 <- func.checkmyfiles(remove.chloromito) else files.1 <- NULL
                                             files <- Interrogate.GOdatabase(...) #Formerly AGIlist
                                             files <- lapply(files, function(a) a <- cbind(ID=a[,1],a))
                                             if(!is.null(files.1)) files <- append(files, files.1)}
      } else {if(use.mylists == TRUE) files.1 <- func.checkmyfiles(remove.chloromito) else files.1 <- NULL
              files <- Interrogate.GOdatabase(...) #Formerly AGIlist
              files <- lapply(files, function(a) a <- cbind(ID=a[,1],a))
              if(!is.null(files.1)) files <- append(files, files.1)}
  	}
  
  if(pairwise.only==FALSE) {
          if(justATG==FALSE) {
            	files.out <- Phyper.this(files, plot.phyper, plot.proportions, plot.heatmap, plot.heatmap.pergene, heatmap.filename, stringent, p.threshold)
            	thresholed.clusters <- as.integer(files.out[[1]])
            	files.out <- as.character(files.out[[2]])} else {
            	  if(length(ATGgrep)!=0) files.out <- ATGgrep else if(length(CATMAgrep)!=0) files.out <- CATMAgrep
            	}
      	  
        if(plot.profile != FALSE) {
            gene <- rbind(droplevels(TimeSeries[TimeSeries$ID %in% files.out,]), droplevels(TimeSeries[TimeSeries$probe %in% files.out,])) #Extract the genes in the genelist from the time series. The 'rbind' is a workaround to not need to use and "if" statement for CATMA probe and ATG
          	#if(justATG==FALSE) gene <- droplevels(gene[gene$Cluster %in% thresholed.clusters, ])
          			
          	func.stderr <- function(a) sqrt(var(a,na.rm=TRUE)/length(na.omit(a)))
            func.aov <- function(input) data.frame('ID'=input$ID[1],'probe'=input$probe[1], 'pvalue'=summary(aov(value ~ va.treat, data=input))[[1]][[5]][[1]], 'Time'=input$Time[1])
            
            gene.agg <- data.frame(aggregate(gene[,'value'], as.list(gene[,c('ID', 'Cluster', 'probe', 'va.treat', 'Time')]), mean), aggregate(gene[,5], as.list(gene[,c(1,2,4,6,8)]), sd)[6], aggregate(gene[,5], as.list(gene[,c(1,2,4,6,8)]), func.stderr)[6])
          	colnames(gene.agg) <- c('ID', 'Cluster', 'probe', 'Treatment', 'Time', 'Mean', 'SD', 'SE')
          	f <- fullmap[!duplicated(fullmap$TAIR),]
          	gene.agg <- droplevels(merge(gene.agg, f, by.x='ID',by.y="TAIR", all.x=TRUE))
          	gene.agg <- gene.agg[order(gene.agg$Treatment, gene.agg$ID, gene.agg$probe, gene.agg$Time),]
            
            #Check if there are gene ID's not in the differentially expressed timeseries (tclusters$ID)
            if(!all(gene.agg$ID %in% tclusters$ID)){
              cat("###Your genelist contains genes (below) that are not differentially expressed:\n")
              cat(as.character(unique(gene.agg[!gene.agg$ID %in% tclusters$ID,'ID'])), '\n')
              cat("Would you like to plot all profiles (a) or just those in the timeseries (t)?")
              if(interactive()) input <- scan(what=character(), n=1, quiet=TRUE)
              if(input=='t') {gene.agg <- droplevels(gene.agg[gene.agg$ID %in% tclusters$ID,])
                              gene <- droplevels(gene[gene$ID %in% tclusters$ID,])}
            } 
                        
            ##ANOVA to add significance *'s to the expression profiles
            anova.s <- split(gene, f=list(gene$probe, as.factor(gene$Time)), drop=TRUE)
            anova.s <- data.frame(do.call(rbind, lapply(anova.s, func.aov)), 'thr'='')
            anova.s$thr <- as.character(anova.s$thr) #Create asterisks if p.value < 0.05
            anova.s[anova.s$pvalue < 0.05, 'thr'] <- '*'
                      
            #aggregate to take only the max value of control/light, which will be used to place the significance asterisks with a geom_text()
            a <- aggregate(gene.agg$Mean, by=list(gene.agg$ID, gene.agg$probe, gene.agg$Time), FUN=max)
            colnames(a) <- c('ID', 'probe', 'Time', 'maxv')
            a <- a[order(a$ID, a$probe, a$Time),]
            anova.s4 <- anova.s[order(anova.s$ID, anova.s$probe, anova.s$Time),]
            anova.s4 <- cbind(anova.s4, a)[, c(1:5,9)]
            anova.s4 <- anova.s4[order(anova.s4$ID, anova.s4$probe, anova.s4$Time),]
            gene.agg <- cbind(gene.agg, anova.s4)[, c(1:8, 17, 19, 20, 9:12, 14)]
            gene.aggexcel <- gene.agg[,c(1,13,2:6)]
            gene.aggexcel <- dcast(gene.aggexcel, Cluster + ID + probe + Symbol + Treatment ~ Time, mean, value.var='Mean')
            write.table(gene.agg, file='gene_expression_values.txt', sep='\t', row.names=F)
            write.table(gene.aggexcel, file='gene_expression_values_hetmapexcel.txt', sep='\t', row.names=F)
            cat('Printing gene expression profiles...\n')
            
            #Expression values heatmap
            g2 <- gene.agg[,c(1:11,13)]
            #g2$Cluster <- as.factor(g2$Cluster)
            #g2$Time <- as.factor(g2$Time)
            #g2$Treatment <- as.factor(g2$Treatment)
            #g2$thr <- as.character(g2$thr)
            index<- g2$Treatment=='Control'
            g2[, "y"] <- factor(paste0(g2[,"Cluster"],' ',g2[,"ID"],'_',g2[,"probe"],'_',g2[,"Symbol"],'_',g2[,"Treatment"]))
            g2[index, "z"] <-paste0(g2[index,"Cluster"],' ',g2[index,"ID"],'_',g2[index,"probe"],'_',g2[index,"Symbol"],"_Control")
            g2[!index, "z"] <-"HL"
            g2 <- g2[order(g2$Cluster, g2$ID, g2$probe, g2$Time,g2$Treatment),]
            sortlev <- rev(unique(as.numeric(g2$y))) #note the reverse order (rev)
            g2$y <- factor(g2$y, levels=levels(g2$y)[sortlev])
            g2$Treatment <- factor(g2$Treatment, levels=c("Light", "Control"))
            g2labels <- rev(g2[!duplicated(g2$y),'z'])
            
            gg <- ggplot(g2, aes(x=Time, y=y)) +
              scale_y_discrete(labels=g2labels) +
              geom_raster(aes(fill=Mean)) + 
              geom_text(data=g2[g2$Treatment=='Light',],aes(label=thr), vjust=0.8) +
              scale_fill_gradient2(low='blue', high='red') +
              xlab('Time (h)') +
              ylab(NULL) +
              theme(panel.background=element_blank(),
                    axis.text=element_text(colour='black'))
            ggsave('heatmap_geneexpressionprofiles.pdf',gg, width=8, height=1+(length(g2labels)/6), limitsize=F)
            
            if(plot.profile != 'heatmap'){
                #Linegraphs expression values
                gene.agg2 <- split(gene.agg, gene.agg$Cluster)
              	gene.agg2 <- lapply(gene.agg2, droplevels)
              	names(gene.agg2) <- levels(gene.agg$Cluster)
                
                #scale_x_continuous(breaks=seq(0,6,0.5),labels=as.character(c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6))) +
                  
              	func.ggplot <- function(a1, v.plotname) {
              		nlev <- nlevels(a1$probe)
              		gg <- ggplot(a1, aes(x=Time, y=Mean)) +
              			geom_line(size=0.7, aes(group=Treatment:ID:probe, linetype=Treatment)) +
              			ggtitle(paste0("Cluster_", levels(a1$Cluster))) +
              			geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), size=0.4, width=0.03) +
              			ylab("Mean expression (log2)") + xlab("Time (h)") +
              		  #ylab("Mean expression (log2)") + xlab("Time (days)") + #This is for the drought timeseries
                    #scale_x_continuous(breaks=seq(0,13,1)) + #This is for the drought timeseries
              		  scale_x_continuous(breaks=seq(0,6,0.5),labels=as.character(c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6))) +
              		  facet_wrap(~ ID + probe + Symbol, scales='free_y', ncol=2) +
              		  theme(panel.background=element_blank(),
              		        panel.grid.major.y = element_line(colour='grey90'),
              		        panel.grid.minor.y = element_line(colour='grey90'),
              		        axis.text = element_text(colour='black')#,
              		        #legend.position = 'none'
                          )+
              		  geom_text(aes(x=Time-0.08, y=(maxv + max(maxv)*0.035), label=thr, group=probe))
              		 #use x=Time-0.1 for the drought timeseries
              		
              		if(nlev %% 2 == 0) h <- (nlev+2/9)*2.2 else h <- (nlev + 12/9)*2.2
                  if(nlev == 1) w <- 10 else w <- 16
              
                  gg1 <- list(gg, h, w)
              		ggsave(paste0("Cluster_", levels(a1$Cluster), "_", v.plotname), gg, width=w, height=h, limitsize=FALSE)
              		cat('.')
                  return(gg1)
              	}
          	          	
          	gg1 <- lapply(gene.agg2, func.ggplot, v.plotname)
          	geneplots <- lapply(gg1, function(a) ggplotGrob(a[[1]]))
            geneplots <- lapply(gg1, function(a) a[[1]])
            heights <- as.integer(do.call(c, lapply(gg1, function(a)a[[2]]))) #Extract the 'h' values defined by the func.ggplot
            cat('printing the multigene_plot.pdf...\n')
            
            pdf('multigene_plot.pdf', w=8, h=sum(heights))
          	do.call(grid.arrange, c(geneplots, clip=T, ncol=1, heights=list(heights))) #Trick here is that widths and heights MUST be lists, so do.call actually works with lists embedded in lists.
            dev.off()
           }
            
        }
  } else phyper.pairwise=TRUE
  if(phyper.pairwise==TRUE) func.phyperpair(files)
cat('... All done!')
}


Phyper.this <- function(files, plot.phyper=FALSE, plot.proportions=FALSE, plot.heatmap=TRUE, plot.heatmap.pergene=FALSE, heatmap.filename='heatmap.pdf', stringent=FALSE, p.threshold=1, phyper.pairwise=FALSE){
		## This script computes the phyper for the overlap of the clusters file and the user-supplied file.
		# In phyper(a,b,c,d,lower.tail=FALSE),
		# a = the number of genes in the overlap (in a mix of black and white balls, white balls drawn with no replacement)
		# b = the number of genes in the cluster (number of white balls)
		# c = the number of genes from the microarray pool thare are not in the cluster ( c = 24000 - b) (number of black balls)
		# d = the number of genes in the user-supplied list (balls drawn)
		#
		## This scripts creates several files:
		# - "RES_mergeAGI.txt" is the alignment of the clusters file with the user-supplied files. All the genes from the clusters file are present, and the genes from the test list that do not overlap are ommited.
		# - "RES_mergedAGIs_listsonly.txt" is the alignment above, showing only the overlapped AGIs.
		# - "RES_phyper_allresults.txt" is the results file, with a column corresponding to the cluster numbers, another to the number of genes that overlap between the clusters file
		#    and the user-supplied file, and the probability of the overlap given by the phyper distribution, for each of the user-supplied files.
		# - "RES_phyper_results_probabilitiesonly.txt". Same as above, but only with the column "cluster" and the probabilities for each of the user-supplied files. It also creates a replicate set of -log10 of the probabilities. 
		#
		## LIMITATIONS
		# - No limit of number of files, BUT they must have a filename "mylist"+number.txt" (e.g. mylist.txt, mylist1.txt, mylist2.txt)
		# - The lists have to be tab-delimited and have column names (generally, it is recommended that the column name of the TAIR column is e.g. the name of the experiment or the mutant).
		# - No limit of number of columns per file, BUT the first column has to contain the TAIR (AGI codes).
		# - The lists to test should NOT have repeated TAIR names. The scripts checks this, and deletes any repeated names before computing the phyper, notifying the user.
		# - Sometimes more that one gene hybridises to a probe, giving as results entries like "ATG145670; ATG243256". The script checks for the presence of ";" and if found it DELETES the reference.
		# - The user should make sure that TAIR names are in uppercase (ATG145590 and not Atg145590), since if they are in different case they are not recognised as equal. Anyway, the scripts converts all the columns to match into uppercase.
		#
		## ADDITIONAL USES
		# If the user creates a tab-delimited file called "Clustered Genesco2.txt", with columns "cluster"(only 1 cluster), "ID" and any other information
		# (e.g. mapping symbols names, expression levels, etc.), this script will create a merged file an compute the phyper of all the other files (myfile1..n)
		# with respect to that file.
		#
		## USAGE: you need
		# - the clusters file (called "Clustered Genesco2.txt", tab-delimited file)
		# - the lists you want to test.
		#
		## Open R, set up the working directory (File>Change dir...) to where the files are, and call this script (File>Source R Code...). 
		####
		#
		
		
		###########PARAMETERS###########
		#In case a particular order or the genes is required for the graphs, type it here.
		gene.levels <- c("LECTIN356", 'TCP19', 'BBX32', 'CRK13', 'AIG2L', 'NUDT7', 'WHY2', 'CGS1', 'HSFA7a', 'PAD4', 'ERF3', 'LOL1')
		use.gene.levels = FALSE
		
		
		###########FUNCTIONS############

		func.phyper <- function(clus,b) {
		  genes.in.cluster <- length(clus)
		  overlaps <- sum(clus %in% files.trim[[b]])
		  k <- length(files.trim[[b]])
		  n <- 24000 - genes.in.cluster
		  prob <- phyper(overlaps,genes.in.cluster,n,k,lower.tail=FALSE)
		  results <- cbind("n.genes"=genes.in.cluster, overlaps, prob, proportion=(2*overlaps*100/(genes.in.cluster+k)))
		  colnames(results)[3] <- c(paste0(names(files.trim)[b]))
		  return(results)                
		}
    

		###########SCRIPT BODY###################
		#Set the output digits to 5
		options(digits=5)
		#
		cat('..')
		numberofclusters <- nlevels(tclusters$Cluster)
		
		# Import the files from the working directory. Convert them in a list of dataframes.
		  files.out <- unique(do.call(c, lapply(files, function(a) as.character(a[[1]]))))
		  len.file <- length(files)
		
   	#Merge all the files together with the timeseries data, and sort the Clusters by name.
		mergefile1 <- tclusters
	  for(x in 1:len.file) mergefile1 <- merge(mergefile1, files[[x]], by="ID", all.x=TRUE, sort=FALSE)
		mergefile1 <- mergefile1[order(mergefile1$Cluster),]
    len.colmer1 <- length(colnames(mergefile1))-2
    if(len.colmer1==1) {
      mergefile2 <- mergefile1[!is.na(mergefile1[,-c(1,2)]),]} else {
		  mergefile2 <- mergefile1[!rowSums(is.na(mergefile1[,-c(1,2)]))==len.colmer1,]  #Delete the rows that are all NA's
      } #This is to prevent an error in rowSums when len.colmer1=1, i.e. when the array has only one column to make the sums.
	  ful <- fullmap[,c(1,6,7)]
		ful <- ful[!duplicated(ful[,1]),]
	  mergefile2 <- merge(mergefile2, ful, by=1, all.x=TRUE, sort=FALSE)
  
    ##Data for the heatmap of merged files
		merge1.1 <- unique(unlist(lapply(lapply(files, function(a) a$ID), as.character))) #Create a vector with the unique AGIcodes of all the lists together. 
		if(all(merge1.1 %in% as.character(tclusters[,2]))==TRUE) { #This is to prevent an error when merging if the whole list of genes is in the timeseries
		  merge1.2 <- tclusters} else {
		    merge1.2 <- rbind(tclusters, data.frame(Cluster=99, ID=merge1.1[!(merge1.1 %in% as.character(tclusters[,2]))])) #Create a 'cluster99' with the AGIcodes not in the timeseries
		  }
		fil.justID.1 <- lapply(files, function(a) a[,c(1,2)])
		merge1.2 <- merge(merge1.2, tclusters[,c(2,2)], by='ID', all.x=TRUE)
    colnames(merge1.2)[3] <- 'TimeSeries'
		cat('.')
		for(x in 1:length(fil.justID.1)) merge1.2 <- merge(merge1.2, by='ID', fil.justID.1[[x]], all.x=TRUE, sort=FALSE)
    merge1.2 <- merge1.2[order(merge1.2$Cluster),c(2,1,3:ncol(merge1.2))]
    colnames(merge1.2)[2] <- 'ALL.DATASETS'
    m.1 <- (!is.na(merge1.2[,-1]))*1
		rownames(m.1) <- merge1.2[,2]
		colnames(merge1.2)[2] <- 'AGIcode'
		merge1.2 <- cbind(merge1.2[c(1,2)],m.1)
		merge1.3 <- melt(merge1.2, id.var=c('Cluster', 'AGIcode'), value.name='val')
		colnames(merge1.3)[3] <- 'Dataset'
		merge1.3 <- merge1.3[order(merge1.3[,1]),]
		merge1.3$Cluster <- as.factor(merge1.3$Cluster)
		merge1.3$val <- as.factor(merge1.3$val)
		
		ggm <- ggplot(merge1.3, aes(Dataset, AGIcode, group=Cluster)) + 
		  geom_tile(aes(fill=val, group=Cluster)) +
		  theme(panel.background = element_blank(), 
		        axis.text.y=element_text(colour="black", size =1),
		        axis.text.x=element_text(angle=90),
		        legend.position='none')
		ggm2 <- ggm + facet_grid(Cluster~., space='free', scales='free')
		
		##	

    # Convert the timeseries in a list of genes per cluster. Trim the mylist files to keep only the TAIR codes.
    t.clus <- rbind(data.frame(Cluster=0, ID=tclusters$ID, stringsAsFactors=FALSE), tclusters)
    t.clus$Cluster <- factor(t.clus$Cluster)
    T.Series <- split(t.clus[, 2], t.clus$Cluster)
    T.Series <- lapply(T.Series, droplevels)
		names(T.Series) <- levels(t.clus$Cluster)
		files.trim <- lapply(files, function(a) a[,2])
		names(files.trim) <- lapply(files, function(a) colnames(a)[2])
		
		###Below, steps to create a matrix of percentajes of matches between genelists, to split those percentajes in quantiles to plot, and to associate a matrix of colours.
		#Only valid if more than 1 genelist is supplied by the user
		#Calculates the percentaje of matches between two genelists (%matches = number of matches / sum of the number of elements in list 1 and 2)
		somedayIdothis=FALSE
    if(somedayIdothis==TRUE) {
			func.for <- function(x) {
				mat1 <- integer()
				for(j in seq_along(files.trim)) {
					mat1 <- c(mat1, (2*length(which(x %in% files.trim[[j]])))*100/sum(length(x), length(files.trim[[j]])))
				}
				return(mat1)
			}
			#Create a matrix of percentajes of matches to plot
			mat <- as.data.frame(do.call(rbind, lapply(files.trim, func.for)))
			colnames(mat) <- names(files.trim)
			mat <- cbind(nam = rownames(mat), mat)
			mat.melt <- melt(mat, id="nam")
			#Define the quantiles of the percentajes, omitting the percentaje "100" (this is to get the best distribution of colours for the plotting, diving the data at regular intervals)
			mat.quantile.range <- quantile(mat.melt$value[mat.melt$value!=100], probs=seq(0,1,0.2))
			#Add the value of "100%" to the quantile range and round the values to the closest upper integer.
			mat.quantile.range <- ceiling(c(mat.quantile.range, 100))
			mod_mat <- matrix(findInterval(mat.melt$value, mat.quantile.range, all.inside = TRUE))
			##Create a palette of colours, omitting the value of "100%" color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(mat.quantile.range) - 2)
			#Add the "100%" as grey colour to the palette #color_palette <- c(color_palette, "black")
			#This is the definitive, manually-created palette
			color_palette <- c("grey90","#3794BF", "#9BC9DF", "#EFC29F", "#DF8640", "black")
			mat.melt <- cbind(mat.melt, mod_mat)
			#Create the labels for the guides on the graph
			label_text <- c(sapply(c(1:5), function(i) paste0(mat.quantile.range[i],"-", mat.quantile.range[i+1]," %")), "100 %")
			
			if(use.gene.levels == TRUE) {
				mat.melt$variable <- factor(mat.melt$variable, levels=gene.levels)
				mat.melt$nam <- factor(mat.melt$nam, levels=gene.levels)
			}
		}
		
		# Compute the phyper
		results.phyper <- data.frame("Cluster" = 0:(nlevels(factor(t.clus$Cluster))-1))
		for(a in seq_along(files.trim)) {results.phyper <- cbind(results.phyper, do.call(rbind, lapply(T.Series, func.phyper, b=a)))}
	  
		# Write the results
		#Extract the columns with the phyper probabilities
		sequence <- seq.int(4, length.out=len.file, by=4)
		results.phyper.probsonly <- results.phyper[,c(1,sequence)]
		thr <- results.phyper.probsonly[,-1] < p.threshold

		#If stringent==TRUE only the clusters that have a phyper probability under the threshold for all the genelists (one or more than one genelist) are selected. If FALSE, all clusters for which at least a genelist is under the threshold are selected.
		if(length(files)>1) {
      if(stringent==TRUE) files.out0 <- results.phyper.probsonly[rowSums(thr)==ncol(thr),] else files.out0 <- results.phyper.probsonly[rowSums(thr)>=1,]   #The first part counts the number or TRUE's. If all columns are TRUE, returns TRUE.
		} else files.out0 <- results.phyper.probsonly[thr,]
         
    thresholed.clusters <- files.out0[,"Cluster"]
      
		mergefile1 <- droplevels(mergefile1[mergefile1$Cluster %in% thresholed.clusters,])
		mergefile2 <- droplevels(mergefile2[mergefile2$Cluster %in% thresholed.clusters,])
    mergefile2 <- droplevels(mergefile2[order(mergefile2$Cluster, mergefile2$ID),])
		results.phyper <- droplevels(results.phyper[results.phyper$Cluster %in% thresholed.clusters,])
		results.phyper.probsonly <- droplevels(results.phyper.probsonly[results.phyper.probsonly$Cluster %in% thresholed.clusters,])

    res.heatmap <- results.phyper.probsonly
    rownames(res.heatmap) <- res.heatmap$Cluster
    res.heatmap <- as.matrix(res.heatmap[,-1])
    #colnames(res.heatmap) <- gsub('prob\\.T\\.', '', colnames(res.heatmap))
       
		if(plot.phyper==TRUE) {
				write.table(mergefile1,file="RES_mergedAGIs.txt", row.names=FALSE, sep="\t")
				write.table(mergefile2,file="RES_mergedAGIs_listsonly.txt", row.names=FALSE, sep="\t")
				write.table(results.phyper, file="RES_phyper_allresults.txt",row.names=FALSE,sep="\t")
				write.table(results.phyper.probsonly, file="RES_phyper_results_probabilitiesonly.txt",row.names=FALSE,sep="\t")
				cat("Analysis finished. The files RES_merge.txt, RES_results.txt and RES_results_probsonly.txt have been created in the working directory.\n")
		}
		#Function to slice the phyper probabilities from the dataframe
		func.sliceprob <- function(vseq, a, b) {
			res <- data.frame("Cluster"=results.phyper$Cluster, results.phyper[,vseq])
			levels.unique <- unique(as.character(res$Cluster))
			res$Cluster <- factor(res$Cluster, levels=levels.unique)
			colnames(res) <- c("Cluster", names(files.trim), paste0(names(files.trim), '.proportion'))
			res.phy <- melt(res[, c(1:(length(vseq)/2+1))], id="Cluster", value.name=a)
      res.pro <- melt(res[, c(1,(length(vseq)/2+2):(2*length(vseq)/2+1))], id="Cluster", value.name=b)
      res <- data.frame(res.phy, round(res.pro[,3],1))
      colnames(res) <- c('Cluster', 'variable', a, b)
			res <- data.frame(res, Prob=(factor("grey70", levels=c("grey40", "grey70", "black"))))
			res[res[,3]<0.001, 5] <- "grey40"
			res[res[,3]<0.0001, 5] <- "black"
			return(res)
		}
		
		#Slice the probability values
		res <- func.sliceprob(c(sequence, sequence+1), "phyper", "overlap")
		res <- cbind(res, vlog10=-log10(res$phyper))
		
		if(use.gene.levels == TRUE) res$variable <- factor(res$variable, levels=gene.levels)
		
    s <- length(thresholed.clusters)*0.6
    if(s<8) s<-8
    
		if(plot.heatmap==TRUE) {
		  cat('Printing heatmap...\n')
		  p.color<- c('grey80', 'yellow', 'orange', "red")
		  b.r <- c(1, 1e-5, 1e-10, 1e-15, 1e-100)
		  ifelse(ncol(res.heatmap)<=10, d.n <- TRUE, d.n <- FALSE)
		  
      if(ncol(res.heatmap)>30) wd1<-2+ncol(res.heatmap)/3 else wd1<- 2+ncol(res.heatmap)/4
      
		  pdf(heatmap.filename, width=wd1)
		  if(dim(res.heatmap)[2]!=1) pheatmap(-log10(res.heatmap + 1e-100), treeheight_row=15, treeheight_col=15, cellwidth=20, cluster_rows=FALSE, color=p.color, breaks= -log10(b.r), legend_breaks=-log10(c(1e-5, 1e-10, 1e-15)), legend_labels=as.character(c(1e-5, 1e-10, 1e-15)), border_color='white',display_numbers=d.n,number_format='%.1f')
      pheatmap(m.1, main='Similarity of datasets (genes in red)', border_color=NA,treeheight_row=15, treeheight_col=15, show_rownames=FALSE, legend=FALSE)
      dev.off()
      if(plot.heatmap.pergene==TRUE) ggsave('heatmap_pergene.pdf', ggm2, dpi=1200, width=(2+0.1*length(files)), height=nrow(merge1.2)*0.02, limitsize=FALSE)
		  }
    
		gg.theme <- theme(panel.background = element_blank(),
			panel.grid.major = element_line(colour="#CCCCCC", linetype="dashed"),
			#panel.border = element_rect(fill = NA, colour="black"),
			axis.line = element_line(colour="black"),
			axis.text = element_text(colour="black", size =s),
			axis.title = element_text(colour="black", size=s),
			strip.text = element_text(size=s, face="bold"))
		
		#Phyper plot per gene list.
	  if(plot.phyper==TRUE) {
				gg<- ggplot(data=res, aes(x=Cluster, y=vlog10)) +
					geom_bar(aes(fill=Prob),stat="identity") +
					scale_fill_identity() +
					#geom_hline(aes(yintercept=-log10(0.001)), linetype="dashed") +
					geom_hline(aes(yintercept=-log10(0.00001)), linetype="dashed") +
					gg.theme +
					ylab("-log10(p-value)") +
				  xlab("Cluster\n (dark grey bars, P<0.001, black bars, P<0.0001;\n dashed line, P=0.00001)\n On the bars, proportion(%) of overlapping genes") +
					facet_grid(variable ~ .)
			  gg <-	gg + geom_text(data=res,aes(y=(vlog10+(max(vlog10[is.finite(vlog10)]))*0.06), x=Cluster, label=overlap))
				cat('Printing phyperplot.pdf...\n')				
				ggsave("phyperplot.pdf", gg, width=0.7*(3+ length(thresholed.clusters)), height=1.5+(len.file*3.3), limitsize=FALSE)
				
				#Plot of the percentaje of genes in each gene list that overlap with the genes in a cluster, with respect to the total number of genes in that list
				gg2 <- ggplot(data=res, aes(x=Cluster, y=overlap)) +
					geom_bar(stat="identity") +
					gg.theme +
					ylab("Overlapping genes per cluster (%)") +
          xlab("Cluster") +
					facet_grid(variable ~ .)
				if(plot.proportions==TRUE){
          cat('Printing phyperplot_proportions.pdf...\n')				
  				ggsave("phyperplot_proportions.pdf", gg2, width=0.7*(3+ length(thresholed.clusters)), height=len.file*3, limitsize=FALSE)
				}
				if(use.gene.levels == TRUE) {
					#Plot for the relationships between the different lists of genes, representing the percentaje of genes that each gene list has with any other list.
					gg3 <- ggplot(data=mat.melt, aes(x=variable, y=nam, fill=factor(mod_mat))) +
						geom_tile() +
						geom_tile(colour="white", show_guide=FALSE) +
						#geom_tile(color="white") +
						#scale_fill_gradient2(low = "white", high = "steelblue") +
						scale_fill_manual(values = color_palette, name="", labels = label_text, guide = guide_legend(reverse=TRUE)) +
						ylab("Gene list") + xlab("Gene list") +
						theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_blank())
					
					ggsave("Percentage_of_Genes_in_common.pdf", gg3, width=nlevels(factor(mat.melt$variable))/2, height=length(levels(mat.melt$variable))/2, limitsize=FALSE)
				}
	  }
	files.out <- list(thresholed.clusters, files.out)
	return(files.out)
}

func.phyperpair <- function(a) {
  func.phyperpairwise <- function(aa) {
    a1 <- aa[[1]]
    a2 <- aa[[2]]
    len.a <- length(a1)
    k <- length(a2)
    overlaps <- sum(a1 %in% a2)
    n <- 24000-len.a
    prob <- phyper(overlaps, len.a, n, k, lower.tail=FALSE)
    results <- data.frame(n.genes1=len.a, n.genes2=k, overlaps, prob, proportion=(100*2*overlaps/(len.a+k)))
    }
  cat('.')
  ##Pairwise comparison of the lists of genes.
  f.pair <- lapply(a, function(a) unique(a[,2]))
  #f.pair <- do.call(rbind, combn(f.pair, 2, func.phyperpairwise, simplify=FALSE))
  #f.pair <- cbind(do.call(rbind, combn(names(files), 2, simplify=FALSE, function(a) data.frame(var1=a[1], var2=a[2]))), f.pair)
  #f.pair.prob <- f.pair[, c(1,2,5)]
  f <- lapply(f.pair, function(a) lapply(f.pair, function(b) list(a,b)))
  f1 <- lapply(f, function(a) lapply(a, func.phyperpairwise))
  f2 <- do.call(rbind, do.call(c, f1)) #It seems that if I do do.call(rbind, do.call(rbind, f1)), the rbind takes the first element of each of the 3 root lists, then the second, then the third, etc (l1e1, l2e1, l3e1, l1e2, l1e2...etc). Doing 'rbind, c' instead, puts the list in a dataframe in the same original order (l1e1, l1e2, l1e3, l2e1, l2e2... etc.)
  d <- data.frame(ID1=rep(names(a), each=length(names(a))), ID2=rep(names(a), length(names(a))), f2)
  f.prob <- acast(d, ID1~ID2, value.var="prob")
  f.proportion <- acast(d, ID1~ID2, value.var="proportion")
  
  fil <- lapply(lapply(a, function(a) a$ID), as.character)
  merge1 <- unique(unlist(fil))
  merge1.justID <- merge1 <- data.frame(ID=merge1, TAIR=merge1)
  fil.justID <- lapply(a, function(a) a[,c(1,2)])
  
  
  cat('.')
  for(x in 1:length(a)) merge1 <- merge(merge1, a[[x]], by="ID", all.x=TRUE, sort=FALSE)
  for(x in 1:length(fil.justID)) merge1.justID <- merge(merge1.justID, by='ID', fil.justID[[x]], all.x=TRUE, sort=FALSE)
  m <- (!is.na(merge1.justID[,-c(1,2)]))*1
  rownames(m) <- merge1.justID[,1]
  
  ##Calculation of all the overlaps between genelists (not just pairwise)
          #'merge1.justID' is the result of merging the lists (a file with a column per genelist, and "NA's" for genes that are not present in one or more of the lists)
          # the 'apply' function below creates a new column with letters that identify the genelists that have a NA.
            c.numbers <- c(1:100)  
            merge1.justID2 <- merge1.justID[,-c(1,2)]        
            table1 <- apply(merge1.justID2,MARGIN=1, FUN=function(a) paste0('ID', paste0(c.numbers[which(!is.na(a))], collapse='.')))
            capture.output(file='alloverlaps_openwithWordpad.txt',append=T,
               cat('\nThe overlaps of all the combinations of genelists are shown below. Check the letters (a,b,c... etc) for correspondence\n\n\n'),
               cat('\nTotal number of genes per genelist:\n'),
               colSums(!is.na(merge1.justID2)),
               cat('\nGenelist name equivalences and number of overlapping genes:\n'),
               data.frame('genelist'=colnames(merge1.justID2), 'ID'= c.numbers[1:ncol(merge1.justID2)]),
               table('\nOverlaps:'=table1)
            )
  
  
  ful <- fullmap[,c(1,6,7)]
  ful <- ful[!duplicated(ful[,1]),]
  merge0 <- merge(merge1, ful, by='TAIR', all.x=TRUE, sort=FALSE)
  write.table(merge0, 'pairwise_merge.txt', sep='\t', row.names=FALSE)
  write.table(d, 'pairwise_phyper.txt', sep='\t', row.names=FALSE)
    
  cat('.')
  if(length(names(a)>30)){
    wd2<- 4+length(names(a))/2
    hd2<- 3+length(names(a))/2} else{
    wd1<- 4+length(names(a))/3
    hd2<- 3+length(names(a))/3}
    
  p.color<- c('grey80', 'yellow', 'orange', "red")
  b.r <- c(1, 1e-5, 1e-10, 1e-15, 1e-100)
  
  pdf('heatmap_phyperpairwise.pdf', width=wd2, height=hd2)
  pheatmap(-log10(f.prob+1e-100), treeheight_row=15, treeheight_col=15, cellwidth=30, cellheight=30,  color=p.color, breaks= -log10(b.r), legend_breaks=-log10(c(1e-5, 1e-10, 1e-15)), display_numbers=TRUE, number_format='%.1f', main='-log10(probabilities)')
  pheatmap(f.proportion, cellwidth=30, treeheight_row=15, treeheight_col=15, cellheight=30, display_numbers=TRUE, number_format='%.1f', main='overlaps(%)')
  pheatmap(m, main='Similarity of datasets (genes in red)', border_color=NA,treeheight_row=15, treeheight_col=15, show_rownames=FALSE, legend=FALSE)
  dev.off()
  
  #Put all the AGI codes of the list in a single file, to use as reference for merging.
  fil <- lapply(lapply(a, function(a) a$ID), as.character)
  merge1 <- unique(unlist(fil))
  merge1.justID <- merge1 <- data.frame(ID=merge1, TAIR=merge1)
  fil.justID <- lapply(a, function(a) a[,c(1,2)])
  
  
  cat('.')
  for(x in 1:length(a)) merge1 <- merge(merge1, a[[x]], by="ID", all.x=TRUE, sort=FALSE)
  for(x in 1:length(fil.justID)) merge1.justID <- merge(merge1.justID, by='ID', fil.justID[[x]], all.x=TRUE, sort=FALSE)
  m <- is.na(merge1.justID[,-c(1,2)])*1
  rownames(m) <- merge1.justID[,1]
  
  ful <- fullmap[,c(1,6,7)]
  ful <- ful[!duplicated(ful[,1]),]
  merge0 <- merge(merge1, ful, by='TAIR', all.x=TRUE, sort=FALSE)
  write.table(merge0, 'pairwise_merge.txt', sep='\t', row.names=FALSE)
  write.table(d, 'pairwise_phyper.txt', sep='\t', row.names=FALSE)
  
  cat('\n...Plotting heatmaps of pairwise probabilities and overlaps...\n')
}

func.time <-function(){
  #This function reads the HighLightCombo and converts it into a TimeSeries format for the Interrogate.timeseries function
  HL <- read.delim('HighLight-Combo3.txt')
  library(reshape2)
  HLc2 <- melt(HL, id.vars=c(1,2), variable.name='var', value.name='val')
  s.split <- do.call(rbind,strsplit(as.character(HLc2$var), split='_', fixed=TRUE))
  HLc3 <- cbind(s.split[,1:3], HLc2)
  HLc3[,1] <- gsub('T', '', HLc3[,1])
  HLc3 <- HLc3[, c(5,4,7,2,3,1)]
  colnames(HLc3) <- c('ID','Probe', 'value','va.treat', 'va.replicate', 'Time')
  #fullmap <- read.delim('TAIR_fullmapping_NCBI_GO_symbol.txt')
  HLc4 <- merge(HLc3, unique(fullmap[,c(1,3)]), by.x='ID', by.y='TAIR', sort=FALSE, all.x=TRUE)
  HLc4 <- merge(HLc4, tclusters, by='ID', all.x=TRUE, sort=FALSE)
  HLc4 <- HLc4[, c(1,8,7,2,3:6)]
  HLc4[is.na(HLc4$Cluster),'Cluster'] <- 99
  HLc4$Cluster <- as.factor(HLc4$Cluster)
  HLc4$Time <- factor(HLc4$Time, levels=c(0:12), labels=seq(0,6,0.5))
  HLc4$Time <- as.numeric(HLc4$Time)
  return(HLc4)
}