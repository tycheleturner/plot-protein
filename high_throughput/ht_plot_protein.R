#!/bin/R
#Tychele N. Turner
#Protein Plotting Script for Snakefile
#Programming Language: R
#Updated 06/09/2016

#Description: This script takes mutation information at the protein level and plots out the mutation above the schematic of the protein. It also plots the domains. 

#NOTE: All files should be referring to the same isoform of the protein. This is imperative for drawing the plot correctly.

#Usage:
##Basic
###Rscript plotProtein.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -p psen1_post_translation_file.txt -l 463

##Advanced
###Rscript plotProtein.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -p psen1_post_translation_file.txt -l 464 -n Disease -t 25 -s yes -z yes -b 50 -c 100

#will have to run the install for optparse if never installed before
#install.packages("optparse")

library("optparse")

option_list <- list(
    make_option(c('-m', '--mutations'), action='store', type='character', default='mutationFile.txt', help='This is the mutation file. It should be a tab-delimited file containing 5 columns (ProteinId, GeneName, ProteinPositionOfMutation, ReferenceAminoAcid, AlternateAminoAcid) NO HEADER FOR NEEDED FOR THIS FILE. (REQUIRED)'),
    make_option(c('-a', '--architecture'), action='store', type='character', default='architectureFile.txt', help='This is the protein architecture file. It should be a tab-delimited file containing 3 columns (architecture_name, start_site, end_site). This file NEEDS the header and it is the same as what was previously written. This information can be downloaded from the HPRD (http://hprd.org/). Although the most recent files are quite old so looking in the web browser you can get much more up to date information. (REQUIRED)'),
    make_option(c('-p', '--posttranslational'), action='store', type='character', default='posttranslationalFile.txt', help='This is the protein post-translational modification file. This is a tab-delimited file with only one column and that is the site. This file NEEDS a header and is as previously written (site). (REQUIRED)'),
    make_option(c('-l', '--length'), action='store', type='character', default=100, help='protein length file (REQUIRED)'),
    make_option(c('-n', '--name'), action='store', type='character', default='Test', help='Name of your query. Default is Test'),
    make_option(c('-t', '--ticksize'), action='store', type='numeric', default=10, help='Size of ticks on x-axis. Default is 10'),
    make_option(c('-s', '--showlabels'), action='store', type='character', default='no', help='Option to show labels. Default is no'),
    make_option(c('-z', '--zoom'), action='store', type='character', default='no', help='Option to zoom in. Default is no'),
    make_option(c('-b', '--zoomstart'), action='store', type='numeric', default=1, help='Starting number for zoom in. Use if zoom option is set to yes. Default is 1'),
    make_option(c('-c', '--zoomend'), action='store', type='numeric', default=10, help='Ending number for zoom in. Use if zoom option is set to yes. Default is 10')
)
opt <- parse_args(OptionParser(option_list = option_list))

set = opt$set #global quality threshold (can vary)
mutationFile <- opt$mutations
proteinArchitectureFile <- opt$architecture
postTranslationalModificationFile <- opt$posttranslational
proteinLength <- opt$length
nameOfYourQuery <- opt$name
tickSize <- opt$ticksize
showLabels <- opt$showlabels
zoomIn <- opt$zoom
if(zoomIn == "yes"){
	zoomStart <- opt$zoomstart
	zoomEnd <- opt$zoomend
}

####################ANALYSIS####################
#Read in the files
var <- read.table(mutationFile, sep="\t")
pa <- read.table(proteinArchitectureFile, sep="\t", header=FALSE)
pt <- read.table(postTranslationalModificationFile, sep="\t", header=FALSE)
plength <- read.table(proteinLength, sep="\t", header=F)
proteinLengthPlot <- as.numeric(plength$V1[1])

############PLOTTING#############
#x is the input data, y is rpt, z is rpa from HPRD
pdf(paste(as.character(var[1,2]), "_", as.character(var[1,1]), "_protein_plot.pdf", sep=""), height=7.5, width=10)
#par(oma=c(8, 1.2, 8, 1.2))
layout(matrix(c(1,2),nrow=1), widths=c(1,3))
par(oma=c(4, 0, 4, 0), mar=c(5, 0, 4, 0) + 0.4)

#stable legend
plot((-30:-15), rep(-1, 16), col="white", type="l", ann=FALSE, bty="n", xaxt="n", yaxt="n", xlim=c(-160, -15), ylim=c(1,-5.5))
	
#query text
text(-100,-2.5,nameOfYourQuery, col="blue", cex=0.9, font=2)

xlimRegion <- c(0, proteinLengthPlot)
	if(zoomIn == "yes") {
          xlimRegion <- c(as.numeric(zoomStart), as.numeric(zoomEnd))
	}




xlimRegion <- c(0, as.numeric(proteinLengthPlot))
	if(zoomIn == "yes") {
          xlimRegion <- c(as.numeric(zoomStart), as.numeric(zoomEnd))
	}
	
plot((1:as.numeric(proteinLengthPlot)), rep(-2, as.numeric(proteinLengthPlot)), type="l", lwd=5, main=paste("Amino Acid Changes in", " ", as.character(var[1,2]), " ", "(", as.character(var[1,1]), ")", sep=""), xlab="Amino Acid Position", ylab="", ylim=c(-1,-4), cex.lab=0.9, cex.main=1, yaxt="n", xlim=xlimRegion, xaxt="n", ann=FALSE, bty="n")

#Plot mutations
points(var[,3], rep(-2.5, length(var[,3])), pch=19, col=rgb(0,0,1,1/2), cex=0.7)

if(showLabels == "yes"){
	#Label mutations
	for(i in 1:nrow(var)){
		text(var[i,3], rep(-2.7, length(var[i,3])), as.character(var[i,3]), col="blue", cex=0.9, srt=90, adj = 0)
	}
}

ticks=seq(0,as.numeric(proteinLengthPlot), by=tickSize) 
axis(side = 1, at = ticks, las=3)

#labels
for(i in 1:length(pt$V1)){
	segments(as.numeric(pt$V1[i]), -2, as.numeric(pt$V1[i]), -2.25, lwd=2, col="black")
	points(as.numeric(pt$V1[i]), -2.25, pch=19, col="deeppink", cex=0.7)
}
for(i in 1:length(pa$V2)){
	rect(as.numeric(pa$V2[i]), -2.05, as.numeric(pa$V3[i]), -1.95, col="lightseagreen")
}
for(i in 1:length(pa$V1)){
	text(median(c(as.numeric(pa$V2[i]), as.numeric(pa$V3[i]))), -1.80, pa$V1[i], cex=1)
}
legend("topright", c("Protein Domain", "Post-Translational Modification"), fill=c("lightseagreen", "deeppink"),  box.col="white", bg="white", cex=1)
dev.off()


