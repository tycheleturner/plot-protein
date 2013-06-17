#!/bin/R
#Tychele N. Turner
#Laboratory of Aravinda Chakravarti, Ph.D.
#Johns Hopkins University School of Medicine
#Protein Plotting Script for Shiny
#Programming Language: R
#Updated 06/15/2013

#Description: This script is the server.R script for the Plot Protein Shiny Application

library(shiny)

# Define server 
shinyServer(function(input, output) {

  output$plot <- renderPlot({
    
########mutation file########
	if (is.null(input$mutationFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    var <- read.table(input$mutationFile$datapath, sep="\t")
########mutation file########
		
########protein architecture file########
	if (is.null(input$proteinArchitectureFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    pa <- read.table(input$proteinArchitectureFile$datapath, sep="\t", header=TRUE)
########protein architecture file########
	
########post translational modification file########
	if (is.null(input$postTranslationalModificationFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    pt <- read.table(input$postTranslationalModificationFile$datapath, sep="\t", header=TRUE)
########post translational modification file########
	
	############PLOTTING#############
	#x is the input data, y is rpt, z is rpa from HPRD
	#pdf(paste(as.character(var[1,2]), "_protein_plot.pdf", sep=""), height=7.5, width=10)
	#par(oma=c(8, 1.2, 8, 1.2))

	xlimRegion <- c(-145, as.numeric(input$proteinLength))

	if(input$zoom) {

		xlimRegion <- c(as.numeric(input$zoomStart), as.numeric(input$zoomEnd))
}

	plot((1:as.numeric(input$proteinLength)), rep(-2, as.numeric(input$proteinLength)), type="l", lwd=5, main=paste("Amino Acid Changes in", " ", as.character(var[1,2]), " ", "(", as.character(var[1,1]), ")", sep=""), xlab="Amino Acid Position", ylab="", xlim=xlimRegion, ylim=c(-1,-4), cex.lab=0.9, cex.main=1, yaxt="n", xaxt="n")
	
	ticks=seq(0,as.numeric(input$proteinLength), by=as.numeric(input$tickSize))
	#Specify the ticks and grids you want	
	axis(side = 1, at = ticks, las=3)


	#Plot mutations
	points(var[,3], rep(-2.5, length(var[,3])), pch=19, col="blue", cex=0.7)

	#Label mutations
	if (input$labels) {
	for(i in 1:nrow(var)){
		text(var[i,3], rep(-2.7, length(var[i,3])), paste(as.character(var[i,4]), as.character(var[i,3]), as.character(var[i,5]), sep=""), col="blue", cex=0.9, srt=90, adj = 0)
	}
}

	#labels
	text(-100,-2.5,input$nameOfQuery, col="blue", cex=0.9)
	for(i in 1:length(pt$site)){
		segments(as.numeric(pt$site[i]), -2, as.numeric(pt$site[i]), -2.25, lwd=2, col="black")
		points(as.numeric(pt$site[i]), -2.25, pch=19, col="deeppink", cex=0.7)
	}

	for(i in 1:length(pa$start_site)){
#		lines(c(as.numeric(pa$start_site[i]):as.numeric(pa$end_site[i])), rep(-2,(as.numeric(pa$end_site[i])- as.numeric(pa$start_site[i]))+1),lwd=10, col=c("purple"))
		rect(as.numeric(pa$start_site[i]), -2.05, as.numeric(pa$end_site[i]), -1.95, col="lightseagreen")
	}
	for(i in 1:length(pa$architecture_name)){
		text(median(c(as.numeric(pa$start_site[i]), as.numeric(pa$end_site[i]))), -1.80, pa$architecture_name[i], cex=1)
	}
	legend("topright", c("Protein Domain", "Post-Translational Modification"), fill=c("lightseagreen", "deeppink"), cex=1)
	#dev.off()


    })
  #end of plot
  
  
  #called when user uploads a file
   output$table <- renderTable({
########mutation file########
	if (is.null(input$mutationFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    var <- read.table(input$mutationFile$datapath, sep="\t")
########mutation file########
		
########protein architecture file########
	if (is.null(input$proteinArchitectureFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    pa <- read.table(input$proteinArchitectureFile$datapath, sep="\t", header=TRUE)
########protein architecture file########
	
########post translational modification file########
	if (is.null(input$postTranslationalModificationFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    pt <- read.table(input$postTranslationalModificationFile$datapath, sep="\t", header=TRUE)
########post translational modification file########
	
	
	#Architecture
	rpa <- pa[order(pa$start_site),]
	for(i in 1:nrow(rpa)){
		rpa$domain_name[i] <- paste(rpa$architecture_name[i], " (", "Domain Number",  i, ") ", sep="")
	}

	var$name <- "Not in Domain"
	for(i in 1:nrow(var)){
		for(j in 1:nrow(rpa)){
			if(as.numeric(var[,3][i]) >= as.numeric(rpa$start_site[j]) & (as.numeric(var[,3][i]) <= as.numeric(rpa$end_site[j])))
			var$name[i] <- rpa$domain_name[j]
		}
	}

	colnames(var) <- c("proteinID", "geneName", "aminoAcidPosition", "refAA", "altAA", "domain")

	#Post Translational Modification
	rpt <- pt[order(pt$site),]

	var$postTranslationSite <- "No"

	for(i in 1:nrow(var)){
		for(j in 1:length(rpt)){
			if(as.numeric(var[,3][i]) == as.numeric(rpt[j]))
			var$postTranslationSite[i] <- "Yes"
		}
	}

	data.frame(var)
	  
  }) # end of filetable 

 })

