#!/bin/R
#Tychele N. Turner
#Laboratory of Aravinda Chakravarti, Ph.D.
#Johns Hopkins University School of Medicine
#Protein Plotting Script for Shiny
#Programming Language: R
#Updated 06/15/2013

#Description: This script is the ui.R script for the Plot Protein Shiny Application

library(shiny)

# Define UI 
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Plot Protein: Visualization of Mutations"),

  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the br()
  # element to introduce extra vertical spacing
  sidebarPanel(
  	 
	fileInput(inputId="mutationFile", label="Mutation File:"),
				 
	fileInput(inputId="proteinArchitectureFile", label="Protein Architecture File:"),

	fileInput(inputId="postTranslationalModificationFile", label="Post Translational Modification File:"),
				 
	textInput("proteinLength", "Length of Protein:", "Data"),
	
	textInput("nameOfQuery", "Name of Query:", "Data"),
	
	textInput("tickSize", "Tick Size:", 10),	

	checkboxInput(inputId = "labels",
      label = strong("Show Labels"),
      value = FALSE),

	checkboxInput(inputId = "zoom",
      label = strong("Zoom In?"),
      value = FALSE),

	textInput("zoomStart", "Zoom Start:", "Data"),	

	textInput("zoomEnd", "Zoom End:", "Data")
		 
  ),

  # Show a tabset that includes a plot and table
  mainPanel(
	h5("Tychele N. Turner, tycheleturner@gmail.com"),
	h5("Laboratory of Aravinda Chakravarti, Ph.D."),
	
    tabsetPanel(
      tabPanel("Plot", plotOutput("plot")), 
      tabPanel("Table", tableOutput("table"))
    )
  )
))
