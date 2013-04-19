library(shiny)

# Define UI for random distribution application 
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

	checkboxInput(inputId = "labels",
      label = strong("Show Labels"),
      value = FALSE),

	checkboxInput(inputId = "zoom",
      label = strong("Zoom In?"),
      value = FALSE),

	textInput("zoomStart", "Zoom Start:", "Data"),	

	textInput("zoomEnd", "Zoom End:", "Data")
		 
  ),

  # Show a tabset that includes a plot, summary, and table view
  # of the generated distribution
  mainPanel(
	h5("Tychele N. Turner, tycheleturner@gmail.com"),
	h5("Laboratory of Aravinda Chakravarti, Ph.D."),
	
    tabsetPanel(
      tabPanel("Plot", plotOutput("plot")), 
      tabPanel("Table", tableOutput("table"))
    )
  )
))
