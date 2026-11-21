# Include the whole code for the app.R
library(uwot)
library(NMF)
TBnanostring <- readRDS("TBnanostring.rds")

# PCA
pca <- prcomp(TBnanostring[, -1])
# UMAP
set.seed(1123)
um <- umap(TBnanostring[, -1], n_neighbors = 15)
# NMF
tb_ng <- sweep(TBnanostring[, -1], 2, apply(TBnanostring[, -1], 2, min))
set.seed(32)
nmf_tb <- nmf(tb_ng, 2, nrun = 20)
W <- basis(nmf_tb)
colnames(W) <- c("NMF1", "NMF2")


library(shiny)
ui <- fluidPage(
  titlePanel("TB NanoString K-means Clustering"),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Dimension Reduction",
                  c("PCA", "UMAP", "NMF")),
      numericInput("clusters", "Cluster count", 2,
                   min = 1, max = 9)
    ),
    
    
    mainPanel(
      plotOutput("plot1"),
      tableOutput("table1"),
      # Create a preformatted text output area for ARI
      verbatimTextOutput("ari")
    )
  )
)

# Define the server function for calculations and outputs
server <- function(input, output, session) {
  selectedData <- reactive({
    switch(input$method,
           PCA = pca$x[, 1:2],  
           UMAP = um,          
           NMF = W)           
  })
  
  # Perform k-means and store the clustering results
  # Recalculate only when the selected method or k changes
  clusters <- reactive({
    set.seed(10)
    kmeans(selectedData(), input$clusters, nstart = 25)
  })
  
  
  output$plot1 <- renderPlot({
    palette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))
    
    # Set plot margins: bottom, left, top, right
    par(mar = c(5.1, 4.1, 0, 1))
    plot(selectedData(),
         # Color each point by its k-means cluster
         col = clusters()$cluster,
         # Assign different shapes based on TB status
         # TB = circle (16), LTBI = triangle (17)
         pch = ifelse(TBnanostring$TB_Status == "TB", 16, 17),
         # Set point size
         cex = 1.2)
    
    # Add cluster centers to the existing plot
    # pch = 4: cross shape
    # cex = 4: symbol size
    # lwd = 4: line width
    points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
    
    # Add a legend to explain the shapes
    legend("topright",
           legend = c("TB", "LTBI"),
           pch = c(16, 17))
  })
  
  # Generate the cross-tabulation table
  output$table1 <- renderTable({
    as.data.frame.matrix(
      table(cluster = clusters()$cluster,
            TB_Status = TBnanostring$TB_Status)
    )
    
    # Add row names as cluster numbers
  })
  
  # Display the ARI
  output$ari <- renderPrint({
    # Measure agreement between clusters and actual TB status
    mclust::adjustedRandIndex(clusters()$cluster,
                              TBnanostring$TB_Status)
  })
}

# Launch the Shiny app using the defined UI and server
shinyApp(ui = ui, server = server)