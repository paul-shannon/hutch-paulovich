library(shiny)
library(DT)
library(shinyBS)
library(yaml)
library(later)
#----------------------------------------------------------------------------------------------------
printf = function (...) print (noquote (sprintf (...)))
#----------------------------------------------------------------------------------------------------
tbl <- get(load("tbl.numa.RData"))
#----------------------------------------------------------------------------------------------------
# https://docs.google.com/document/d/1Qwj9-vj8Q7b0GWLCs5UmHirQqGkejMex5w11E-CXLvU/edit?usp=sharing
#----------------------------------------------------------------------------------------------------
mainTableTab <- function()
{
   destinations <- list(" " = "goNowhere",
                        "All Tabs" = "allTabs",
                        "GeneCards" = "genecards",
                        "HomoloGene" = "homologene",
                        "PubMed (curated)" = "pubmedCurated",
                        "PubMed (gene + longevity)" = "pubmedLongevity",
                        "Orthologs" = "orthologs"
                        )
   treatmentOptions <- c("none", "ATMi", "DMSO")
   treatmentSelector <- checkboxGroupInput("treatmentCheckBoxGroup",
                                           label = h3("Treatments"),
                                           choices = treatmentOptions,
                                           selected=treatmentOptions)
   radiationOptions <-  c("mock", "1Gy", "2Gy", "5Gy", "10Gy", "mock")
   radiationSelector <- checkboxGroupInput("radiationCheckBoxGroup",
                                           label = h3("Radiation"),
                                           choices = radiationOptions,
                                           selected=radiationOptions)

   tab <- tabPanel(title="Gene Table", value="byGeneTab",
             sidebarLayout(
                 sidebarPanel(
                     treatmentSelector,
                     radiationSelector,
                     width=3),
                 mainPanel(
                    DTOutput("mainTable"),
                    style="margin-top:5px;"
                    ) # mainPanel
               ) # sideBarLayout
             ) # tabPanel

   tab

} # mainGeneTab
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

   tags$head(tags$script(type="module", src="https://unpkg.com/x-frame-bypass")),
   titlePanel("NUMA Proteomics Demo"),
   tabsetPanel(type="tabs", id="lcSnpsTabs",
               mainTableTab()
               )
) # ui
#----------------------------------------------------------------------------------------------------
server <- function(session, input, output) {

   #observeEvent(input$treatmentCheckBoxGroup, ignoreInit=TRUE, {
   #   newValue <- input$treatmentCheckBoxGroup
   #   printf("newValue: %s", paste(newValue, collapse=" "))
   #   })

   output$mainTable <- DT::renderDataTable({
       treatmentsToInclude <- input$treatmentCheckBoxGroup
       tbl <- subset(tbl, treatment %in% treatmentsToInclude)
       radiationLevelsToInclude <- input$radiationCheckBoxGroup
       tbl <- subset(tbl, radiation %in% radiationLevelsToInclude)
       DT::datatable(tbl,
                     rownames=FALSE,
                     options=list(dom='<lfi<t>>',
                                  scrollX=TRUE,
                                  autoWidth=TRUE,
                                  columnDefs=list(list(width="10%", targets=c(0,1)),
                                                  list(width="40%", targets=c(2,3)),
                                                  list(width="10%", targets=4)),
                                  #lengthMenu = c(3,5,10,50),
                                  #pageLength = 5,
                                  paging=FALSE),
                     selection="single")

      })

} # server
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
runApp(shinyApp(ui=ui, server=server), port=9004)
#shinyApp(ui=ui, server=server)
