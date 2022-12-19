library(shiny)
source("ui_load.R", local=TRUE)

ui = fluidPage(
  # website title and icon
  tags$head(tags$link(rel="shortcut icon",
                      href=knitr::image_uri("PES.ico"),
                      type="image/x-icon")),
  titlePanel(div("Photoelectron Spectroscopy",
                 style="padding-bottom:0.8rem;"),
             windowTitle="Photoelectron Spectroscopy"),
  
  # sidebaar panel with selection and description
  sidebarPanel(
    selectInput("atom", "Atom:",
                choices = choice_input_list(),
                selected = 1),
    # would like to reverse this to match the graph
    sliderInput("Range", "Range:", min=0, max=2, value=c(0, 2), step=0.1),
    htmlOutput("avee"),
    
    # description of PES
    HTML("<p>Photoelectron spectroscopy refers to the energy measurement of 
         electrons emitted from solids, gases, or liquids by the photoelectric 
         effect, to determine the binding energies of electrons in a substance. 
         The term refers to various techniques, depending on whether the 
         ionization energy is provided by an X-ray photon, an extreme 
         ultraviolet (EUV) photon, or an ultraviolet (UV) photon. Regardless of 
         the incident photon beam, however, all photoelectron spectroscopy 
         revolves around the general theme of surface analysis by measuring the 
         ejected electrons (Wikipedia).<p>"),
    
    # taken from https://css-tricks.com/examples/hrs/
    HTML("<hr style='border:0; height:1px; background:#333;
         background-image:linear-gradient(to right, #ccc, #999, #ccc);'>"),
    checkboxInput("use_gaus_peaks", "Use gaussian-like peaks (experimental)"),
    uiOutput("std_render")
  ),
  
  # output plots and table
  mainPanel(
    align = "center",
    plotOutput("plotxy"),
    tableOutput("detail")
  )
)
