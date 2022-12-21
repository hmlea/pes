source("PES_functions.R", local=TRUE)

shinyServer(function(input, output, session){
  # observe the change in the slider
  observe({
    atomchoice = as.numeric(input$atom)
    maxlimit = ceiling(max(na.omit(binding_energies[be_rows, atomchoice+1])[[1]])
                       * 1.6022e-19 * 6.022e23/1e6)
    # update the slider range
    updateSliderInput(session, "Range", max = round(maxlimit*1.1))
    updateSliderInput(session, "Range", value = c(0, round(maxlimit*1.1)))
  })
  
  # reactively get change when input atom is changed
  spectrum = reactive({
    # added suppressWarnings to prevent "NAs introduced by coercion" warning as
    # it doesn't cause problems
    out = suppressWarnings(simulate_PES(el.input = (as.numeric(input$atom))))
  })
  
  # render the plots
  output$plotxy = renderPlot({
    # get PES data
    pes_sim_data = spectrum()
    
    # if experimental use gaus peaks is NOT selected
    if(!input$use_gaus_peaks) {
      # plot the right PES with the name and symbol as the title (histogram)
      plot_PES(pes_data=pes_sim_data,
              x_low_limit=input$Range[1], 
              x_high_limit=input$Range[2])
    } else {
      # plot the right PES with the name and symbol as the title (gaussian)
      plot_gaussian_PES(pes_data=pes_sim_data,
                    std=std_reactive(),
                    x_low_limit=input$Range[1],
                    x_high_limit=input$Range[2])
    }
  })
  
  # renders the output table of electron energies
  output$detail = renderTable({
    create_table(pes_data = spectrum(), display_console = FALSE)
  })
  
  # output function to display the avee
  output$avee = renderUI({
    HTML(paste("<b>Average valence electron energy:</b><i>",
               round((spectrum()$avee * 1.6022e-19 * 6.022e23/1e6), digits = 2),
               "(MJ/mol)</i></br></br>",
               sep = " "))
  })
  
  # store the last std value and react after 300 ms on change of std
  last_std = reactiveValues(val=0.04)
  std_reactive = debounce(reactive({input$std_val}), 300)
  
  # observe the change of std input and update the last std value
  observeEvent(input$std_val, {
    last_std$val = input$std_val
  })
  
  # render an experimental ui when the box is checked
  output$std_render = renderUI({
    if(input$use_gaus_peaks) {
      # reset the bottom margin of the checkbox (10px)
      tags$style(".checkbox { margin-bottom: 10px }")
      numericInput("std_val", "Standard Deviation:",
                   value=last_std$val, min=0.01, step=0.01)
    } else {
      # format bottom margin of checkbox to have no extra margin (-10px)
      tags$style(".checkbox { margin-bottom: -10px }")
    }
  })
})
    
