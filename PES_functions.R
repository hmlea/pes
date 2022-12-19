library(stringr)

# load data of binding energies and names for elements through xenon
# binding_energies = read_csv("binding_energies.csv", show_col_types=FALSE,
                            # col_types=cols(.default="d", "Element"="c"))
binding_energies = read.csv("binding_energies.csv")
el_names = read.csv("elements.csv")

# rows in binding_energies that contain binding energies or number of electrons
# binding energy is in eV (electron volts)
be_data_nrow = nrow(binding_energies) - 1
be_rows = 1:(be_data_nrow/2) + 1
ne_rows = be_rows + (be_data_nrow/2)

# function to simulate a photoelectron spectrum
simulate_PES = function(el.input){
  # filter element params (atomic no, name, symbol) to get atomic number
  if(!is.na(suppressWarnings(as.numeric(el.input)))) {
    # if atomic number is given
    at.no = el.input
  } else if(nchar(el.input) < 3) {
    # if symbol is given
    at.no = which(tolower(el_names$symbol) == tolower(el.input))
  } else {
    # if name is given
    at.no = which(tolower(el_names$name) == tolower(el.input))
  }
  
  # if at.no is null raise an error
  if(length(at.no) != 1) {
    stop(paste0("Element \"", el.input, "\" does not exist in the data"))
  }
  
  # get symbol and name from atomic number
  el.sym = as.character(el_names[at.no, 1])
  el.name = as.character(el_names[at.no, 2])
 
  # column index for element is one greater than atomic number
  element = at.no + 1
  
  # reads binding energies in eV; converts to MJ/mol; omits all NA
  # 1 eV = 1.6022e-19 J; 1 mol = 6.022e23; 1 MJ = 1e6 J
  conv = 1.6022e-19 * 6.022e23 / 1e6
  bind_energy = as.vector(na.omit(binding_energies[be_rows, element])) * conv
  
  # reads number of electrons for each orbital; omits all NA
  num_elec = as.vector(na.omit(binding_energies[ne_rows, element]))
  
  # find the outermost orbitals [refactor]
  last_orbital_cell = length(na.omit(binding_energies[be_rows, element])[[1]])+1
  last_orbital = str_sub(binding_energies[last_orbital_cell, 1][[1]], start=-2)
  last_orbital_count = binding_energies[last_orbital_cell+11, element][[1]]
  last_orbital_num = str_sub(last_orbital, start=-2, end=-2)
  last_orbital_letter = str_sub(last_orbital, start=-1)

  # use outermost orbitals to calculate average valence electron
  # energy (avee) in eV (for now) [refactor]
  if (last_orbital_letter == "s") {
    avee = binding_energies[last_orbital_cell, element][[1]]
  } else if (last_orbital_letter == "p") {
    if (last_orbital_num == 2 || last_orbital_num == 3) relative_s_loc = 1
    else relative_s_loc = 2
    
    num_s_elec = binding_energies[last_orbital_cell+(11-relative_s_loc), element][[1]]
    total_s_elec_energy = binding_energies[last_orbital_cell-relative_s_loc, element][[1]] * num_s_elec
    
    num_p_elec = binding_energies[last_orbital_cell+11, element][[1]]
    total_p_elec_energy = binding_energies[last_orbital_cell, element][[1]] * num_p_elec
    
    avee = (total_s_elec_energy + total_p_elec_energy) / (num_s_elec + num_p_elec)
  } else if (last_orbital_letter == "d") {
    num_s_elec = binding_energies[last_orbital_cell+10, element][[1]]
    total_s_elec_energy = binding_energies[last_orbital_cell-1, element][[1]] * num_s_elec
    
    num_d_elec = binding_energies[last_orbital_cell+11, element][[1]]
    total_d_elec_energy = binding_energies[last_orbital_cell, element][[1]] * num_d_elec
    
    avee = (total_s_elec_energy + total_d_elec_energy) / (num_s_elec + num_d_elec)
  }
  
  # output a list of items
  output = list("atomic_number" = at.no,
                "element_symbol" = el.sym,
                "element_name" = str_to_title(el.name),
                "binding_energy" = na.omit(bind_energy), 
                "number_electrons" = na.omit(num_elec),
                "avee" = avee)
}

# plot the PES data simulated using a histogram-type plot
plot_PES = function(pes_data, x_low_limit=NULL, x_high_limit=NULL,
                   y_limit=NULL, plot_title=NULL) {
  # if x_low_limit is null, use 0
  if(is.null(x_low_limit)) {
    x_low_limit = 0
  }
  # if x_high_limit is null, use maximum binding energy
  if(is.null(x_high_limit)){
    x_high_limit = max(pes_data$binding_energy)
  }
  # if y_limit is null, use the maximum number of electrons
  if(is.null(y_limit)) {
    max_num_electrons = max(pes_data$number_electrons)
    y_limit = ifelse(max_num_electrons %% 2 == 0,
                   max_num_electrons + 2,
                   max_num_electrons + 1)
  }
  # if plot_title is null, create a title; format: Name (Symbol)
  if(is.null(plot_title)) {
    plot_title = paste0(pes_data$element_name, " (",
                        pes_data$element_symbol, ")")
  }
  
  # plot the spectrum
  plot(x=pes_data$binding_energy, y=pes_data$number_electrons, 
       type="h", xlim=c(x_high_limit, x_low_limit), ylim=c(0, y_limit), 
       xlab="Binding Energy (MJ/mol)", ylab="Electron Count", 
       lwd=2, main=plot_title)
}

# function to create gaussian peaks from the simulated PES data
# resembles a PES a little more than the histogram (std = peak width)
# WORK IN PROGRESS
create_gaussian_peaks = function(pes_data, std, x_interval,
                        lower_bound=NULL, upper_bound=NULL) {
  # display nothing while std is null
  if(is.null(std)) {
    plot(0, 0, type="n")
    return()
  }
  
  # get data for gaus peak creation
  positions = pes_data$binding_energy
  counts = pes_data$number_electrons
  
  # if bounds are null
  if(is.null(lower_bound)) lower_bound = 0
  if(is.null(upper_bound)) {
    max_pos = max(positions)
    upper_bound = ceiling(max_pos * 1.1)
  }
  
  # create the initial x and y sequences
  x_seq = seq(lower_bound, upper_bound, 0.01)
  y_val = rep(0, length(x_seq))
  
  # for each electron
  for(i in seq_along(positions)) {
    # get its position (relative binding energy) and count (num of electrons)
    cur_pos = positions[[i]]
    cur_count = counts[[i]]
    
    # create a peak by evaluating the gaussian expression
    str_expr = paste0(cur_count, "*exp(-(((x-", cur_pos, ")^2)/(2*(", std, "^2))))")
    peak = eval(parse(text=str_expr), env=list(x=x_seq))
    # add the peak to the baseline y_val of 0
    y_val = y_val + peak
  }
  
  # output a dataframe of x and y values (ready to plot pretty much)
  output = data.frame("x"=x_seq, "y"=y_val)
}

# plots the simulated PES data using the gausian function above
# WORK IN PROGRESS
plot_gaussian_PES = function(pes_data, std=0.04, x_interval=0.01,
                         x_low_limit=NULL, x_high_limit=NULL,
                         y_limit=NULL, plot_title=NULL) {
  # if x_low_limit is null, use 0
  if(is.null(x_low_limit)) {
    x_low_limit = 0
  }
  # if x_high_limit is null, use maximum binding energy
  if(is.null(x_high_limit)){
    x_high_limit = max(pes_data$binding_energy)
  }
  # if y_limit is null, use the maximum number of electrons
  if(is.null(y_limit)) {
    max_num_electrons = max(pes_data$number_electrons)
    y_limit = ifelse(max_num_electrons %% 2 == 0,
                     max_num_electrons + 2,
                     max_num_electrons + 1)
  }
  # if plot_title is null, create a title; format: Name (Symbol)
  if(is.null(plot_title)) {
    plot_title = paste0(pes_data$element_name, " (",
                        pes_data$element_symbol, ")")
  }  

  # get the gausian peak data
  gaus_peaks = create_gaussian_peaks(pes_data, std, x_interval,
                                 lower_bound=x_low_limit,
                                 upper_bound=x_high_limit)
  
  # plot the spectrum
  plot(gaus_peaks, type = "l", xlim = c(x_high_limit, x_low_limit),
       ylim=c(0, y_limit), xlab="Binding Energy (MJ/mol)",
       ylab="Electron Count", main=plot_title)
}

# create a table of electron binding energies
create_table = function(pes_data, display_console=TRUE){
  # get a list of all orbitals (currently in binding_energy data)
  atomic_orbital = unlist(lapply(as.character(binding_energies$Element[be_rows]), 
                                 function(x) substr(x, nchar(x)-1, nchar(x)+1)))
  
  # get orbitals, their binding energies, and electron counts
  number_orbital = length(pes_data$binding_energy)
  peaks = round(pes_data$binding_energy, digits = 2)
  counts = as.integer(pes_data$number_electrons)
  
  # create the table data frame
  peak_table = data.frame(atomic_orbital[1:number_orbital], peaks, counts)
  colnames(peak_table) = c("Atomic Orbital", "Binding Energy (MJ/mol)",
                           "Number of Electrons")
  
  # display (if need be) and output the table
  if (display_console == TRUE) show(peak_table)
  output = peak_table
}

