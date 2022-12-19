# create a list for use in the ui atom selection input
choice_input_list = function() {
  element_names = read.csv("elements.csv")
  choice_list = as.list(as.numeric(row.names(element_names)))
  names(choice_list) = element_names$symbol
  choice_list
}

