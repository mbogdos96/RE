# Here you can store the theme for the ggplot functions as well as the 
# functions for replacing the names of ligands from the csv files
# so that they are plotted correctly on plots and tables

# Define the replace_backbone_labels function
replace_backbone_labels <- function(x) {
  ifelse(x == "CF3-NMs-Ph", 
         "CF[3]*-NMs-H",
         ifelse(x == "CF3-NMs-OMe", 
                "CF[3]*-NMs-OMe",
                ifelse(x == "CF3-NMs-NMe2", 
                       "CF[3]*-NMs-NMe[2]",
                       ifelse(x == "F4-NMs-Ph", 
                              "F[4]*-NMs-H",
                              ifelse(x == "F4-NMs-OMe", 
                                     "F[4]*-NMs-OMe",
                                     ifelse(x == "F4-NMs-NMe2", 
                                            "F[4]*-NMs-NMe[2]",
                                            ifelse(x == "F4-NClMs-OMe", 
                                                   "F[4]*-N^Cl*Ms-OMe",
                                                   ifelse(x == "F4-NTf-OMe", 
                                                          "F[4]*-NTf-OMe", 
                                                          x)
                                            )
                                     )
                              )
                       )
                )
         )
  )
}


# Function for creating subscripts where necessary for ligand labels
# Define the replace_ligand_expression function
# This may or may not work better if you are using it in something like 
# ggrepel
replace_ligand_expression_og <- function(x) {
  ifelse(x == "PAd2Bu", 
         "PAd[2]*Bu",
         ifelse(x == "PMes3", 
                "PMes[3]",
                ifelse(x == "PtBu2Cy", 
                       "P^t*Bu[2]*Cy",
                       ifelse(x == "PtBu2Np", 
                              "P^t*Bu[2]*Np",
                              ifelse(x == "PtBu3", 
                                     "P^t*Bu[3]", x)))))
}

#Define a theme for the plots
themething <- theme(axis.line = element_line(linewidth = 1), 
                    panel.background=element_rect(fill="white"), 
                    axis.ticks = element_line(linewidth = 1), 
                    axis.text = element_text(size = rel(1.1)), 
                    axis.title = element_text(size = rel(1.1)),
                    plot.title = element_text(size = rel(1.75),
                                              face="bold"),
                    plot.title.position = "panel")
