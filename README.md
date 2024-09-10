
V1.1 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8133971.svg)](https://doi.org/10.5281/zenodo.8133971)

v0.9 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3051109.svg)](https://doi.org/10.5281/zenodo.3051109)

# Generator of Root ANAtomy in R [GRANAR]

This package contains functions to generate and save root cross section anatomy, based on global parameters, such as the mean size of cells and number of cell layers. 

```{r}
# Load input
params <- read_param_xml(path = system.file("extdata", "root_monocot.xml", package = "granar"))
# Generate anatomy
root = create_anatomy(parameters = params)
# Visualize the simulation output
plot_anatomy(root)
# Write the simulation output
write_anatomy_xml(sim = root, path = system.file("extdata", "current_root.xml", package = "granar"))
```
