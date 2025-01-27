# Generator of Root ANAtomy in R (GRANAR)

v1.1 --> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8133971.svg)](https://doi.org/10.5281/zenodo.8133971)

GRANAR is an open-source R package designed to generate root cross-section anatomical networks with ease. Using global parameters, such as mean cell size and number of cell layers, GRANAR simplifies the creation and visualization of virtual monocotyledon root anatomies. This tool is particularly useful for researchers studying root hydraulic conductivity and the impact of anatomical features on water transport.

---

## Features

- **Generate Root Anatomies**: Create root cross-sections using customizable parameters.
- **Visualize Results**: Easily plot and explore generated anatomies.
- **Save & Share Outputs**: Export anatomical data for further analysis.
- **Integration with Hydraulic Models**: Seamlessly use GRANAR outputs with tools like [MECHA](https://github.com/MECHARoot/MECHA) for hydraulic property estimations.

## Installation

Install the latest version of GRANAR directly from GitHub using `devtools`:

```R
# Install devtools if not already installed
install.packages("devtools")

# Install GRANAR
devtools::install_github("granar/granar")
```


## Getting Started

Below is a quick example to demonstrate how to use GRANAR:

```R
# Load GRANAR
library(granar)

# Load input parameters
params <- read_param_xml(path = system.file("extdata", "root_monocot.xml", package = "granar"))

# Generate root anatomy
root <- create_anatomy(parameters = params)

# Visualize the generated anatomy
plot_anatomy(root)

# Save the generated anatomy to an XML file
write_anatomy_xml(sim = root, path = "current_root.xml")
```


## Documentation

For detailed documentation and usage examples, please visit the [GRANAR website](http://granar.github.io).

## Support

If you encounter any issues, have questions, or wish to suggest enhancements, feel free to [raise an issue](https://github.com/granar/granar/issues/new?template=Blank+issue).

## Contributing

Contributions are welcome! If you'd like to contribute to GRANAR, please:

1. Fork the repository.
2. Make your changes in a new branch.
3. Submit a pull request with a clear description of your changes.

## Citation

If you use GRANAR in your research, please cite the following paper:

> Heymans A., Couvreur V., LaRue T., Paez-Garcia A., Lobet G. (2020). GRANAR, a computational tool to better understand the functional importance of monocotyledon root anatomy. *Plant Physiology*, 182(2), 707-720. doi: [10.1104/pp.19.00617](https://doi.org/10.1104/pp.19.00617)



## License

GRANAR is released under a GNU General Public Licence v3.0, which means that distribution and use and binary forms, with or without modification, are permitted under the GNU General Public License v3 and provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

See the [LICENSE](https://github.com/granar/granar/blob/main/LICENSE) file for details.

## Versioning

>  v0.9 --> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3051109.svg)](https://doi.org/10.5281/zenodo.3051109)
