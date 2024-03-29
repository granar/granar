% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/granar_metadata.R
\name{granar_metadata}
\alias{granar_metadata}
\title{Get GRANAR metadata out of .xml}
\usage{
granar_metadata(path = NULL)
}
\arguments{
\item{path}{Path to the XML file containing granar output. Default = NULL.}
}
\description{
To get the metadata out of the xml file created with write_anatomy_xml.
}
\examples{
metadata <- granar_metadata(path = system.file("extdata", "current_root.xml", package = "granar"))

}
\keyword{anatomy}
\keyword{root}
