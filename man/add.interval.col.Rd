% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/001-add.interval.col.R
\name{add.interval.col}
\alias{add.interval.col}
\title{Add columns for calculations within PKNCA intervals}
\usage{
add.interval.col(name, FUN, values, depends = c(), desc = "",
  datatype = c("interval", "individual", "population"))
}
\arguments{
\item{name}{The column name}

\item{FUN}{The function to run (as a character string)}

\item{values}{Valid values for the column}

\item{depends}{Character vector of columns that must be run before
this column.}

\item{desc}{A human-readable description of the parameter (<=40
characters to comply with SDTM)}

\item{datatype}{The type of data}
}
\description{
Add columns for calculations within PKNCA intervals
}
