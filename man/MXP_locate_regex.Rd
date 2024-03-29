\name{MXP_locate_regex}
\alias{MXP_locate_regex}
\title{
MXP Locate regex
}
\description{
Locate indices of the pattern in the string
}
\usage{
MXP_locate_regex(string, pattern, ignore.case = FALSE, perl = FALSE, fixed = FALSE,
useBytes = FALSE)
}
\arguments{
  \item{string}{
a string as character
}
  \item{pattern}{
a pattern to screen
}
  \item{ignore.case}{
ignore.case
}
  \item{perl}{
perl
}
  \item{fixed}{
fixed
}
  \item{useBytes}{
useBytes
}
}
\details{
This function returns `NULL` when no matches are detected for the pattern.
}
\value{
A 2-column matrix of location indices. The first and second columns represent start and end positions, respectively.
}
\examples{
pattern <- "Cl"
string <- "NaCl.5HCl"
Location_Cl <- MXP_locate_regex(string, pattern)
}