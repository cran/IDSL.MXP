peak2list <- function(path = getwd(), MSfile = "") {
  ##
  MSfile <- paste0(path, "/", MSfile)
  if(!file.exists(MSfile)) {
    MSfile <- substr(MSfile, 1, (nchar(MSfile) - 1))
  }
  ##
  msFormat <- strsplit(MSfile, "[.]")[[1]]
  msFormat <- tolower(msFormat[length(msFormat)])
  msFormat <- gsub("/", "", msFormat)
  ##
  if ((msFormat == "mzml") | (msFormat == "mzxml")) {
    ##
    xmlData <- read_xml(MSfile)
    ##
    scanTable <- getScanTable(xmlData, msFormat)
    ##
    spectraList <- getSpectra(xmlData, msFormat)
    ##
    p2l <- list(scanTable, spectraList)
    ##
    names(p2l) <- c("scanTable", "spectraList")
    ##
  } else if ((msFormat == "cdf")) {
    ##
    p2l <- getNetCDF(MSfile)
    ##
  } else {
    stop("The MSfile is not consistent with the IDSL.MXP package!")
  }
  return(p2l)
}