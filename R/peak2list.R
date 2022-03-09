peak2list <- function(path = getwd(), MSfile = "") {
  MSfile <- paste0(path, "/", MSfile)
  ##
  msFormat <- strsplit(MSfile, "[.]")[[1]]
  msFormat <- tolower(msFormat[length(msFormat)])
  msFormat <- gsub("/", "", msFormat)
  ##
  if ((msFormat == "mzml") | (msFormat == "mzxml")) {
    ##
    xmlData <- read_xml(MSfile)
    ##
    peakTable <- getPeakTable(xmlData, msFormat)
    ##
    spectraList <- getSpectra(xmlData, msFormat)
    ##
    p2l <- list(peakTable, spectraList)
    ##
    names(p2l) <- c("peakTable", "spectraList")
    ##
    return(p2l)
    ##
  } else {
    stop("The MSfile is not consistent with the IDSL.MXP package!")
  }
}