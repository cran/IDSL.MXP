getPeakTable <- function(xmlData, msFormat) {
  peakTable <- NA
  if (tolower(msFormat) == "mzml") {
    spectrumNodes <- xml_find_all(xmlData, '//d1:spectrum')
    xmlNSspectrumNodes <- xml_ns(spectrumNodes)
    ############################################################################
    positiveScanNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="positive scan"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    negativeScanNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="negative scan"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    centroidSpectrumNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="centroid spectrum"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    profileSpectrumNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="profile spectrum"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    basePeakMZNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="base peak m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    basePeakIntensityNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="base peak intensity"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    totIonCurrentNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="total ion current"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    msLevelNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="ms level"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    lowMZNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="lowest observed m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    highMZNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="highest observed m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    retentionTimeNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    filterStringNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:cvParam[@name="filter string"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    injectionTimeNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:cvParam[@name="ion injection time"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    scanWindowLowerLimitNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:scanWindowList/d1:scanWindow/d1:cvParam[@name="scan window lower limit"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    scanWindowUpperLimitNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:scanWindowList/d1:scanWindow/d1:cvParam[@name="scan window upper limit"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    isolationWindowTargetMZNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:isolationWindow/d1:cvParam[@name="isolation window target m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    isolationWindowLowerOffsetNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:isolationWindow/d1:cvParam[@name="isolation window lower offset"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    isolationWindowUpperOffsetNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:isolationWindow/d1:cvParam[@name="isolation window upper offset"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    precursorScanNumNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor', ns = xmlNSspectrumNodes, flatten = FALSE)
    precursorMZNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:selectedIonList/d1:selectedIon/d1:cvParam[@name="selected ion m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    precursorChargeNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:selectedIonList/d1:selectedIon/d1:cvParam[@name="charge state"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    precursorIntensityNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:selectedIonList/d1:selectedIon/d1:cvParam[@name="peak intensity"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    collisionEnergyNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:activation/d1:cvParam[@name="collision energy"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    peakTable <- do.call(rbind, lapply(1:length(spectrumNodes), function(i) {
      S <- xml_attrs(spectrumNodes[[i]])
      ##
      seqNum <- as.numeric(S[1]) + 1
      ##
      spectrumId <- S[2]
      if (length(spectrumId) == 0) {
        spectrumId <- NA
      }
      ##
      x_acquisitionNum <- MXP_locate_regex(spectrumId, "=")
      x_acquisitionNum <- x_acquisitionNum[nrow(x_acquisitionNum), 2]
      acquisitionNum <- substr(spectrumId, (x_acquisitionNum + 1), nchar(spectrumId))
      if (length(acquisitionNum) == 0) {
        acquisitionNum <- NA
      }
      ##
      peaksCount <- S[3]
      if (length(peaksCount) == 0) {
        peaksCount <- 0
      }
      ##########################################################################
      positiveScan <- xml_has_attr(positiveScanNodes[[i]], "value")
      ##
      negativeScan <- xml_has_attr(negativeScanNodes[[i]], "value")
      ##
      if (length(negativeScan) > 0) {
        if (negativeScan) {
          polarity <- 0
        }
      } else if (length(positiveScan) > 0) {
        if (positiveScan) {
          polarity <- 1
        }
      } else {
        polarity <- -1
      }
      ##########################################################################
      centroidSpectrum <- xml_has_attr(centroidSpectrumNodes[[i]], "value")
      ##
      profileSpectrum <- xml_has_attr(profileSpectrumNodes[[i]], "value")
      ##
      if (length(centroidSpectrum) > 0) {
        if (centroidSpectrum) {
          centroided <- TRUE
        } else {
          centroided <- FALSE
        }
      } else if (length(profileSpectrum) > 0) {
        centroided <- FALSE
      } else {
        centroided <- FALSE
      }
      ##########################################################################
      basePeakMZ <- xml_attr(basePeakMZNodes[[i]], "value")
      if (length(basePeakMZ) == 0) {
        basePeakMZ <- NA
      }
      ##
      basePeakIntensity <- xml_attr(basePeakIntensityNodes[[i]], "value")
      if (length(basePeakIntensity) == 0) {
        basePeakIntensity <- NA
      }
      ##
      totIonCurrent <- xml_attr(totIonCurrentNodes[[i]], "value")
      if (length(totIonCurrent) == 0) {
        totIonCurrent <- NA
      }
      ##
      msLevel <- xml_attr(msLevelNodes[[i]], "value")
      if (length(msLevel) == 0) {
        msLevel <- 0
      }
      ##
      lowMZ <- xml_attr(lowMZNodes[[i]], "value")
      if (length(lowMZ) == 0) {
        lowMZ <- NA
      }
      ##
      highMZ <- xml_attr(highMZNodes[[i]], "value")
      if (length(highMZ) == 0) {
        highMZ <- NA
      }
      ##########################################################################
      retentionTime <- xml_attr(retentionTimeNodes[[i]], "value")
      if (length(retentionTime) == 0) {
        retentionTime <- NA
      }
      ##
      filterString <- xml_attr(filterStringNodes[[i]], "value")
      if (length(filterString) == 0) {
        filterString <- NA
      }
      ##
      injectionTime <- xml_attr(injectionTimeNodes[[i]], "value")
      if (length(injectionTime) == 0) {
        injectionTime <- NA
      }
      ##########################################################################
      scanWindowLowerLimit <- xml_attr(scanWindowLowerLimitNodes[[i]], "value")
      if (length(scanWindowLowerLimit) == 0) {
        scanWindowLowerLimit <- NA
      }
      ##
      scanWindowUpperLimit <- xml_attr(scanWindowUpperLimitNodes[[i]], "value")
      if (length(scanWindowUpperLimit) == 0) {
        scanWindowUpperLimit <- NA
      }
      ##########################################################################
      isolationWindowTargetMZ <- xml_attr(isolationWindowTargetMZNodes[[i]], "value")
      if (length(isolationWindowTargetMZ) == 0) {
        isolationWindowTargetMZ <- NA
      }
      ##
      isolationWindowLowerOffset <- xml_attr(isolationWindowLowerOffsetNodes[[i]], "value")
      if (length(isolationWindowLowerOffset) == 0) {
        isolationWindowLowerOffset <- NA
      }
      ##
      isolationWindowUpperOffset <- xml_attr(isolationWindowUpperOffsetNodes[[i]], "value")
      if (length(isolationWindowUpperOffset) == 0) {
        isolationWindowUpperOffset <- NA
      }
      ##########################################################################
      precursorScanNum <- xml_attr(precursorScanNumNodes[[i]], "spectrumRef")
      if (length(precursorScanNum) > 0) {
        if (!is.na(precursorScanNum)) {
          x_scan <- MXP_locate_regex(precursorScanNum, "scan=")
          precursorScanNum <- substr(precursorScanNum, (x_scan[2] + 1), nchar(precursorScanNum))
        } else {
          precursorScanNum <- NA
        }
      } else {
        precursorScanNum <- NA
      }
      ##
      precursorMZ <- xml_attr(precursorMZNodes[[i]], "value")
      if (length(precursorMZ) == 0) {
        precursorMZ <- NA
      }
      ##
      precursorCharge <- xml_attr(precursorChargeNodes[[i]], "value")
      if (length(precursorCharge) == 0) {
        precursorCharge <- NA
      }
      ##
      precursorIntensity <- xml_attr(precursorIntensityNodes[[i]], "value")
      if (length(precursorIntensity) == 0) {
        precursorIntensity <- NA
      }
      ##########################################################################
      collisionEnergy <- xml_attr(collisionEnergyNodes[[i]], "value")
      if (length(collisionEnergy) == 0) {
        collisionEnergy <- NA
      }
      ##########################################################################
      c(seqNum, acquisitionNum, msLevel, polarity, peaksCount, totIonCurrent, retentionTime,
        basePeakMZ, basePeakIntensity, collisionEnergy, lowMZ, highMZ, precursorScanNum, precursorMZ,
        precursorCharge, precursorIntensity, injectionTime, filterString, spectrumId, centroided,
        isolationWindowTargetMZ, isolationWindowLowerOffset, isolationWindowUpperOffset,
        scanWindowLowerLimit, scanWindowUpperLimit)
    }))
    ############################################################################
    peakTable <- data.frame(peakTable)
    colnames(peakTable) <- c("seqNum", "acquisitionNum", "msLevel", "polarity", "peaksCount", "totIonCurrent", "retentionTime",
                             "basePeakMZ", "basePeakIntensity", "collisionEnergy", "lowMZ", "highMZ", "precursorScanNum", "precursorMZ",
                             "precursorCharge", "precursorIntensity", "injectionTime", "filterString", "spectrumId", "centroided",
                             "isolationWindowTargetMZ", "isolationWindowLowerOffset", "isolationWindowUpperOffset",
                             "scanWindowLowerLimit", "scanWindowUpperLimit")
    ##
    peakTable$acquisitionNum <- as.numeric(peakTable$acquisitionNum)
    ############################################################################
  } else if (tolower(msFormat) == "mzxml") {
    ##
    scanNodes <- xml_find_all(xmlData, '//d1:scan')
    ##
    peakTable <- do.call(rbind, lapply(1:length(scanNodes), function(i) {
      S <- xml_attrs(scanNodes[[i]])
      ##########################################################################
      nameS <- names(S)
      ##########################################################################
      x_num <- which(nameS == "num")
      if (length(x_num) > 0) {
        seqNum <- S[x_num]
      } else {
        seqNum <- 0
      }
      ##########################################################################
      x_scanType <- which(nameS == "scanType")
      if (length(x_scanType) > 0) {
        scanType <- S[x_scanType]
      } else {
        scanType <- NA
      }
      ##########################################################################
      x_peaksCount <- which(nameS == "peaksCount")
      if (length(x_peaksCount) > 0) {
        peaksCount <- S[x_peaksCount]
      } else {
        peaksCount <- 0
      }
      ##########################################################################
      x_polarity <- which(nameS == "polarity")
      if (length(x_polarity) > 0) {
        pol <- S[x_polarity]
        if (pol == "-" | grepl("N", pol, ignore.case = TRUE)) {
          polarity <- 0
        } else if (pol == "+" | grepl("P", pol, ignore.case = TRUE)) {
          polarity <- 1
        }
      } else {
        polarity <- -1
      }
      ##########################################################################
      x_centroided <- which(nameS == "centroided")
      if (length(x_centroided) > 0) {
        cent <- as.numeric(S[x_centroided])
        if (cent == 1) {
          centroided <- TRUE
        } else if (cent == 0) {
          centroided <- FALSE
        }
      } else {
        centroided <- FALSE
      }
      ##########################################################################      
      x_basePeakMZ <- which(nameS == "basePeakMZ")
      if (length(x_basePeakMZ) > 0) {
        basePeakMZ <- S[x_basePeakMZ]
      } else {
        basePeakMZ <- NA
      }
      ##########################################################################
      x_basePeakIntensity <- which(nameS == "basePeakIntensity")
      if (length(x_basePeakIntensity) > 0) {
        basePeakIntensity <- S[x_basePeakIntensity]
      } else {
        basePeakIntensity <- NA
      }
      ##########################################################################
      x_totIonCurrent <- which(nameS == "totIonCurrent")
      if (length(x_totIonCurrent) > 0) {
        totIonCurrent <- S[x_totIonCurrent]
      } else {
        totIonCurrent <- 0
      }
      ##########################################################################
      x_msLevel <- which(nameS == "msLevel")
      if (length(x_msLevel) > 0) {
        msLevel <- S[x_msLevel]
      } else {
        msLevel <- 0
      }
      ##########################################################################
      x_lowMZ <- which(nameS == "lowMZ")
      if (length(x_lowMZ) > 0) {
        lowMZ <- S[x_lowMZ]
      } else {
        lowMZ <- NA
      }
      ##########################################################################
      x_highMZ <- which(nameS == "highMZ")
      if (length(x_highMZ) > 0) {
        highMZ <- S[x_highMZ]
      } else {
        highMZ <- NA
      }
      ##########################################################################
      x_retentionTime <- which(nameS == "retentionTime")
      if (length(x_retentionTime) > 0) {
        retentionTime <- S[x_retentionTime]
      } else {
        retentionTime <- NA
      }
      ##########################################################################
      x_filterString <- which(nameS == "filterString")
      if (length(x_filterString) > 0) {
        filterString <- S[x_filterString]
      } else {
        filterString <- NA
      }
      ##########################################################################
      x_injectionTime <- which(nameS == "injectionTime")
      if (length(x_injectionTime) > 0) {
        injectionTime <- S[x_injectionTime]
      } else {
        injectionTime <- NA
      }
      ##########################################################################
      x_scanWindowLowerLimit <- which(nameS == "scanWindowLowerLimit")
      if (length(x_scanWindowLowerLimit) > 0) {
        scanWindowLowerLimit <- S[x_scanWindowLowerLimit]
      } else {
        scanWindowLowerLimit <- NA
      }
      ##########################################################################
      x_scanWindowUpperLimit <- which(nameS == "scanWindowUpperLimit")
      if (length(x_scanWindowUpperLimit) > 0) {
        scanWindowUpperLimit <- S[x_scanWindowUpperLimit]
      } else {
        scanWindowUpperLimit <- NA
      }
      ##########################################################################
      x_isolationWindowTargetMZ <- which(nameS == "isolationWindowTargetMZ")
      if (length(x_isolationWindowTargetMZ) > 0) {
        isolationWindowTargetMZ <- S[x_isolationWindowTargetMZ]
      } else {
        isolationWindowTargetMZ <- NA
      }
      ##########################################################################
      x_isolationWindowLowerOffset <- which(nameS == "isolationWindowLowerOffset")
      if (length(x_isolationWindowLowerOffset) > 0) {
        isolationWindowLowerOffset <- S[x_isolationWindowLowerOffset]
      } else {
        isolationWindowLowerOffset <- NA
      }
      ##########################################################################
      x_isolationWindowUpperOffset <- which(nameS == "isolationWindowUpperOffset")
      if (length(x_isolationWindowUpperOffset) > 0) {
        isolationWindowUpperOffset <- S[x_isolationWindowUpperOffset]
      } else {
        isolationWindowUpperOffset <- NA
      }
      ##########################################################################
      x_precursorScanNum <- which(nameS == "precursorScanNum")
      if (length(x_precursorScanNum) > 0) {
        precursorScanNum <- S[x_precursorScanNum]
      } else {
        precursorScanNum <- NA
      }
      ##########################################################################
      x_precursorMZ <- which(nameS == "precursorMZ")
      if (length(x_precursorMZ) > 0) {
        precursorMZ <- S[x_precursorMZ]
      } else {
        precursorMZ <- NA
      }
      ##########################################################################
      x_precursorCharge <- which(nameS == "precursorCharge")
      if (length(x_precursorCharge) > 0) {
        precursorCharge <- S[x_precursorCharge]
      } else {
        precursorCharge <- NA
      }
      ##########################################################################
      x_precursorIntensity <- which(nameS == "precursorIntensity")
      if (length(x_precursorIntensity) > 0) {
        precursorIntensity <- S[x_precursorIntensity]
      } else {
        precursorIntensity <- NA
      }
      ##########################################################################
      x_collisionEnergy <- which(nameS == "collisionEnergy")
      if (length(x_collisionEnergy) > 0) {
        collisionEnergy <- S[x_collisionEnergy]
      } else {
        collisionEnergy <- NA
      }
      ##########################################################################
      c(seqNum, msLevel, polarity, peaksCount, totIonCurrent, retentionTime, basePeakMZ,
        basePeakIntensity, collisionEnergy, lowMZ, highMZ, precursorScanNum, precursorMZ,
        precursorCharge, precursorIntensity, injectionTime, filterString, scanType, centroided,
        isolationWindowTargetMZ, isolationWindowLowerOffset, isolationWindowUpperOffset,
        scanWindowLowerLimit, scanWindowUpperLimit)
    }))
    ##
    peakTable <- data.frame(peakTable)
    colnames(peakTable) <- c("seqNum", "msLevel", "polarity", "peaksCount", "totIonCurrent", "retentionTime", "basePeakMZ",
                             "basePeakIntensity", "collisionEnergy", "lowMZ", "highMZ", "precursorScanNum", "precursorMZ",
                             "precursorCharge", "precursorIntensity", "injectionTime", "filterString", "scanType", "centroided",
                             "isolationWindowTargetMZ", "isolationWindowLowerOffset", "isolationWindowUpperOffset",
                             "scanWindowLowerLimit", "scanWindowUpperLimit")
    ##
    peakTable$retentionTime <- as.numeric(gsub("[a-zA-Z]", "", peakTable$retentionTime, ignore.case = TRUE))/60
    ##
  } else {
    stop("The MSfile is not consistent with the IDSL.MXP package!")
  }
  ##############################################################################
  peakTable$seqNum <- as.numeric(peakTable$seqNum)
  peakTable$msLevel <- as.numeric(peakTable$msLevel)
  peakTable$polarity <- as.numeric(peakTable$polarity)
  peakTable$peaksCount <- as.numeric(peakTable$peaksCount)
  peakTable$totIonCurrent <- as.numeric(peakTable$totIonCurrent)
  peakTable$retentionTime <- as.numeric(peakTable$retentionTime)
  peakTable$basePeakMZ <- as.numeric(peakTable$basePeakMZ)
  peakTable$basePeakIntensity <- as.numeric(peakTable$basePeakIntensity)
  peakTable$lowMZ <- as.numeric(peakTable$lowMZ)
  peakTable$highMZ <- as.numeric(peakTable$highMZ)
  peakTable$precursorScanNum <- as.numeric(peakTable$precursorScanNum)
  peakTable$precursorMZ <- as.numeric(peakTable$precursorMZ)
  peakTable$precursorIntensity <- as.numeric(peakTable$precursorIntensity)
  peakTable$injectionTime <- as.numeric(peakTable$injectionTime)
  peakTable$isolationWindowTargetMZ <- as.numeric(peakTable$isolationWindowTargetMZ)
  peakTable$isolationWindowLowerOffset <- as.numeric(peakTable$isolationWindowLowerOffset)
  peakTable$isolationWindowUpperOffset <- as.numeric(peakTable$isolationWindowUpperOffset)
  peakTable$scanWindowLowerLimit <- as.numeric(peakTable$scanWindowLowerLimit)
  peakTable$scanWindowUpperLimit <- as.numeric(peakTable$scanWindowUpperLimit)
  ##############################################################################
  return(peakTable)
}