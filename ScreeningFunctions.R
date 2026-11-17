##Function to obtain the scan numbers for peaks of interest in an XIC ##
scanMatch<-function(chromatogram, intensityCutoff, scanIndex, MSlevel){ #MSlevel = "Ms" or "Ms2"
  peakIndex<-which(chromatogram[[1]]$intensities>intensityCutoff) #Extract the vector position of the detected peaks with intensity greater than the provided value
  peakTimes<-chromatogram[[1]]$times[peakIndex] #Find the retention times (in minutes) which correspond to the selected peak intensities

  
  scanList<-scanIndex %>% #Find the scan numbers corresponding to the retention times collected above
    filter(MSOrder==MSlevel & StartTime %in% peakTimes) %>%
    pull(var=scan)
  
  if(length(scanList)<3){
    print(paste("Only", length(scanList), "MS2 scan(s) matched your criteria. Try lowering the intensity cutoff (but no lower than 1E3, or 1000).", sep=" "))
  }else if(length(scanList)>10){
    print(paste("There were", length(scanList), "MS2 scans which matched your criteria. I recommend raising the intensity cutoff and running this again.", sep=" "))
  }else{
    print(paste("There were", length(scanList), "MS2 scans which matched your criteria.", sep=" "))
  }
  
  return(scanList) #Return a vector of scan numbers, to be used with the readSpectrum() function
}

## Function to retrieve the intensities of input mz from a given spectrum ##
peakFindSpec<-function(spectrum, targetmz, ppm){ 
  mztable<-tibble(mz=spectrum$mZ, intensity=spectrum$intensity)
  returntable<-tibble(mz=NULL, intensity=NULL)
  for(i in 1:length(targetmz)) {
    mzupp<-targetmz[i]+(targetmz[i]*ppm*(10^-6))
    mzlow<-targetmz[i]-(targetmz[i]*ppm*(10^-6))
    mztablenew<-mztable %>%
      filter(mz>mzlow & mz<mzupp)
    if(nrow(mztablenew)>1){
      print(c("The mass tolerance is too wide for mz #", i))
      print("Please reduce the ppm value and try again")
      break
    }else if(nrow(mztable)<1){next}
    returntable<-bind_rows(returntable,mztablenew)
  }
  return(returntable)
}

## Function to collect potential precursor matches for a given list of MS2 scans containing the the fragment of interest
precursorFind<-function(rawdata, chromatogram, MS2_scanIDs, scantable){
  precursorTable<-tibble(frag_scanID=NULL, frag_RTmins=NULL, frag_IC=NULL, prec_scanID=NULL, prec_mz=NULL, prec_IC=NULL)
  MS1_scans<-filter(scantable, MSOrder=="Ms")
  fragChrom<-tibble(RTminutes=chromatogram$times, IC=chromatogram$intensities)
  for(j in 1:length(MS2_scanIDs)){
    print(c(j, "out of", length(MS2_scanIDs)))
    min_intensity<-fragChrom$IC[which(fragChrom$RTminutes==scantable$StartTime[MS2_scanIDs[j]])]
    
    frag_RT<-scantable$StartTime[MS2_scanIDs[j]]
    nearestMS1<-which.min(abs(MS1_scans$StartTime-frag_RT))
    nearestMS1<-MS1_scans$scan[nearestMS1]
    
    DIA_windowcenter<-scantable$precursorMass[MS2_scanIDs[j]]
    DIA_minmz<-DIA_windowcenter-5
    DIA_maxmz<-DIA_windowcenter+5
    
    MS1_spectrum<-readSpectrum(rawdata, scan=nearestMS1)
    MS1_shortlist<-tibble(frag_scanID=MS2_scanIDs[j], frag_RTmins=round(frag_RT, 4), frag_IC=min_intensity, prec_scanID=nearestMS1, prec_mz=MS1_spectrum[[1]]$mZ, prec_IC=MS1_spectrum[[1]]$intensity)
    MS1_shortlist<-MS1_shortlist %>%
      filter(between(prec_mz, DIA_minmz, DIA_maxmz)) %>%
      filter(prec_IC>=min_intensity)
    
    precursorTable<-bind_rows(precursorTable, MS1_shortlist)
  }
  
  if(nrow(precursorTable)<1){
    print("There were no potential precursors detected.")
  }else{
    print(paste("There were", nrow(precursorTable), "potential precursors detected.", sep=" "))
  }
  
  return(precursorTable)
}


##Function to do preliminary filtering of precursor results

#The idea with this function is to get rid of the most obvious mismatches, since it would be impossible to do this sorting step properly without some kind of visual confirmation
#I'm just trying to narrow down the list from hundreds to tens
#The only way I can think of is to obtain the chromatograph for each precursor mz, then check for a narrow, isolated peak at the correct retention time
  #e.g. contiguous intensity values should span no more than plus or minus 15 seconds from the fragment retention time
  #Could also check how many intensity values are greater than the peak from the desired retention time (not a dealbreaker, but if there are 10 or more higher peaks, then perhaps that could be another reason to exclude that match)
  #OR check for a local maximum within +- 1 seconds of the fragment RT, then see if that intensity is within the top 20 or so (to allow for distant unrelated peaks)

RTfilter<-function(rawdata, precList){

  precDelete<-NULL
  for(k in 1:nrow(precList)){
    if(k%%10==0){print(paste0(k, " out of ", nrow(precList)))}
    ms1Chrom<-readChromatogram(rawdata, mass=c(precList$prec_mz[k]), tol=20, type="xic", filter="ms")
    ms1Peaks<-tibble(RTminutes=round(ms1Chrom[[1]]$times, 4), IC=ms1Chrom[[1]]$intensities)
    RT<-precList$frag_RTmins[k]
    RTmin<-max(c(RT-0.0333, min(ms1Peaks$RTminutes)))
    RTmax<-min(c(RT+0.0333, max(ms1Peaks$RTminutes)))
    ICmin<-sort(ms1Peaks$IC, decreasing=TRUE)[10]
    ICmax<-max(ms1Peaks$IC)
    
    ms1Filtered<-ms1Peaks[localMaxima(ms1Peaks$IC),] %>%
      filter(between(RTminutes, RTmin, RTmax)) %>%
      filter(between(IC, ICmin, ICmax))
    
    if(nrow(ms1Filtered) < 1){
      precDelete<-c(precDelete, k)
    }
  }
  
  precList<-precList[-precDelete,]
  print("Finished!")
  return(precList)
}


##Function to obtain local maxima (source: https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima)

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}
