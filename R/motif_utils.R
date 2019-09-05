motifDb_to_tfbs_pwmMatrixList <- function(motifDb, pseudocounts = 0.8){
  # returns TFBSTools-PFMatrixList from motifDb entries
  motifs <- universalmotif::convert_motifs(motifDb, out_class = "universalmotif")
  
  tfbs_list <- lapply(motifs, function(motif){
    universalmotif::convert_motifs(motif, out_class = "TFBSTools-PFMatrix")
  })
  
  tfbs_PWMatrixList <- do.call(TFBSTools::PWMatrixList, lapply(tfbs_list, TFBSTools::toPWM, pseudocounts = pseudocounts))
}


compute_background_nucleotide_frequency <- function(background_regions, 
                                                    BSgenome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3){
  # Take background GRanges regions, and calculate background nucleotide frequencies
  # for use with motifMatchr(bg = bg_freq)
  background <- Biostrings::getSeq(BSgenome,
                     names = background_regions)
  
  bg_freq <- Biostrings::oligonucleotideFrequency(background, width = 1, as.prob = T, simplify.as = "collapsed")
  
  return(bg_freq) 
}

get_sequence_from_regions <- function(regions, genome){
   chrNames <- seqnames(regions) 
   startPos <- start(regions)
   endPos <- end(regions)
   
   feature_names <- paste0(chrNames,":",startPos,"-",endPos)
   
   sequences <- Biostrings::getSeq(genome, regions)
   names(sequences) <- feature_names
   return(sequences)
}
#library(MotifDb)
#library(motifmatchr)
#library(TFBSTools)
#library(BSgenome.Dmelanogaster.UCSC.dm3)
#
## Need to call "subset"
#flyFactor.motifDb <- subset(MotifDb, dataSource == "FlyFactorSurvey")

write_fasta_from_GRanges <- function(regions, resize_bp, 
                                     name = deparse(substitute(regions)), 
                                     genome =  BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3,
                                     outdir){
  # deparse(substitute(regions)) is a hacky way of using the variable name in the filename
  
   nseq <- length(regions)
   fasta_name <- paste0(outdir, name, "_", nseq, "_features", ".fasta") 
  
   center_flank <- resize(regions, 1, "center") %>%
              promoters(., resize_bp, resize_bp)
   
   seq <- get_sequence_from_regions(center_flank, genome)
   
   R.utils::mkdirs(outdir)
   Biostrings::writeXStringSet(seq, fasta_name, format = "fasta")
   return(NULL)
}

write_fasta_from_list <- function(regionsList, resize_bp, genome, outdir) {
  lapply(names(regionsList), function(name){
    # write fasta files of all regions named after the type of feature
   regions <- regionsList[[name]] 
   write_fasta_from_GRanges(regions, resize_bp, name, outdir, genome) 
  })
  return(NULL)
}

write_fasta_from_list_dm3 <- function(regionsList, resize_bp, genome, outdir) {
  lapply(names(regionsList), function(name){
    # write fasta files of all regions named after the type of feature
   regions <- regionsList[[name]] 
   write_fasta_from_GRanges(regions, resize_bp, name, outdir, BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3) 
  })
  return(NULL)
}

write_fasta_of_region <- function(regionsList, 
                      outdir,
                      resize_bp = 200, 
                      genome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3){
    # takes regions, resizes to symmetric about resize_bp distance of center
    # regionsList is named list of types of features
    
    # TODO: check that resize_bp is not negative 
  
    # check outdir ends with /, else append /:
    if (stringr::str_sub(outdir, -1) != "/") {
      outdir <- paste0(outdir, "/")  
    }
     
    R.utils::mkdirs(outdir)
    fasta_names <- lapply(names(regionsList), function(name){
      # write fasta files of all regions named after the type of feature
      regions <- regionsList[[name]] 
      nseq <- length(regions)
      fasta_name <- paste0(outdir, name, "_", nseq, "_features", ".fasta") 
     
      if (resize_bp > 0) {
        center_flank <- resize(regions, 1, "center") %>%
          promoters(., resize_bp, resize_bp)
      } else {
        center_flank <- regions
      }
      
      seq <- get_sequence_from_regions(center_flank, genome)
      
      Biostrings::writeXStringSet(seq, fasta_name, format = "fasta")
      return(fasta_name) 
    })
    
    names(fasta_names) <- names(regionsList)
    return(fasta_names) 
}


write_fasta_from_overlaps <- function(overlappingPeaks, by, outdir, resize_bp, 
         genome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3,
         summit_pos_colname = "summitPos"){
  # Write fasta sequence of each set of unique&merged peaks from chip-peak anno overlap
  # where `by` = colname to compare merged peaks by
  #
  # for proper motif analysis of ov need to do the following:
  # make overlap object, dropping all metadata except grp/peakDatasource & summitPos
  # export uniquePeaks (split on grp/peak_type)
  # export merged peaks that overlap with original peaks
  # merged peak exports will use summits from both datasets, so there will be 2 files one using set1 other with set2 summits
  # 
  
  if (class(overlappingPeaks) != "overlappingPeaks"){ stop("Error: overlappingPeaks must be of class 'overlappingPeaks'")} 
  
  uniquePeaks <- overlappingPeaks$uniquePeaks %>% 
    split(., mcols(.)[[by]])
  
  
  mergedPeaks <- overlappingPeaks$mergedPeaks
  allPeaks <- overlappingPeaks$all.peaks
  
  mergedPeaksFromSource <- lapply(names(allPeaks), function(name){
    sourcePeaks <- allPeaks[[name]]
    
    peaks_ov_merged <- subsetByOverlaps(sourcePeaks, mergedPeaks)
    return(peaks_ov_merged)
  }) 
  
  names(mergedPeaksFromSource) <- paste0("mergedPeaks_", names(allPeaks), ".summits")
  
  # compute summits 
  uniquePeaks_summits <- lapply(uniquePeaks, computeSummits, summitColumn = summit_pos_colname)
  names(uniquePeaks_summits) <- paste0("uniquePeaks_", names(uniquePeaks_summits), ".summits") 
  
  mergedPeaksFromSource_summits <- lapply(mergedPeaksFromSource, computeSummits, summitColumn = summit_pos_colname) 
  
  # extend about summits and write regions
  write_fasta_of_region(uniquePeaks_summits, outdir, resize_bp, genome)
  write_fasta_of_region(mergedPeaksFromSource_summits, outdir, resize_bp, genome)
  return(NULL) # try returning file paths or something
}

compute_letter_frequency <- function(peaks, letters, as.prob = F,
                               genome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3){
  
  # takes GRanges as input, returns vector of letter frequencies within each region of the GR object
  # ex: computing GC content letters  = "GC". 
  # output will be reported as fraction of the sequence with G or C
  freq <- Biostrings::getSeq(genome, peaks) %>% 
    Biostrings::letterFrequency(., letters, as.prob = as.prob) %>% 
    as.vector
  
  if (length(freq) != length(peaks)) {stop("ERROR: output vector is not same length as input peak list")}
  else {return(freq)}
  
}


compute_gc_content <- function(peaks, 
                               genome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3){
  # special case of compute_letter_frequency
  compute_letter_frequency(peaks, "GC", as.prob = T, genome)
}


countMotifMatches <- function(peaks, pwm, colname, minScore = "80%", resize_bp = 200,
                              genome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3){
  # counts motif matches inside peaks, returns peaklist w/ column of counts appended to it (defined by `colname`).
  # requires: seqPattern
  
  motifHits <- peaks %>% 
    GRanges %>% 
    computeSummits() %>% 
    resize(resize_bp, "center") %>% 
    Biostrings::getSeq(genome, .) %>% 
    seqPattern::motifScanHits(., motifPWM = pwm, minScore = minScore)
  
  motifHits_count <- motifHits %>% 
    dplyr::count(sequence) %>% 
    dplyr::rename(!!colname := n)
 
  out_column <- rlang::sym(colname) # ORIGINAL
  
  peaks_countAnno <- peaks %>% 
    data.frame %>% 
    dplyr::mutate(seqID = seq(1,nrow(.))) %>% 
    dplyr::left_join(., motifHits_count, by = c("seqID" = "sequence")) %>% 
    dplyr::select(-seqID) %>% #dplyr::select(colname) 
    dplyr::mutate(!!out_column := ifelse(is.na(!!out_column), 0, !!out_column)) %>% # ORIGIINAL LINE
    dplyr::mutate(!!out_column := factor(!!out_column)) #%>% dplyr::select(colname)
  return(peaks_countAnno)
}

motifToPWM <- function(motif){
  # helper function for taking input from motifStack::importMatrix() and output pwm. 
  # Useful for creating input to countMotifMatches.
  motif@mat %>% 
    PWMEnrich::toPWM() %>% 
    .@pwm
}











