#' Example ChIP-seq peaks
#'
#' 10 ChIP-seq peaks from GSE141738
#'
#' A small number of transcription factor ChIP-seq peaks as a GRanges object,
#' taken from [GSE141738](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141738)
#'
#' @source https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141738&format=file&file=GSE141738%5Funion%5Fchip%5Fpeaks%5Fannotated%2Ecsv%2Egz
"example_peaks"

#' Annotated Transcription Factor ChIP-seq summits
#'
#' ChIP-seq summit positions on Drosophila melanogaster chromosome 3 for the
#' transcription factor E93 using a union set of peaks from third-instar larval
#' wings ("Early") and 24 hour Pupal ("Late") wings.
#'
#' E93 is a transcription factor normally present only in Late wings. An
#' experimental perturbation precociously expressed E93 during Early stages.
#' Binding profiles between Early and Late E93 were compared. Peaks are
#' annotated as whether they are bound in Early wings only ("ectopic"), both
#' Early and Late wings ("entopic"), or only bound in Late
#' wings ("orphan").
#'
#' DNA elements can be made "open" or "closed" in response to binding of
#' transcription factors like E93. Accessibility of E93 binding sites before and
#' after E93 expression was measured using FAIRE-seq. ChIP peaks are annotated
#' by how their accessibility changes in response to E93 binding . Peaks can
#' become more open ("Increasing"), more closed ("Decreasing"), or unchanged in
#' accessibility ("Static"). These experiments demonstrate a causal relationship
#' between E93 binding and both opening and closing of DNA elements.
#'
#' @format A GRanges object of ChIP summit position with 2 metadata columns
#'  \describe{
#' \item{peak_binding_description}{Binding profiles between Early and Late E93
#' were compared. Peaks are annotated as whether they are bound in Early wings
#' only ("ectopic"), both Early and Late wings ("entopic"), or only bound in
#' Late wings ("orphan").}
#'   \item{e93_sensitive_behavior}{change in chromatin accessibility in response to E93 binding: Increasing, Decreasing, or Static}
#'  }
#'
#'
#' @source https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141738&format=file&file=GSE141738%5Funion%5Fchip%5Fpeaks%5Fannotated%2Ecsv%2Egz
"example_chip_summits"

#' Example runDreme() output
#'
#' Result when running dreme using 100bp window around example_chip_summits
#' using "Decreasing" sites as foreground, and "Static" sites as background.
#'
#' @format a runDreme results data.frame
#'
"example_dreme"

#' Example runDreme() output after passing to runTomTom()
#'
#' Result when running `runTomTom(example_dreme)` using FlyFactorSurvey as database.
#'
#' @format a runDreme results data.frame with runTomTom results columns
#'
"example_dreme_tomtom"

#' runDreme() output for example_chip_summits split by binding description
#'
#' See vignette("integrative_analysis", package = "dremeR") for details
#'
#'
#' @format a runDreme results data.frame
"example_dreme_by_binding"

#' runDreme() output for example_chip_summits split by response behavior
#'
#' See `vignette("integrative_analysis", package = "dremeR")` for details
#'
#' @format a runDreme results data.frame
#'
"example_dreme_by_sens_vs_static"

#' Example runTomTom() output
#'
#' Result when running `runTomTom(example_dreme$motif)` using FlyFactorSurvey as database
#'
#' @format a data.frame
#'
"example_tomtom"

#' Example runAme() output
#'
#' Result when running AME using 100bp window around `example_chip_summits` for
#' "Increasing" and "Decreasing" sites, using "Static" as background.
#'
#' @format A list object of AME results data.frames
#' \describe{
#'  \item{Increasing}{`runAme()` Results object for Increasing sites vs Static sites}
#'  \item{Decreasing}{`runAme()` Results object for Decreasing sites vs Static sites}
#' }
#'
#' @examples
#' # Data can be combined into 1 large data.frame using:
#' # where the "behavior" column will hold the "Increasing"/"Decreasing" information
#' dplyr::bind_rows(example_ame, .id = "behavior")
"example_ame"

#' runAme() output for example_chip_summits split by binding description
#'
#' AME was run for "ectopic", "entopic", and "orphan" sites using shuffled background.
#'
#' see `vignette("integrative_analysis", package = "dremeR")` for details.
#'
#' @format a list of runAme() results data.frames
#' @examples
#' # Data can be combined into 1 large data.frame using:
#' dplyr::bind_rows(example_ame_large, .id = "binding_type")
"example_ame_large"

#' Example runFimo() output
#'
#' Run using 100bp windows around `example_chip_summits`, using E93 motif as database.
#'
#' @format A GRanges object of E93 motif positions within 100bp windows of `example_chip_summits`
#'
"example_fimo"
