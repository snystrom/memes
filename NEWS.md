# memes 1.7.1
* Fixed typos describing the `parse_genomic_coord` parameter in `runFimo` and `runMeme` documentation & error messages for clarity. (Nice catch, Soren Heidelbach!)
* Added CITATION file
# memes 1.7.0
* Fixed a bug in `parseStreme` causing motif ids and list object names to be out of sync.
* Changed `parseStreme` and `runStreme` output motif ids slightly to always pad leading 0's to the numeric ids, consistent with how other MEME Suite tools behave. For instance, a previous motif may be named "m1-ATGC", will now be named "m01-ATGC". This may cause slight code breakage if filtering on specific hard-coded ids. This change also ensures correct sorting order of motif outputs consistent with the behavior of the rest of the package.
# memes 1.5.2
* added a warning in `runFimo` in case `bfile` argument matches an existing local file.
# memes 1.4.1
* added an informative error message during `importMeme` when `parse_genomic_coords` fails with custom fasta import.
# memes 1.2.5
* fixed error in `ame_compare_heatmap_methods` that triggers when plotting without providing a `group` argument.
# memes 1.2.4
* updated NAMESPACE to fix R CMD CHECK note.
# memes 1.2.3
* fixed a bug in `importTomTomXML` where `tomtom` list column would contain missing data if `tomtom` was run using multiple database sources as input.
# memes 1.2.2
* fixed a bug in `runStreme` causing failures on data import for STREME version >= 5.4.1
# memes 1.2.1
* fixed a bug in `runStreme` causing failures when using BStringSetLists as input
# memes 1.1.3
* fixed a bug in `runTomTom` where setting `norc = TRUE` failed on data import
# memes 1.1.1
* `runFimo` now returns `NULL` and prints a message when `text = FALSE` and FIMO detects no matches instead of throwing a cryptic error message
# memes 0.99.11
* Add support for STREME with `runStreme()`. STREME will supercede DREME in a future MEME Suite release.
# memes 0.99.10
* Fixed a bug where paths weren't correctly expanded when used as `database` entry under certain conditions
# memes 0.99.8
* Removed inline `r` call in integrative_analysis vignette to fix issue on bioc build machine
# memes 0.99.7
* Version bump to force pkg rebuild
# memes 0.99.6
* added list S3 method for `plot_sequence_heatmap` so now named lists of sequences can be passed natively to this function.
  * updated ChIP-seq vignette to demonstrate this

# memes 0.99.5
* added `plot_sequence_heatmap` for making heatmaps of sequence lists
* Added significantly more explanation to the ChIP-seq vignette
* renamed `ame_plot_heatmap` -> `plot_ame_heatmap` for consistency 

# memes 0.1.2
* `runFimo()` `skip_matched_sequence` default is now `FALSE`. Set this to `TRUE` if fimo takes a long time to run, then use `add_sequence()` to add it back if needed.
* `runTomTom()` `dist` default is now `ed` (changed from `pearson`).

# memes 0.1.0
* Removed `as_universalmotif_df()`, `as_universalmotif()`, and `update_motifs()`.
  * These functions are replaced by `universalmotif::to_df()`, `universalmotif::to_list()`, and `universalmotif::update_motifs()`
* `runDreme` and `runTomTom` results are now returned in `universalmotif_df` format (behaves just like a data.frame)
  * The `motif` column of `universalmotif_df` objects can no longer be called directly in `universalmotif` operations like `view_motifs(df$motif)`. Use `to_list()` for this behavior instead.
  * To support this change, the `pvalue`, `evalue`, and `qvalue` columns are renamed `pval`, `eval`, and `qval`. The same is true for tomtom output columns `match_pvalue` -> `match_pval`, `best_match_pvalue` -> `best_match_pval`, etc.
  * Updated example datasets to use `unviversalmotif_df` type
* `ame_plot_heatmap` ranking issue is resolved, plots now sort correctly
* Added `remove_duplicate_motifs` and `has_duplicate_motifs` for detecting and removing duplicated matrices in a universalmotif object or data.frame
* Overhauled the Tidying Motifs vignette for more extensive EDA and a demo of deduplication
  * Updated the `flyFactorSurvey_cleaned.meme` example database to reflect new changes to the vignette

