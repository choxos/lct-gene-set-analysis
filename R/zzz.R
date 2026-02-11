.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "lctGSA v", utils::packageVersion("lctGSA"), "\n",
    "Linear Combination Test for Gene Set Analysis\n",
    "Documentation: ?lctGSA or browseVignettes('lctGSA')\n",
    "\nNote: qvalue package (Bioconductor) is recommended for FDR correction.\n",
    "Install with: BiocManager::install('qvalue')"
  )
}
