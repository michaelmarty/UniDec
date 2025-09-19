# Platform-specific DLL file extension
dllname <- switch(Sys.info()[["sysname"]],
                 "Windows" = "unideclib.dll",
                 "Linux"   = "unideclib.so",
                 "Darwin"  = "unideclib.dylib")

# Print DLL name
cat("DLL name is", dllname, "\n")

# Set working directory to script's directory
get_local_path <- function() {
  # RStudio case
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(path)) return(dirname(path))
  }
  # Command line Rscript
  args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, args)
  if (length(match)) {
    file <- sub(needle, "", args[match])
    return(dirname(file))
  }
  # Fallback to working directory
  return(getwd())
}
current_path <- get_local_path()
setwd(current_path)
cat("Current path is", current_path, "\n")


# Function to search for DLL, if analogous
start_at_iso <- function(dllname, guesspath = "bin") {
  # Custom search for DLL as in Python (add details as needed)
  candidate <- file.path(guesspath, dllname)
  if (file.exists(candidate)) candidate else NULL
}

# Locate the DLL
default_dll_path <- start_at_iso(dllname)
if (is.null(default_dll_path)) {
  cat("DLL not found anywhere\n")
} else {
  cat("Using DLL at", default_dll_path, "\n")
}

# Check if DLL path is valid
if (!file.exists(default_dll_path)) {
  stop("DLL file not found at the specified path: ", default_dll_path)
}

# Add folder with DLL to system PATH
dll_dir <- dirname(normalizePath(default_dll_path))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), dll_dir, sep = .Platform$path.sep))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/rtools45/usr/bin", sep=";"))
cat("Updated system PATH to include:", dll_dir, "\n")

# Load the DLL
dyn.load(default_dll_path)

library(Rcpp)
# Set the path to your C++ file as in the src directory of the current path
sourceCpp("src/rcpp_unidec.cpp") # path to your cpp file
loadModule("unidec_mod", TRUE)

cfg <- new(Config)

inp <- new(Input)
inp$dataMZ <- c(100, 200, 300)
inp$nztab <- c(1L, 2L, 3L)
inp$testmasses <- c(199.5, 201.2)



