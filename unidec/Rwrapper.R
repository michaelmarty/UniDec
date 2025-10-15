setClass("DeconResult", slots = list(
  mass = "numeric",
  intensity = "numeric"
  # Add other slots as needed
))

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
    # Print file name for debugging
    cat("Script file is", file, "\n")

    return(dirname(file))
  }
  # Otherwise look for arg equal to -f exactly
  match <- which(args == "-f")
  # match <- grep("-f", args)
  # Print matc

  if (length(match)) {
    file <- args[match + 1]
    # Print file name for debugging
    cat("Script file is", file, "\n")

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
  ## Custom search for DLL as in Python (add details as needed)
  candidate <- file.path(guesspath, dllname)
  if (file.exists(candidate)) {
    candidate
  } else {
    # Try searching in unidec/bin as a fallback
    fallback_candidate <- file.path("unidec", "bin", dllname)
    if (file.exists(fallback_candidate)) fallback_candidate else NULL
  }
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

# Register the custom plugin, WILL NEED TO ADJUST THIS TO YOUR PATH (Fix that in future)
Rcpp::registerPlugin("unideclib", function() {
  list(
    env = list(
      PKG_LIBS = "-LC:/Python/UniDecDev/unidec/bin -lunideclib"
    )
  )
})

library(Rcpp)
# Set the path to your C++ file as in the src directory of the current path
sourceCpp("src/rcpp_unidec.cpp") # path to your cpp file
# cat("Test0")
# loadModule("unidec_mod", TRUE)
#
# UDConfig <- function(){
#   cfg <- new(Config)
#   cfg$SetDefaultConfig()
#
#   cfgptr <- create_config()
#
#   return(cfg)
# }
#
# cfg <- UDConfig()
# inp <- new(InputVec)
# cat("Test1")

cfg <- create_config()
inp <- create_inputvec()

print(typeof(cfg))
#
# cat("Test2")

# Data file path "C:\Data\UniDecTest\BSA.txt"
datapath <- file.path("C:", "Data", "UniDecTest", "BSA.txt")
print(datapath)
# Open data file in R as x y columns
data <- read.table(datapath, header = TRUE)
x <- data[[1]]
y <- data[[2]]

# Print first few rows of data for verification
# print(head(data))

#
# x <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
# y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

set_input_data(inp, cfg, x, y)

# cat("Test4")

decon <- run_unidec(cfg, inp)

# Plot decon
x <- decon@mass
y <- decon@intensity

plot(x, y, type = "l", xlab = "Mass", ylab = "Intensity", main = "Deconvolution Result")





