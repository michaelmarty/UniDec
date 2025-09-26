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

# Register the custom plugin
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
cat("Test0")
loadModule("unidec_mod", TRUE)
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
cat("Test1")

cfg <- create_config()
inp <- create_inputvec()

print(typeof(cfg))

cat("Test2")
x <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

set_input_data(inp, cfg, x, y)

cat("Test4")

test <- run_unidec(cfg, inp)

print("Test5")


