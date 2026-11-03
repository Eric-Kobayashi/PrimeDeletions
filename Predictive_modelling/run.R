# Run from any working directory: pixi run reproduce --data-dir PATH [--model PATH].
args <- commandArgs(trailingOnly=TRUE)
usage <- paste("Usage: pixi run reproduce --data-dir PATH [--model SAVED_MODEL.json]",
               "[--reference ARCHIVED_PREDICTIONS.csv] [--output PATH] [--figures-only]")
if ("--help" %in% args || "-h" %in% args) { cat(usage,"\n"); quit(status=0) }
get_arg <- function(name, default=NULL) {
  idx <- which(args == name)
  if (length(idx) == 0) return(default)
  if (length(idx) != 1 || idx == length(args) || startsWith(args[idx+1],"--"))
    stop("Missing or repeated value for ",name)
  args[idx+1]
}
allowed <- c("--data-dir","--output","--model","--reference","--figures-only")
unknown <- args[startsWith(args,"--") & !args %in% allowed]
if (length(unknown)) stop("Unknown option: ",paste(unknown,collapse=", "))
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
script_dir <- dirname(normalizePath(script_path, mustWork=TRUE))
output <- get_arg("--output", "output")
dir.create(output, recursive=TRUE, showWarnings=FALSE)
output <- normalizePath(output, mustWork=TRUE)
python <- Sys.getenv("RETICULATE_PYTHON", unset=Sys.which("python"))
if (!nzchar(python)) stop("Run inside the Pixi environment to select Python.")
Sys.setenv(RETICULATE_PYTHON=python, RETICULATE_CONDA=Sys.which("micromamba"),
           MAMBA_ROOT_PREFIX=file.path(output,".mamba"),
           XDG_CACHE_HOME=file.path(output,".cache"),
           MPLBACKEND="Agg", MPLCONFIGDIR=file.path(output,".matplotlib"))
if (!"--figures-only" %in% args) {
  data_dir <- get_arg("--data-dir")
  if (is.null(data_dir)) stop(usage)
  model <- get_arg("--model")
  mode <- if (is.null(model)) "train" else "prepare"
  python_args <- c(file.path(script_dir,"model.py"),mode,"--data-dir",data_dir,"--output",output)
  if (!is.null(model)) python_args <- c(python_args,"--model",model)
  reference <- get_arg("--reference")
  if (!is.null(reference)) python_args <- c(python_args,"--reference",reference)
  status <- system2(python,args=vapply(python_args,shQuote,character(1)))
  if (status != 0) stop("Model preparation failed; inspect the metrics before plotting.")
}
reticulate::use_python(python, required=TRUE)
module <- reticulate::import_from_path("model",path=script_dir,convert=FALSE)
context <- module$load_plot_inputs(output)
for (panel in c("i","j","k")) source(file.path(script_dir,paste0("figure2",panel,".R")))
plot_figure2i(file.path(output,"test_predictions.csv"),file.path(output,"figure2i.pdf"))
plot_figure2j(context,file.path(output,"figure2j.pdf"))
plot_figure2k(context,file.path(output,"figure2k.pdf"))
cat("Wrote Figure 2i, 2j and 2k to ",output,"\n")
