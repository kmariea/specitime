args <- commandArgs(trailingOnly = FALSE)

script_path <- sub(
  "^--file=",
  "",
  args[grepl("^--file=", args)]
)

app_dir <- normalizePath(
  dirname(script_path),
  winslash = "/",
  mustWork = TRUE
)

shiny::runApp(
  appDir = app_dir,
  launch.browser = TRUE
)