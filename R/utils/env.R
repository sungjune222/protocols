library(dotenv)
library(rprojroot)

.root <- rprojroot::find_root(rprojroot::has_file("pixi.toml"))
.env_path <- file.path(.root, ".env")

get_env <- function(key) {
  dotenv::load_dot_env(file = .env_path)

  path <- Sys.getenv(key, unset = NA)
  if (is.na(path) || path == "") {
    stop(paste("Environment variable not set:", key))
  }
  path
}

get_env_dir <- function(key) {
  path <- get_env(key)

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  path
}
