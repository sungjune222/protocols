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
  root_path <- get_env("ROOT_DIR")
  path <- get_env(key)
  full_path <- file.path(root_path, path)

  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
  }
  full_path
}
