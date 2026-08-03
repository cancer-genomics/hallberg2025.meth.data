## R does not search parent directories for .Rprofile -- without this file,
## Rscript invoked from hallberg2025.meth.data/data-raw/ (as this pipeline's own
## _targets.R header instructs) falls back to ~/.Rprofile, which sets a
## personal library incompatible with this project's renv library and
## crashes on dyn.load() (observed 2026-07-04: rlang.so built for a
## different R minor version). Source the project root's .Rprofile with
## cwd temporarily switched to the root, so its `file.exists("renv/activate.R")`
## check (relative-path, not here()-based) resolves correctly, then restore
## cwd so this pipeline's own here("..","..",...) -relative paths still work.
local({
  original_wd <- getwd()
  setwd(file.path(original_wd, "..", ".."))
  source(".Rprofile")
  setwd(original_wd)
})
