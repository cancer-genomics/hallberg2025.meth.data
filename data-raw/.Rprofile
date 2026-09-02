## R does not search parent directories for .Rprofile -- without this file,
## Rscript invoked from hallberg2025.meth.data/data-raw/ (as this pipeline's own
## _targets.R header instructs) falls back to ~/.Rprofile, which sets a
## personal library incompatible with this project's renv library and
## crashes on dyn.load() (observed 2026-07-04: rlang.so built for a
## different R minor version). Source the project root's .Rprofile with
## cwd temporarily switched to the root, so its `file.exists("renv/activate.R")`
## check (relative-path, not here()-based) resolves correctly, then restore
## cwd so this pipeline's own here("..","..",...) -relative paths still work.
##
## The source() is guarded because the root .Rprofile is not always there.
## build_public.sh excludes .Rprofile from the public release, so in CI this
## package is checked out beside a project root that has none, and an
## unconditional source() aborts the run. Skipping it there is correct, not
## merely tolerable: the ~/.Rprofile hazard above is a laptop hazard, and the
## CI container has no ~/.Rprofile and sets R_LIBS explicitly (see the env:
## block in check.yml). On the laptop the root .Rprofile does exist and is
## still sourced exactly as before.
local({
  original_wd <- getwd()
  setwd(file.path(original_wd, "..", ".."))
  if (file.exists(".Rprofile")) source(".Rprofile")
  setwd(original_wd)
})
