# One-shot setup script for the FuseViewer.
#
# Run from your own terminal:
#
#     julia --project=@v1.11 fuseviewer/setup.jl
#
# A small splash window opens, prefilled with the data directory. Edit the
# path if you want a different folder containing `dd_sim.json` (and optionally
# `dd_exp.json`) and press Enter; once loaded the splash is replaced by the
# full viewer. Julia exits when you close the window.
#
# For an interactive REPL alongside the window, add -i:
#
#     julia -i --project=@v1.11 fuseviewer/setup.jl
#
# At the REPL:
#
#     FuseViewer.open_dir("/path/to/another/shot")    # programmatic switch
#     FuseViewer.start()                              # back to the splash

using Revise

@info "Loading FUSE / IMAS / GLMakie ..."
@time using FUSE
@time using IMAS
@time using IMASdd
@time using GLMakie

const DATA_DIR = dirname(@__DIR__)

@info "Tracking fuseviewer.jl with Revise ..."
includet(joinpath(@__DIR__, "fuseviewer.jl"))

@info "Opening splash (no data loaded yet) ..."
FuseViewer.start(; default_dir=DATA_DIR)

@info """

Splash is open. Type a directory path and press Enter to load it.

  FuseViewer.open_dir("/path/to/shot")    # programmatic load
  FuseViewer.start()                      # back to splash
"""

# Block until the window is closed (no-op when Julia is started with -i).
isinteractive() || FuseViewer.wait_for_close()
