# FuseViewer

Interactive Makie GUI for inspecting time-dependent [FUSE](https://github.com/ProjectTorreyPines/FUSE.jl) data: equilibrium, profiles, magnetics, core sources, core transport and NBI, with a time slider and simulation-vs-experiment overlays.

## Installation

FuseViewer is an unregistered Julia package (Julia 1.11+). Add it to an environment that already has IMAS and IMASdd, such as the one you run FUSE from:

```julia
using Pkg
Pkg.develop(path="/path/to/fuseviewer")
```

Its dependencies are GLMakie, Makie, IMAS and IMASdd.

## Usage

The viewer loads a directory containing `dd_sim.json` and, optionally, `dd_exp.json`.

```julia
using FuseViewer

FuseViewer.open_dir("/path/to/shot")     # load the JSONs and open the viewer
FuseViewer.launch(dd_sim, dd_exp)        # view dd objects already in memory
FuseViewer.start()                       # splash window with a directory box
```

`open_dir` also binds the loaded data to `dd_sim` and `dd_exp` in `Main`.

In a notebook, ask for the figure instead of a window. The inline render is static, so the slider and toggles are lost:

```julia
fig = FuseViewer.open_dir("/path/to/shot"; display_window=false)
```

From a non-interactive script, block until the window is closed:

```bash
julia -e 'using FuseViewer; FuseViewer.start(); FuseViewer.wait_for_close()'
```

## Development

Load Revise before FuseViewer and edits to `src/FuseViewer.jl` are picked up in the running session. Changes inside a panel builder apply on the next plot update; layout changes need the window closed and `launch`/`open_dir` called again.

```julia
using Revise
using FuseViewer
```

Run the tests with `Pkg.test("FuseViewer")`.
