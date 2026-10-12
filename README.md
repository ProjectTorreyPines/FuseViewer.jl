# FuseViewer

Interactive Makie GUI for inspecting time-dependent [FUSE](https://github.com/ProjectTorreyPines/FUSE.jl) data: equilibrium, profiles, magnetics, core sources, core transport and NBI, with a time slider and simulation-vs-experiment overlays.

## Installation

FuseViewer (Julia 1.11+) is registered in [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl). If you have FUSE installed you already have that registry; otherwise add it first:

```julia
using Pkg
pkg"registry add https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"
```

Then:

```julia
Pkg.add("FuseViewer")
```

To work on the viewer itself, use `Pkg.develop("FuseViewer")` instead, which clones an editable copy to `~/.julia/dev/FuseViewer`.

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

These open an interactive GLMakie window and return `nothing`, from the REPL or from a Jupyter notebook. To render inline in a notebook instead, ask for the figure. The inline render is static, so the slider and toggles are lost; choose what it shows with `time` and `view`:

```julia
fig = FuseViewer.open_dir("/path/to/shot"; display_window=false)
fig = FuseViewer.launch(dd_sim, dd_exp; display_window=false, time=2.1, view="NBI")
```

The views are Summary, Magnetics, Sources, NBI, EC and Transport.

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

## License

Apache License 2.0. See [LICENSE](LICENSE).
