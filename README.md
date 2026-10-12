# FuseViewer

Interactive Makie GUI for inspecting time-dependent [FUSE](https://github.com/ProjectTorreyPines/FUSE.jl) data: equilibrium, profiles, magnetics, core sources, core transport and NBI, with a time slider and simulation-vs-experiment overlays.

## Requirements

Julia 1.11 with these packages available in the default `@v1.11` environment:

- FUSE, IMAS, IMASdd
- GLMakie
- Revise

## Data

The viewer loads a directory containing `dd_sim.json` and, optionally, `dd_exp.json`.

## Usage

Place (or clone) this directory inside the shot directory, then:

```bash
./fuseviewer/run                 # data directory = the one containing fuseviewer/
./fuseviewer/run /path/to/shot   # or name it explicitly
```

A splash window opens with the data directory prefilled. Press Enter to load it; the splash is replaced by the full viewer. Julia exits when the window is closed.

For a REPL alongside the window, start Julia with `-i`:

```bash
julia -i --project=@v1.11 fuseviewer/setup.jl
```

```julia
FuseViewer.open_dir("/path/to/another/shot")   # switch shots
FuseViewer.start()                             # back to the splash
FuseViewer.launch(dd_sim, dd_exp)              # view dd objects already in memory
```

From a notebook:

```julia
include("fuseviewer/fuseviewer.jl")
FuseViewer.open_dir("/path/to/shot")                              # GLMakie window
fig = FuseViewer.open_dir("/path/to/shot"; display_window=false)  # inline, static
```

## Files

| File | Purpose |
| --- | --- |
| `fuseviewer.jl` | The `FuseViewer` module |
| `setup.jl` | Loads FUSE/IMAS/GLMakie, tracks `fuseviewer.jl` with Revise, opens the splash |
| `run` | Shell launcher for `setup.jl` |

## Development

`setup.jl` loads `fuseviewer.jl` through Revise, so edits are picked up in a `julia -i` session. Changes inside a panel builder apply on the next plot update; layout changes need the window closed and `launch`/`open_dir` called again.
