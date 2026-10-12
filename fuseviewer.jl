# FuseViewer — interactive Makie GUI for inspecting time-dependent FUSE data.
#
# Workflow (see setup.jl):
#
#   Your terminal:   julia --project=@v1.11 fuseviewer/setup.jl
#   At the prompt:   FuseViewer.launch(dd_sim, dd_exp)
#
# From a notebook (IJulia/Pluto) there is no splash to click through — name
# the directory directly:
#
#   include("fuseviewer/fuseviewer.jl")
#   FuseViewer.open_dir("/path/to/shot")                       # GLMakie window
#   fig = FuseViewer.open_dir("/path/to/shot"; display_window=false)  # inline
#
# `display_window=false` skips `display(fig)` (and the closeall that precedes
# it) and hands back the Figure, so the notebook renders it with whatever
# backend is active. Note that a static inline render loses the slider and
# toggles; keep the default GLMakie window for interactive work.
#
# There is no dependency on a native file dialog: directories are typed into
# a Textbox or passed as an argument, so the module loads anywhere Makie does.
#
# Revise picks up edits to this file. After a layout change, close the
# window and call `launch` again. After an in-place change inside a panel
# builder, no relaunch is needed — the new code runs on the next plot
# update.

module FuseViewer

using GLMakie
using Makie: Point2f
using IMAS
using IMASdd: json2imas

export launch, open_dir, start, wait_for_close

# Tracks the most-recently displayed screen so wait_for_close() can follow
# the splash → viewer transition that open_dir() performs.
const _active_fig    = Ref{Any}(nothing)
const _active_screen = Ref{Any}(nothing)

# HCD source categories: identifier.index values and display properties.
# Matches IMASdd.index_hcd_sources = [2, 3, 4, 5, 14].
const _HCD_CATEGORIES = (
    (idx=2,  name="NBI",    color=:tomato),
    (idx=3,  name="EC",     color=:cornflowerblue),
    (idx=4,  name="LH",     color=:orchid),
    (idx=5,  name="IC",     color=:mediumseagreen),
    (idx=14, name="pellet", color=:darkorange),
)
const _HCD_INDEXES = Set(c.idx for c in _HCD_CATEGORIES)

# ---------------------------------------------------------------------------
# Data extraction helpers (pure functions of dd + time)
# ---------------------------------------------------------------------------

"Build a single Point2f vector that draws every flux contour for `eqt`,
using NaN points as polyline separators so we can render with one `lines!`."
function flux_surface_points(eqt::IMAS.equilibrium__time_slice;
                             nlevels_in::Int=11, nlevels_out::Int=4)
    pts = Point2f[]
    isempty(eqt.profiles_2d) && return pts

    fw = try
        IMAS.first_wall(IMAS.top_dd(eqt).wall)
    catch
        (r=Float64[], z=Float64[])
    end

    ψ1d = eqt.profiles_1d.psi
    ψaxis, ψbnd = ψ1d[1], ψ1d[end]
    Δψ = ψbnd - ψaxis
    levels_in = range(ψaxis, ψbnd, nlevels_in)
    levels_out = ψbnd .+ (Δψ / nlevels_in) .* (1:nlevels_out)
    levels = vcat(collect(levels_in), collect(levels_out))

    for ψlev in levels
        for fs in IMAS.flux_surface(eqt, ψlev, :any, fw.r, fw.z)
            isempty(fs.r) && continue
            for (r, z) in zip(fs.r, fs.z)
                push!(pts, Point2f(r, z))
            end
            push!(pts, Point2f(NaN, NaN))
        end
    end
    return pts
end

function boundary_points(eqt::IMAS.equilibrium__time_slice)
    bnd = eqt.boundary.outline
    if ismissing(bnd, :r) || ismissing(bnd, :z) || isempty(bnd.r)
        return Point2f[]
    end
    return Point2f.(bnd.r, bnd.z)
end

function maxis_point(eqt::IMAS.equilibrium__time_slice)
    ma = eqt.global_quantities.magnetic_axis
    return Point2f(ma.r, ma.z)
end

function xpoints(eqt::IMAS.equilibrium__time_slice)
    pts = Point2f[]
    for xp in eqt.boundary.x_point
        push!(pts, Point2f(xp.r, xp.z))
    end
    return pts
end

function wall_points(dd::IMAS.dd)
    fw = IMAS.first_wall(dd.wall)
    isempty(fw.r) && return Point2f[]
    return Point2f.(fw.r, fw.z)
end

"Pull a 1D profile vector from cp1d at `time0` and return Point2f(rho, value*scale)."
function profile_points(dd::IMAS.dd, time0::Real, accessor::Function;
                        scale::Real=1.0)
    try
        cp1d = dd.core_profiles.profiles_1d[Float64(time0)]
        rho = cp1d.grid.rho_tor_norm
        y = accessor(cp1d)
        (y === nothing || isempty(y)) && return Point2f[]
        return Point2f.(rho, y .* scale)
    catch
        return Point2f[]
    end
end

"Pull a 1D profile vector from eqt.profiles_1d at `time0` (rho_tor_norm vs accessor(eqt))."
function eq_profile_points(dd::IMAS.dd, time0::Real, accessor::Function;
                           scale::Real=1.0)
    try
        eqt = dd.equilibrium.time_slice[Float64(time0)]
        rho = eqt.profiles_1d.rho_tor_norm
        y = accessor(eqt)
        (y === nothing || isempty(y)) && return Point2f[]
        return Point2f.(rho, y .* scale)
    catch
        return Point2f[]
    end
end

"""
Sample a profile line across `times` (downsampled by `sample_step`) and return
robust `(lo, hi)` y-limits with a 5% margin. Used to pin axes once at launch
so the slider doesn't bounce the y-range from frame to frame.

Diagnostic scatter is excluded — the limits track the simulation/experiment
line only.
"""
function global_ylims(accessor::Function, dds, times;
                      from_eq::Bool=false, scale::Real=1.0,
                      sample_step::Int=4, margin::Real=0.05)
    pts_fn = from_eq ? eq_profile_points : profile_points
    lo = Inf; hi = -Inf
    for dd in dds
        dd === nothing && continue
        for t in times[1:sample_step:end]
            p = pts_fn(dd, t, accessor; scale)
            for pt in p
                y = pt[2]
                isnan(y) && continue
                lo = min(lo, y); hi = max(hi, y)
            end
        end
    end
    (isfinite(lo) && isfinite(hi)) || return nothing
    rng = hi - lo
    rng <= 0 && (rng = abs(hi) + 1.0)
    return (lo - margin * rng, hi + margin * rng)
end

"""
Pull a diagnostic measurement (Thomson scattering or CER) at `time0` and
return scatter positions ready for `scatter!`. Wraps `IMAS.getdata` and
returns empty on any failure.

`what` is one of `:t_e`, `:n_e` (Thomson), or `:t_i`, `:n_imp`, `:ω_tor` (CER).
"""
function diag_points(dd::IMAS.dd, time0::Real, what::Symbol;
                     time_averaging::Float64=0.05, scale::Real=1.0)
    try
        result = IMAS.getdata(Val(what), dd, Float64(time0), time_averaging)
        rho = collect(result.rho)
        vals = Float64[d.val for d in result.data]
        mask = .!isnan.(vals) .& .!isnan.(rho)
        return Point2f.(rho[mask], vals[mask] .* scale)
    catch
        return Point2f[]
    end
end

"Look up the impurity charge from cp1d.ion[2].element[1].z_n; default 1.0 if missing."
function impurity_z(dd::IMAS.dd, time0::Real)
    try
        cp1d = dd.core_profiles.profiles_1d[Float64(time0)]
        return Float64(cp1d.ion[2].element[1].z_n)
    catch
        return 1.0
    end
end

# 0D time traces (fetched once; not time-dependent on the slider)
function trace_ip(dd::IMAS.dd)
    try
        return dd.core_profiles.time, dd.core_profiles.global_quantities.ip ./ 1e6
    catch
        return Float64[], Float64[]
    end
end

function trace_wmhd(dd::IMAS.dd)
    try
        t = dd.equilibrium.time
        w = [eqt.global_quantities.energy_mhd / 1e6 for eqt in dd.equilibrium.time_slice]
        return t, w
    catch
        return Float64[], Float64[]
    end
end

function trace_ne_line(dd::IMAS.dd)
    try
        t = dd.summary.time
        v = dd.summary.line_average.n_e.value
        if !isempty(t) && !isempty(v) && length(t) == length(v)
            return t, v ./ 1e19
        end
    catch
    end
    # Fall back to interferometer: average across channels of n_e_line_average.
    try
        chs = dd.interferometer.channel
        isempty(chs) && return Float64[], Float64[]
        t = chs[1].n_e_line_average.time
        acc = zeros(Float64, length(t))
        cnt = 0
        for ch in chs
            try
                v = ch.n_e_line_average.data
                length(v) == length(t) || continue
                acc .+= v
                cnt += 1
            catch
            end
        end
        cnt == 0 && return Float64[], Float64[]
        return t, (acc ./ cnt) ./ 1e19
    catch
        return Float64[], Float64[]
    end
end

# ---------------------------------------------------------------------------
# Panel builders
# ---------------------------------------------------------------------------

const COL_SIM = :steelblue
const COL_EXP = :black
const COL_DIAG = :firebrick

function build_eq_panel!(fig_pos, t_obs, dd_sim, dd_exp, show_sim, show_exp)
    ax = Axis(fig_pos;
              aspect=DataAspect(),
              xlabel="R [m]", ylabel="Z [m]",
              title=lift(t -> "Equilibrium  —  t = $(round(t; digits=4)) s", t_obs))

    wall_sim = wall_points(dd_sim)
    isempty(wall_sim) || lines!(ax, wall_sim; color=:gray30, linewidth=1.5)

    fs_sim = lift(t -> flux_surface_points(dd_sim.equilibrium.time_slice[Float64(t)]), t_obs)
    bnd_sim = lift(t -> boundary_points(dd_sim.equilibrium.time_slice[Float64(t)]), t_obs)
    ma_sim = lift(t -> [maxis_point(dd_sim.equilibrium.time_slice[Float64(t)])], t_obs)
    xp_sim = lift(t -> xpoints(dd_sim.equilibrium.time_slice[Float64(t)]), t_obs)

    lines!(ax, fs_sim;  color=COL_SIM, linewidth=0.6, visible=show_sim)
    lines!(ax, bnd_sim; color=COL_SIM, linewidth=2.0, label="sim", visible=show_sim)
    scatter!(ax, ma_sim; color=COL_SIM, marker=:cross, markersize=12, visible=show_sim)
    scatter!(ax, xp_sim; color=COL_SIM, marker=:xcross, markersize=14, visible=show_sim)

    if dd_exp !== nothing
        fs_exp = lift(t -> flux_surface_points(dd_exp.equilibrium.time_slice[Float64(t)]), t_obs)
        bnd_exp = lift(t -> boundary_points(dd_exp.equilibrium.time_slice[Float64(t)]), t_obs)
        lines!(ax, fs_exp;  color=COL_EXP, linewidth=0.5, linestyle=:dot, visible=show_exp)
        lines!(ax, bnd_exp; color=COL_EXP, linewidth=1.5, label="exp", visible=show_exp)
    end

    axislegend(ax; position=:rt, framevisible=false)
    return ax
end

"""
Generic profile panel: line for sim and exp from a `cp1d -> y` accessor (or
an `eqt -> y` accessor when `from_eq=true`), with optional diagnostic scatter
overlays sourced from `dd_exp`.

`diags` is a vector of NamedTuples like
    (what=:t_e, label="Thomson", color=:firebrick, marker=:circle, norm=nothing)
where `what` is passed to `IMAS.getdata` and `norm`, when not `nothing`, is
a `(dd, t) -> Float64` callback for a per-time multiplicative factor
(e.g. impurity Z for `:n_imp` -> ne).
"""
function build_profile_panel!(fig_pos, title, ylabel,
                              t_obs, dd_sim, dd_exp,
                              show_sim, show_exp,
                              accessor;
                              scale::Real=1.0,
                              from_eq::Bool=false,
                              diags=(),
                              diag_dt::Float64=0.05,
                              ylims=nothing)
    ax = Axis(fig_pos; xlabel="ρ_tor_norm", ylabel=ylabel, title=title)

    pts_fn = from_eq ? eq_profile_points : profile_points

    sim_pts = lift(t -> pts_fn(dd_sim, t, accessor; scale), t_obs)
    lines!(ax, sim_pts; color=COL_SIM, linewidth=2.0, label="sim", visible=show_sim)

    if dd_exp !== nothing
        exp_pts = lift(t -> pts_fn(dd_exp, t, accessor; scale), t_obs)
        lines!(ax, exp_pts; color=COL_EXP, linewidth=2.0,
               linestyle=:dash, label="exp", visible=show_exp)

        for d in diags
            spec = (
                what    = d.what,
                label   = get(d, :label, string(d.what)),
                color   = get(d, :color, COL_DIAG),
                marker  = get(d, :marker, :circle),
                norm    = get(d, :norm, nothing),
            )
            diag_obs = lift(t_obs) do t
                norm = spec.norm === nothing ? 1.0 : try
                    Float64(spec.norm(dd_exp, t))
                catch
                    1.0
                end
                diag_points(dd_exp, t, spec.what;
                            time_averaging=diag_dt, scale=scale * norm)
            end
            scatter!(ax, diag_obs; color=spec.color, markersize=5,
                     marker=spec.marker, strokewidth=0,
                     label=spec.label, visible=show_exp,
                     inspectable=false)
        end
    end

    if ylims !== nothing
        ylims!(ax, ylims...)
    end
    axislegend(ax; position=:rt, framevisible=false)
    return ax
end

function build_trace_panel!(fig_pos, title, ylabel, t_obs,
                            dd_sim, dd_exp, show_sim, show_exp, fetch::Function)
    ax = Axis(fig_pos; xlabel="time [s]", ylabel=ylabel, title=title)

    n_labeled = 0
    t_s, y_s = fetch(dd_sim)
    if !isempty(t_s)
        lines!(ax, t_s, y_s; color=COL_SIM, linewidth=1.5, label="sim", visible=show_sim)
        n_labeled += 1
    end
    if dd_exp !== nothing
        t_e, y_e = fetch(dd_exp)
        if !isempty(t_e)
            lines!(ax, t_e, y_e; color=COL_EXP, linewidth=1.5,
                   linestyle=:dash, label="exp", visible=show_exp)
            n_labeled += 1
        end
    end

    vlines!(ax, t_obs; color=:gray, linewidth=1.0)
    if n_labeled == 0
        text!(ax, 0.5, 0.5; text="(no data)", space=:relative,
              align=(:center, :center), color=:gray)
    else
        axislegend(ax; position=:rt, framevisible=false)
    end
    return ax
end

# ---------------------------------------------------------------------------
# Profile accessors
# ---------------------------------------------------------------------------

# core_profiles
acc_te(cp1d)       = cp1d.electrons.temperature
acc_ti(cp1d)       = cp1d.t_i_average
acc_ne(cp1d)       = cp1d.electrons.density_thermal
acc_jtor(cp1d)     = cp1d.j_tor
function acc_omega(cp1d)
    isempty(cp1d.ion) && return Float64[]
    try
        return cp1d.ion[1].rotation_frequency_tor
    catch
        return Float64[]
    end
end

# equilibrium.profiles_1d
acc_q(eqt) = eqt.profiles_1d.q

"q95 at `time0`, or NaN if unavailable."
function q95_at(dd::IMAS.dd, time0::Real)
    try
        v = dd.equilibrium.time_slice[Float64(time0)].global_quantities.q_95
        return ismissing(v) ? NaN : Float64(v)
    catch
        return NaN
    end
end

# ---------------------------------------------------------------------------
# Magnetics helpers
# ---------------------------------------------------------------------------

"Vector of (R, Z) Point2f for each magnetic sensor of the given kind."
function magnetics_positions(dd::IMAS.dd, what::Symbol)
    pts = Point2f[]
    if what == :b_field_pol_probe
        for probe in dd.magnetics.b_field_pol_probe
            try
                push!(pts, Point2f(probe.position.r, probe.position.z))
            catch
            end
        end
    elseif what == :flux_loop
        for loop in dd.magnetics.flux_loop
            try
                r = sum(p.r for p in loop.position) / length(loop.position)
                z = sum(p.z for p in loop.position) / length(loop.position)
                push!(pts, Point2f(r, z))
            catch
            end
        end
    end
    return pts
end

"Per-sensor signal value at `time0` (`getavg` weighted-window average)."
function magnetics_values(dd::IMAS.dd, time0::Real, what::Symbol;
                          time_averaging::Float64=0.05)
    try
        result = IMAS.getavg(Val(what), dd, Float64(time0), time_averaging)
        return Float64[d.val for d in result.data]
    catch
        return Float64[]
    end
end

"""
Maximum |sim − exp| residual for a sensor kind across (downsampled) times.
Used to set a symmetric colorrange for the diverging colormap on the
sensor map.
"""
function global_err_range(dd_sim::IMAS.dd, dd_exp::IMAS.dd, what::Symbol;
                          sample_step::Int=8)
    times = collect(dd_sim.equilibrium.time)
    isempty(times) && return 1.0
    step = max(1, length(times) ÷ sample_step)
    abs_max = 0.0
    for t in times[1:step:end]
        m = magnetics_values(dd_exp, t, what)
        s = magnetics_values(dd_sim, t, what)
        n = min(length(m), length(s))
        n == 0 && continue
        for i in 1:n
            d = abs(s[i] - m[i])
            isfinite(d) && (abs_max = max(abs_max, d))
        end
    end
    return abs_max == 0.0 ? 1.0 : abs_max
end

"""
Global (min, max) of sim & exp sensor values for a sensor kind across
(downsampled) times. Used to fix the y-range of the index panels so the
data doesn't leave the view while scanning time. Returns `nothing` if no
finite values are found.
"""
function global_value_range(dd_sim, dd_exp, what::Symbol; sample_step::Int=8)
    lo = Inf; hi = -Inf
    for dd in (dd_sim, dd_exp)
        dd === nothing && continue
        times = try collect(dd.equilibrium.time) catch; Float64[] end
        isempty(times) && continue
        step = max(1, length(times) ÷ sample_step)
        for t in times[1:step:end], v in magnetics_values(dd, t, what)
            isfinite(v) || continue
            lo = min(lo, v); hi = max(hi, v)
        end
    end
    return isfinite(lo) && isfinite(hi) ? (lo, hi) : nothing
end

# ---------------------------------------------------------------------------
# core_sources helpers
# ---------------------------------------------------------------------------

"""
Snapshot of source identifier names at the current global time. Used once at
view-build to determine how many lines to draw and what to label them.
"""
function source_names(dd::IMAS.dd)
    names = String[]
    try
        for source in dd.core_sources.source
            push!(names, try
                string(source.identifier.name)
            catch
                "source_$(length(names) + 1)"
            end)
        end
    catch
    end
    return names
end

"Identifier indices parallel to `source_names`; returns 0 for unknown sources."
function source_indices(dd::IMAS.dd)
    indices = Int[]
    try
        for source in dd.core_sources.source
            push!(indices, try source.identifier.index catch; 0 end)
        end
    catch
    end
    return indices
end

"Profile points (rho, value*scale) for a single core_source at `time0`."
function source_profile_points(source, time0::Real, accessor::Function;
                               scale::Real=1.0)
    try
        cs1d = source.profiles_1d[Float64(time0)]
        rho = cs1d.grid.rho_tor_norm
        y = accessor(cs1d)
        (y === nothing || isempty(y)) && return Point2f[]
        return Point2f.(rho, y .* scale)
    catch
        return Point2f[]
    end
end

"Total source contribution at `time0` via `IMAS.total_sources`."
function total_source_points(dd::IMAS.dd, time0::Real, accessor::Function;
                             scale::Real=1.0)
    try
        total = IMAS.total_sources(dd; time0=Float64(time0))
        rho = total.grid.rho_tor_norm
        y = accessor(total)
        (y === nothing || isempty(y)) && return Point2f[]
        return Point2f.(rho, y .* scale)
    catch
        return Point2f[]
    end
end

"Aggregate of all sources with identifier.index in `include_indexes` at `time0`."
function hcd_category_points(dd::IMAS.dd, time0::Real, accessor::Function,
                              include_indexes::Vector{Int}; scale::Real=1.0)
    try
        agg = IMAS.total_sources(dd; time0=Float64(time0), include_indexes=include_indexes)
        rho = agg.grid.rho_tor_norm
        y = accessor(agg)
        (y === nothing || isempty(y)) && return Point2f[]
        return Point2f.(rho, y .* scale)
    catch
        return Point2f[]
    end
end

# Source channel accessors (cs1d → vector)
acc_src_pe(cs1d)   = cs1d.electrons.energy
acc_src_pi(cs1d)   = cs1d.total_ion_energy
acc_src_part(cs1d) = cs1d.electrons.particles
acc_src_jpar(cs1d) = cs1d.j_parallel
acc_src_mom(cs1d)  = cs1d.momentum_tor

"""
Integrated magnitude of a source's contribution at `time0`. Volume integral
by default; pass `integration=:area` for j_parallel (gives current in A
rather than W/A per m³).
"""
function source_integral(source, time0::Real, accessor::Function;
                         integration::Symbol=:volume)
    try
        cs1d = source.profiles_1d[Float64(time0)]
        y = accessor(cs1d)
        (y === nothing || isempty(y)) && return 0.0
        coord = integration == :area ? cs1d.grid.area : cs1d.grid.volume
        return abs(IMAS.trapz(coord, y))
    catch
        return 0.0
    end
end

"Maximum |integrated source contribution| across (downsampled) `times`."
function source_max_integral(source, accessor::Function, times;
                             integration::Symbol=:volume, sample_step::Int=8)
    isempty(times) && return 0.0
    step = max(1, length(times) ÷ sample_step)
    m = 0.0
    for t in times[1:step:end]
        m = max(m, source_integral(source, t, accessor; integration))
    end
    return m
end

# ---------------------------------------------------------------------------
# core_transport helpers
# ---------------------------------------------------------------------------

# Model identifier indices that are aggregates / placeholders rather than
# physics models; we skip them when iterating per-model.
# Matches IMAS index_2_name__core_transport__model:
#   0=unspecified, 1=combined, 2=transport_solver, 25=unknown
const _CT_SKIP_INDEXES = Set([0, 1, 2, 25])

"Display names of each core_transport model."
function transport_model_names(dd::IMAS.dd)
    names = String[]
    try
        for model in dd.core_transport.model
            push!(names, try
                nm = string(model.identifier.name)
                isempty(nm) ? "model_$(length(names)+1)" : nm
            catch
                "model_$(length(names)+1)"
            end)
        end
    catch
    end
    return names
end

"identifier.index parallel to `transport_model_names`."
function transport_model_indices(dd::IMAS.dd)
    indices = Int[]
    try
        for model in dd.core_transport.model
            push!(indices, try model.identifier.index catch; 0 end)
        end
    catch
    end
    return indices
end

"Flux profile (rho, value*scale) for one core_transport model at `time0`."
function transport_flux_points(model, time0::Real, accessor::Function;
                                scale::Real=1.0)
    try
        ct1d = model.profiles_1d[Float64(time0)]
        rho = ct1d.grid_flux.rho_tor_norm
        y = accessor(ct1d)
        (y === nothing || isempty(y)) && return Point2f[]
        return Point2f.(rho, y .* scale)
    catch
        return Point2f[]
    end
end

"Sum of all transport-model fluxes at `time0` via `IMAS.total_fluxes`."
function total_transport_flux_points(dd::IMAS.dd, time0::Real, accessor::Function;
                                       scale::Real=1.0)
    try
        rho_ref = dd.core_profiles.profiles_1d[Float64(time0)].grid.rho_tor_norm
        total = IMAS.total_fluxes(dd, rho_ref; time0=Float64(time0))
        rho = total.grid_flux.rho_tor_norm
        y = accessor(total)
        (y === nothing || isempty(y)) && return Point2f[]
        return Point2f.(rho, y .* scale)
    catch
        return Point2f[]
    end
end

"Source-implied flux at `time0`: cumulative `inside` / surface area."
function source_implied_flux_points(dd::IMAS.dd, time0::Real, inside_accessor::Function;
                                      scale::Real=1.0)
    try
        total = IMAS.total_sources(dd; time0=Float64(time0))
        rho = total.grid.rho_tor_norm
        surf = total.grid.surface
        inside = inside_accessor(total)
        (inside === nothing || isempty(inside) || isempty(surf)) && return Point2f[]
        flux = similar(inside, Float64)
        for i in eachindex(inside)
            flux[i] = surf[i] > 0 ? inside[i] / surf[i] : NaN
        end
        return Point2f.(rho, flux .* scale)
    catch
        return Point2f[]
    end
end

# Transport flux accessors (ct1d → flux vector)
acc_qe_flux(ct1d)   = ct1d.electrons.energy.flux
acc_qi_flux(ct1d)   = ct1d.total_ion_energy.flux
acc_part_flux(ct1d) = ct1d.electrons.particles.flux
acc_mom_flux(ct1d)  = ct1d.momentum_tor.flux

# Source-implied flux accessors (cs1d → cumulative-inside vector)
acc_qe_inside(cs1d)   = cs1d.electrons.power_inside
acc_qi_inside(cs1d)   = cs1d.total_ion_power_inside
acc_part_inside(cs1d) = cs1d.electrons.particles_inside
acc_mom_inside(cs1d)  = cs1d.torque_tor_inside

# ---------------------------------------------------------------------------
# Layout helper
# ---------------------------------------------------------------------------

"""
Tear down everything inside a GridLayout so a new view can be built into it.

`Base.empty!(::GridLayout)` only clears the layout's bookkeeping — the
underlying Block scenes stay in the figure and overlay the next view. We
have to walk the children and `delete!` each Block (and recurse into
sub-GridLayouts).
"""
function _clear_layout!(gl::GridLayout)
    items = copy(gl.content)
    for item in items
        try
            child = item.content
            if child isa GridLayout
                _clear_layout!(child)
                Makie.GridLayoutBase.remove_from_gridlayout!(item)
            else
                delete!(child)
            end
        catch err
            @warn "FuseViewer: layout cleanup failure" exception=err
        end
    end
    # Reset any Relative()/Fixed() col/row size constraints set by the
    # previous view, so the next view starts with a uniform layout.
    try
        nr, nc = size(gl)
        for c in 1:nc; colsize!(gl, c, Auto()); end
        for r in 1:nr; rowsize!(gl, r, Auto()); end
    catch
    end
end

# ---------------------------------------------------------------------------
# Directory picker / reloader
# ---------------------------------------------------------------------------

"Normalize a user-typed directory path (trims whitespace, expands `~`)."
_clean_dir(s::AbstractString) = abspath(expanduser(strip(String(s))))

"Textbox validator: accept a path only once it names an existing directory."
function _is_loadable_dir(s::AbstractString)
    str = strip(String(s))
    return !isempty(str) && isdir(_clean_dir(str))
end

"""
    start(; default_dir=pwd())

Open a small splash window with a directory entry box, prefilled with
`default_dir`. Type or paste the path of a folder containing `dd_sim.json`
(and optionally `dd_exp.json`) and press Enter; the splash is then replaced
by the full viewer via `open_dir`. The box stays red until the text names an
existing directory.

Nothing here needs a native file dialog, so the splash works over X11
forwarding and from a notebook. From a notebook or the REPL you can skip it
entirely and call `open_dir("/path/to/shot")` directly.
"""
function start(; default_dir::AbstractString=pwd())
    fig = Figure(; size=(720, 250))

    Label(fig[1, 1], "FUSE Viewer";
          font=:bold, fontsize=22, halign=:center, tellwidth=false)
    Label(fig[2, 1], "Directory containing dd_sim.json (and dd_exp.json) — press Enter to load:";
          halign=:center, tellwidth=false)

    dir0 = String(default_dir)
    tb = Textbox(fig[3, 1];
                 stored_string=dir0, displayed_string=dir0,
                 placeholder="/path/to/shot",
                 validator=_is_loadable_dir,
                 width=620, tellwidth=false)

    # Load failures (missing dd_sim.json, unreadable JSON) surface here rather
    # than only in the log, since the splash may be the only thing on screen.
    status = Label(fig[4, 1], ""; halign=:center, tellwidth=false, color=:firebrick)

    on(tb.stored_string) do dir
        isempty(strip(dir)) && return
        status.text[] = ""
        try
            open_dir(dir)
        catch err
            status.text[] = sprint(showerror, err)
            @error "FuseViewer: load failed" exception=(err, catch_backtrace())
        end
    end

    rowgap!(fig.layout, 14)
    _active_screen[] = display(fig)
    _active_fig[] = fig
    return fig
end

"""
    open_dir(dir::AbstractString;
             sim_file="dd_sim.json", exp_file="dd_exp.json", kwargs...)

Load `dd_sim.json` and `dd_exp.json` from `dir`, replace any open viewer
window, and launch a fresh one. Also rebinds `Main.dd_sim` and `Main.dd_exp`
so the new datasets are addressable from the REPL.

`exp_file=nothing` skips the experimental dataset.

Loading 1+ GB JSONs blocks the UI for minutes; expect a frozen window
between click and the new figure appearing.
"""
function open_dir(dir::AbstractString;
                  sim_file::AbstractString="dd_sim.json",
                  exp_file::Union{Nothing,AbstractString}="dd_exp.json",
                  display_window::Bool=true,
                  kwargs...)
    dir = _clean_dir(dir)
    isdir(dir) || error("not a directory: $dir")

    sim_path = joinpath(dir, sim_file)
    isfile(sim_path) || error("missing $sim_file in $dir")

    exp_path = exp_file === nothing ? nothing : joinpath(dir, exp_file)
    if exp_path !== nothing && !isfile(exp_path)
        @warn "FuseViewer: $exp_file not found in $dir — proceeding sim-only"
        exp_path = nothing
    end

    # Only in windowed mode: a notebook asking for the Figure back has no
    # window to replace, and closeall() would kill unrelated GLMakie figures.
    if display_window
        @info "FuseViewer: closing existing windows"
        GLMakie.closeall()
    end

    @info "FuseViewer: loading $sim_path"
    new_sim = json2imas(sim_path)

    new_exp = if exp_path === nothing
        nothing
    else
        @info "FuseViewer: loading $exp_path"
        json2imas(exp_path)
    end

    # Drop old Main bindings before assigning new ones so the GC can free them
    Core.eval(Main, :(dd_sim = $new_sim))
    Core.eval(Main, :(dd_exp = $new_exp))
    GC.gc()

    return launch(new_sim, new_exp; display_window, kwargs...)
end

# ---------------------------------------------------------------------------
# Views — each builds its panels into a host GridLayout (`view_box`).
# Switching views tears down `view_box` and rebuilds.
# ---------------------------------------------------------------------------

function build_summary_view!(view_box::GridLayout,
                             t_obs, dd_sim, dd_exp,
                             show_sim, show_exp,
                             times)
    # equilibrium (left)
    build_eq_panel!(view_box[1, 1], t_obs, dd_sim, dd_exp, show_sim, show_exp)

    # profile grid (right) — pin ylims from the line data across all times
    rgrid = GridLayout(view_box[1, 2])

    dds_for_lim = (dd_sim, dd_exp)
    ylim_te    = global_ylims(acc_te,    dds_for_lim, times; scale=1e-3)
    ylim_ti    = global_ylims(acc_ti,    dds_for_lim, times; scale=1e-3)
    ylim_ne    = global_ylims(acc_ne,    dds_for_lim, times; scale=1e-19)
    ylim_omega = global_ylims(acc_omega, dds_for_lim, times; scale=1e-3)
    ylim_jtor  = global_ylims(acc_jtor,  dds_for_lim, times; scale=1e-6)

    build_profile_panel!(rgrid[1, 1], "Te [keV]", "T",
                         t_obs, dd_sim, dd_exp, show_sim, show_exp,
                         acc_te; scale=1e-3, ylims=ylim_te,
                         diags=((what=:t_e, label="Thomson"),))

    build_profile_panel!(rgrid[1, 2], "Ti [keV]", "T",
                         t_obs, dd_sim, dd_exp, show_sim, show_exp,
                         acc_ti; scale=1e-3, ylims=ylim_ti,
                         diags=((what=:t_i, label="CER"),))

    build_profile_panel!(rgrid[2, 1], "ne [10¹⁹ m⁻³]", "n",
                         t_obs, dd_sim, dd_exp, show_sim, show_exp,
                         acc_ne; scale=1e-19, ylims=ylim_ne,
                         diags=(
                             (what=:n_e, label="Thomson"),
                             (what=:n_imp, label="CER (× Z_imp)",
                              color=:darkorange, marker=:diamond,
                              norm=(dd, t) -> impurity_z(dd, t)),
                         ))

    build_profile_panel!(rgrid[2, 2], "ω_tor [krad/s]", "ω",
                         t_obs, dd_sim, dd_exp, show_sim, show_exp,
                         acc_omega; scale=1e-3, ylims=ylim_omega,
                         diags=((what=:ω_tor, label="CER"),))

    build_profile_panel!(rgrid[3, 1], "J_tor [MA/m²]", "J",
                         t_obs, dd_sim, dd_exp, show_sim, show_exp,
                         acc_jtor; scale=1e-6, ylims=ylim_jtor)

    q_ax = build_profile_panel!(rgrid[3, 2], "q", "q",
                                t_obs, dd_sim, dd_exp, show_sim, show_exp,
                                acc_q; from_eq=true)

    # q-axis tracks 1.1 × q95(t), respecting sign convention.
    function q_lims(t)
        qs = Float64[]
        v = q95_at(dd_sim, t); isnan(v) || push!(qs, v)
        if dd_exp !== nothing
            v = q95_at(dd_exp, t); isnan(v) || push!(qs, v)
        end
        isempty(qs) && return (0.0, 6.0)
        if all(>(0), qs)
            return (0.0, 1.1 * maximum(qs))
        elseif all(<(0), qs)
            return (1.1 * minimum(qs), 0.0)
        else
            ext = 1.1 * maximum(abs, qs)
            return (-ext, ext)
        end
    end
    on(t_obs) do t
        try; ylims!(q_ax, q_lims(t)...); catch; end
    end
    ylims!(q_ax, q_lims(t_obs[])...)

    # traces
    tgrid = GridLayout(view_box[2, 1:2])
    build_trace_panel!(tgrid[1, 1], "Ip [MA]", "Ip",
                       t_obs, dd_sim, dd_exp, show_sim, show_exp, trace_ip)
    build_trace_panel!(tgrid[1, 2], "⟨ne⟩ [10¹⁹ m⁻³]", "n",
                       t_obs, dd_sim, dd_exp, show_sim, show_exp, trace_ne_line)
    build_trace_panel!(tgrid[1, 3], "W_MHD [MJ]", "W",
                       t_obs, dd_sim, dd_exp, show_sim, show_exp, trace_wmhd)

    # equilibrium gets a wide slice; trace row stays compact
    colsize!(view_box, 1, Relative(0.36))
    rowsize!(view_box, 1, Relative(0.78))
    return view_box
end

function build_magnetics_view!(view_box::GridLayout,
                               t_obs, dd_sim, dd_exp,
                               show_sim, show_exp)

    src = dd_exp !== nothing ? dd_exp : dd_sim
    probe_pts = magnetics_positions(src, :b_field_pol_probe)
    loop_pts  = magnetics_positions(src, :flux_loop)

    # ----- Sensor map (left) ----------------------------------------------
    # Nested grid: map | probe colorbar | loop colorbar
    left = GridLayout(view_box[1, 1])
    map_ax = Axis(left[1, 1];
                  aspect=DataAspect(),
                  xlabel="R [m]", ylabel="Z [m]",
                  title=lift(t -> "Sensor errors (sim − exp)  —  t = $(round(t; digits=4)) s", t_obs))

    wp = wall_points(src)
    isempty(wp) || lines!(map_ax, wp; color=:gray30, linewidth=1.2)

    if dd_sim !== nothing
        fs = lift(t -> flux_surface_points(dd_sim.equilibrium.time_slice[Float64(t)]), t_obs)
        bnd = lift(t -> boundary_points(dd_sim.equilibrium.time_slice[Float64(t)]), t_obs)
        lines!(map_ax, fs; color=(:gray60, 0.5), linewidth=0.4)
        lines!(map_ax, bnd; color=:gray40, linewidth=1.2)
    end

    have_both = (dd_sim !== nothing) && (dd_exp !== nothing)
    cb_col = 2
    probe_emax = 0.0
    loop_emax  = 0.0

    if have_both
        probe_emax = global_err_range(dd_sim, dd_exp, :b_field_pol_probe)
        loop_emax  = global_err_range(dd_sim, dd_exp, :flux_loop)

        if !isempty(probe_pts)
            probe_err = lift(t_obs) do t
                m = magnetics_values(dd_exp, t, :b_field_pol_probe)
                s = magnetics_values(dd_sim, t, :b_field_pol_probe)
                n = min(length(m), length(s), length(probe_pts))
                n == 0 ? Float64[] : s[1:n] .- m[1:n]
            end
            scatter!(map_ax, probe_pts; color=probe_err, colormap=:RdBu,
                     colorrange=(-probe_emax, probe_emax),
                     marker=:cross, markersize=14,
                     label="Bθ probes ($(length(probe_pts)))")
            Colorbar(left[1, cb_col]; colormap=:RdBu,
                     colorrange=(-probe_emax, probe_emax),
                     label="Bθ residual [T]", width=14, ticklabelsize=10)
            cb_col += 1
        end
        if !isempty(loop_pts)
            loop_err = lift(t_obs) do t
                m = magnetics_values(dd_exp, t, :flux_loop)
                s = magnetics_values(dd_sim, t, :flux_loop)
                n = min(length(m), length(s), length(loop_pts))
                n == 0 ? Float64[] : s[1:n] .- m[1:n]
            end
            scatter!(map_ax, loop_pts; color=loop_err, colormap=:RdBu,
                     colorrange=(-loop_emax, loop_emax),
                     marker=:circle, markersize=10,
                     label="Flux loops ($(length(loop_pts)))")
            Colorbar(left[1, cb_col]; colormap=:RdBu,
                     colorrange=(-loop_emax, loop_emax),
                     label="Ψ residual [Wb]", width=14, ticklabelsize=10)
            cb_col += 1
        end
    else
        isempty(probe_pts) || scatter!(map_ax, probe_pts;
                                       color=:gray, marker=:cross, markersize=14,
                                       label="Bθ probes ($(length(probe_pts)))")
        isempty(loop_pts) || scatter!(map_ax, loop_pts;
                                      color=:gray, marker=:circle, markersize=10,
                                      label="Flux loops ($(length(loop_pts)))")
    end

    isempty(probe_pts) && isempty(loop_pts) ||
        axislegend(map_ax; position=:rt, framevisible=false)

    # ----- Sim/exp by sensor index (right) --------------------------------
    rgrid = GridLayout(view_box[1, 2])

    function build_index_panel!(pos, what::Symbol, title::AbstractString,
                                ylabel::AbstractString)
        ax = Axis(pos; xlabel="sensor #", ylabel=ylabel, title=title)

        n_labeled = 0
        if dd_sim !== nothing
            sim_pts = lift(t_obs) do t
                v = magnetics_values(dd_sim, t, what)
                Point2f.(eachindex(v), v)
            end
            scatterlines!(ax, sim_pts; color=COL_SIM, markersize=6, linewidth=1.0,
                          label="sim", visible=show_sim)
            n_labeled += 1
        end
        if dd_exp !== nothing
            exp_pts = lift(t_obs) do t
                v = magnetics_values(dd_exp, t, what)
                Point2f.(eachindex(v), v)
            end
            scatterlines!(ax, exp_pts; color=COL_EXP, marker=:diamond, markersize=6,
                          linewidth=1.0, label="exp", visible=show_exp)
            n_labeled += 1
        end
        # Fix the y-range to the global span across all times so points don't
        # leave the view while scanning the time slider.
        rng = global_value_range(dd_sim, dd_exp, what)
        if rng !== nothing
            lo, hi = rng
            if hi > lo
                pad = 0.05 * (hi - lo)
                ylims!(ax, lo - pad, hi + pad)
            else
                pad = lo == 0 ? 1.0 : 0.1 * abs(lo)
                ylims!(ax, lo - pad, hi + pad)
            end
        end
        n_labeled > 0 && axislegend(ax; position=:rt, framevisible=false)
        return ax
    end

    build_index_panel!(rgrid[1, 1], :b_field_pol_probe,
                       "Bθ probes  —  sim & exp by sensor", "Bθ [T]")
    build_index_panel!(rgrid[2, 1], :flux_loop,
                       "Flux loops  —  sim & exp by sensor", "Ψ [Wb]")

    colsize!(view_box, 1, Relative(0.45))
    return view_box
end

# ---------------------------------------------------------------------------
# NBI view
# ---------------------------------------------------------------------------

"Name of each NBI unit."
function nbi_unit_names(dd::IMAS.dd)
    names = String[]
    try
        for unit in dd.nbi.unit
            push!(names, try
                isempty(unit.name) ? "unit_$(length(names)+1)" : unit.name
            catch
                "unit_$(length(names)+1)"
            end)
        end
    catch
    end
    return names
end

"""
Launched power [W] vs time for one NBI unit. Prefers the `power_launched`
signal (which carries its own time); falls back to the prescribed waveform
`power.reference`, which shares the IDS-level `nbi.time` passed in as `time`.
Returns empty if data missing or inconsistent.
"""
function nbi_power_trace(unit, time=Float64[])
    # Preferred: power_launched signal (data + its own time)
    try
        t = collect(unit.power_launched.time)
        p = collect(unit.power_launched.data)
        if !isempty(t) && length(t) == length(p)
            return t, p
        end
    catch
    end
    # Fallback: prescribed waveform power.reference on the homogeneous nbi.time
    try
        p = collect(unit.power.reference)
        t = collect(time)
        if !isempty(p) && length(t) == length(p)
            return t, p
        end
    catch
    end
    return Float64[], Float64[]
end

"Total NBI launched power [MW] vs time, summed across units (with interpolation if grids differ)."
function nbi_total_power_trace(dd::IMAS.dd)
    try
        time = collect(dd.nbi.time)
        traces = [nbi_power_trace(u, time) for u in dd.nbi.unit]
        traces = [(t, p) for (t, p) in traces if !isempty(t)]
        isempty(traces) && return Float64[], Float64[]
        ref_t, p1 = traces[1]
        total = copy(p1)
        for (t, p) in traces[2:end]
            if t == ref_t
                total .+= p
            else
                total .+= IMAS.interp1d(t, p).(ref_t)
            end
        end
        return ref_t, total ./ 1e6
    catch
        return Float64[], Float64[]
    end
end

"x, y of a circle of radius `R` centered at the origin (top-view helper)."
function _circle_xy(R; n::Int=200)
    θ = range(0, 2π; length=n)
    return R .* cos.(θ), R .* sin.(θ)
end

"""
Per-unit-name NBI beamline geometry from `dd.nbi.unit[:].beamlets_group[1]`:
`(R, ϕ, Rtan, dir)` — source major radius & toroidal angle, tangency radius,
and injection direction (±1). Empty if no geometry (e.g. dd_sim stores only
power/energy; dd_exp carries the machine description).
"""
function nbi_beam_geometry(dd)
    geom = Dict{String,NamedTuple}()
    dd === nothing && return geom
    try
        names = nbi_unit_names(dd)
        for (i, u) in enumerate(dd.nbi.unit)
            try
                isempty(u.beamlets_group) && continue
                g = u.beamlets_group[1]
                geom[names[i]] = (R=Float64(g.position.r), ϕ=Float64(g.position.phi),
                                  Rtan=Float64(g.tangency_radius), dir=Float64(g.direction))
            catch
            end
        end
    catch
    end
    return geom
end

"""
Draw the NBI beamlines in a top (toroidal) view: wall + plasma boundary as
concentric circles, each beam as a straight chord tangent to its tangency
radius. `on_names` get a solid colored line, the rest a faint gray one.
`geom` is the per-name geometry dict; wall/plasma come from `geo_dd`.
"""
function build_nbi_topview!(pos, t_obs, geom, geo_dd, names_sim, color_for, on_names)
    ax = Axis(pos; title="top view", titlesize=11, titlegap=2,
              xlabel="X [m]", ylabel="Y [m]", aspect=DataAspect())
    tnow = Makie.to_value(t_obs)

    # wall as inner/outer circles (R at the midplane)
    r_out = 2.5
    try
        fw = IMAS.first_wall(geo_dd.wall)
        if !isempty(fw.r)
            r_in, r_out = extrema(fw.r)
            lines!(ax, _circle_xy(r_in)...;  color=:black, linewidth=1.2)
            lines!(ax, _circle_xy(r_out)...; color=:black, linewidth=1.2)
        end
    catch
    end

    # plasma boundary (snapshot at t_obs) as inner/outer circles + magnetic axis
    try
        eqt = geo_dd.equilibrium.time_slice[tnow]
        rb = extrema(eqt.boundary.outline.r)
        lines!(ax, _circle_xy(rb[1])...; color=(:dodgerblue, 0.9), linewidth=1.0)
        lines!(ax, _circle_xy(rb[2])...; color=(:dodgerblue, 0.9), linewidth=1.0)
        RA = eqt.global_quantities.magnetic_axis.r
        lines!(ax, _circle_xy(RA)...; color=(:gray, 0.6), linewidth=0.8, linestyle=:dash)
    catch
    end

    # beamlines: chord through the source tangent to the tangency-radius circle
    for (i, nm) in enumerate(names_sim)
        haskey(geom, nm) || continue
        g = geom[nm]
        (g.R <= 0 || g.Rtan > g.R) && continue
        α  = asin(clamp(g.Rtan / g.R, -1.0, 1.0))
        ur = (-cos(g.ϕ), -sin(g.ϕ))          # inward radial
        ut = (-sin(g.ϕ),  cos(g.ϕ))          # toroidal (+ϕ̂)
        d  = (cos(α) * ur[1] + g.dir * sin(α) * ut[1],
              cos(α) * ur[2] + g.dir * sin(α) * ut[2])
        S  = (g.R * cos(g.ϕ), g.R * sin(g.ϕ))
        L  = 2 * sqrt(max(g.R^2 - g.Rtan^2, 0.0))
        P  = (S[1] + L * d[1], S[2] + L * d[2])
        if nm in on_names
            lines!(ax, [S[1], P[1]], [S[2], P[2]]; color=color_for(i), linewidth=2.0)
        else
            lines!(ax, [S[1], P[1]], [S[2], P[2]]; color=(:gray, 0.4), linewidth=0.8)
        end
    end

    lim = r_out * 1.12
    limits!(ax, -lim, lim, -lim, lim)
    return ax
end

"""
NBI view — one panel per unit showing launched power vs time (sim solid,
exp dashed), a total panel, and a top-view schematic of the beamlines.
"""
function build_nbi_view!(view_box::GridLayout, t_obs, dd_sim, dd_exp, show_sim, show_exp)
    if dd_sim === nothing || isempty(dd_sim.nbi.unit)
        Label(view_box[1, 1], "No NBI data in dd_sim"; halign=:center, fontsize=14)
        return view_box
    end

    palette = Makie.wong_colors()
    color_for(i) = palette[mod1(i, length(palette))]

    names_sim = nbi_unit_names(dd_sim)
    names_exp = dd_exp !== nothing ? nbi_unit_names(dd_exp) : String[]

    # Build a lookup: exp unit name → unit (for matching to sim by name)
    exp_by_name = Dict{String,Any}()
    if dd_exp !== nothing
        for (i, unit) in enumerate(dd_exp.nbi.unit)
            exp_by_name[names_exp[i]] = unit
        end
    end

    # IDS-level time bases (power.reference shares nbi.time, has no own time)
    t_sim = try collect(dd_sim.nbi.time) catch; Float64[] end
    t_exp = dd_exp !== nothing ? (try collect(dd_exp.nbi.time) catch; Float64[] end) : Float64[]

    # Beam panels + total as a 3-column small-multiples grid on the left; a
    # top-view beamline schematic occupies a narrow right-hand column.
    # With ~7 rows, a full set of per-panel axis decorations (title + "t [s]" +
    # tick labels on every cell) eats the vertical space and collapses the plot
    # areas, so we share them: x labels only on the lowest panel of each column,
    # y labels only on the left column, and one common y-scale across the beams.
    n_units  = length(dd_sim.nbi.unit)
    n_panels = n_units + 1           # beams + total, all in the main grid
    ncols    = min(n_panels, 3)
    grid     = GridLayout(view_box[1, 1])

    pos_of(idx) = (div(idx - 1, ncols) + 1, mod(idx - 1, ncols) + 1)  # idx → (row, col)

    # lowest occupied row in each column → the only panel there that keeps its
    # x-axis labels (the rest of the column shares them)
    bottom_row = Dict{Int,Int}()
    for idx in 1:n_panels
        r, c = pos_of(idx)
        bottom_row[c] = max(get(bottom_row, c, 0), r)
    end

    # `nbi.time` only covers the analysis window (~0.3 s) while the per-unit
    # power traces span the full shot (~4 s), so derive the shared x-range from
    # the actual traces, not from nbi.time.
    beam_axes = Axis[]
    beam_max  = Dict{String,Float64}()  # peak MW per unit → on/off in the top view
    gymax = Ref(0.0)              # shared y-max across all beam panels
    gtmin = Ref(Inf)             # shared x-range, tracked from real traces
    gtmax = Ref(-Inf)
    function track_t!(t)
        isempty(t) && return
        gtmin[] = min(gtmin[], minimum(t))
        gtmax[] = max(gtmax[], maximum(t))
    end

    function add_beam!(idx, title, color, unit_sim, unit_exp)
        r, c = pos_of(idx)
        ax = Axis(grid[r, c]; title=title, titlesize=11, titlegap=2,
                  xlabel="t [s]", ylabel="Power [MW]")
        n = 0
        pk = 0.0
        t, p = nbi_power_trace(unit_sim, t_sim)
        if !isempty(t)
            pmw = p ./ 1e6
            lines!(ax, t, pmw; color=color, label="sim", visible=show_sim)
            pk = max(pk, maximum(pmw))
            track_t!(t)
            n += 1
        end
        if unit_exp !== nothing
            t_e, p_e = nbi_power_trace(unit_exp, t_exp)
            if !isempty(t_e)
                pemw = p_e ./ 1e6
                lines!(ax, t_e, pemw; color=color, linestyle=:dash,
                       label="exp", visible=show_exp)
                pk = max(pk, maximum(pemw))
                track_t!(t_e)
                n += 1
            end
        end
        beam_max[title] = pk
        gymax[] = max(gymax[], pk)
        vlines!(ax, t_obs; color=:gray, linewidth=1.0)
        # share decorations
        if r != bottom_row[c]
            ax.xticklabelsvisible = false
            ax.xlabelvisible = false
        end
        if c != 1
            ax.yticklabelsvisible = false
            ax.ylabelvisible = false
        end
        # one legend (on the first panel) is enough to show solid=sim / dashed=exp
        idx == 1 && n > 0 && axislegend(ax; position=:lt, framevisible=false,
                                        labelsize=8, padding=(3, 3, 2, 2))
        push!(beam_axes, ax)
        return ax
    end

    for (i, unit) in enumerate(dd_sim.nbi.unit)
        exp_unit = get(exp_by_name, names_sim[i], nothing)
        add_beam!(i, names_sim[i], color_for(i), unit, exp_unit)
    end

    # Total panel — back in the main grid (bottom-right). Different magnitude,
    # so it keeps its own y-scale and labels.
    trow, tcol = pos_of(n_panels)
    ax_tot = Axis(grid[trow, tcol]; title="NBI — total", titlesize=11, titlegap=2,
                  xlabel="t [s]", ylabel="Power [MW]")
    ymax_tot = 0.0
    t_tot, p_tot = nbi_total_power_trace(dd_sim)
    if !isempty(t_tot)
        lines!(ax_tot, t_tot, p_tot; color=:black, linewidth=2,
               label="sim", visible=show_sim)
        ymax_tot = maximum(p_tot)
        track_t!(t_tot)
    end
    vlines!(ax_tot, t_obs; color=:gray, linewidth=1.0)
    ylims!(ax_tot, 0, ymax_tot > 0 ? ymax_tot * 1.15 : 2.0)  # was hardcoded 0–2 MW, which clipped the total

    # ---- right-hand column: top-view beamline schematic (full height) --------
    on_thresh = max(0.05, 0.02 * gymax[])
    on_names  = Set(nm for (nm, m) in beam_max if m > on_thresh)
    geom = nbi_beam_geometry(dd_exp)             # dd_exp carries the machine geometry
    isempty(geom) && (geom = nbi_beam_geometry(dd_sim))
    build_nbi_topview!(view_box[1, 2], t_obs, geom, dd_sim, names_sim, color_for, on_names)

    # shared scales: one y-range for the beam panels, one x-range (the real
    # data span) for every panel including the total
    gylim = gymax[] > 0 ? gymax[] * 1.15 : 1.0
    have_x = gtmin[] < gtmax[]
    for ax in beam_axes
        ylims!(ax, 0, gylim)
        have_x && xlims!(ax, gtmin[], gtmax[])
    end
    have_x && xlims!(ax_tot, gtmin[], gtmax[])

    rowgap!(grid, 4)
    colgap!(grid, 10)
    colsize!(view_box, 2, Relative(0.22))

    return view_box
end

# ---------------------------------------------------------------------------
# EC view
# ---------------------------------------------------------------------------

"Name of each EC launcher beam."
function ec_beam_names(dd)
    names = String[]
    try
        for b in dd.ec_launchers.beam
            push!(names, try
                isempty(b.name) ? "beam_$(length(names)+1)" : b.name
            catch
                "beam_$(length(names)+1)"
            end)
        end
    catch
    end
    return names
end

"Launched power [W] vs time for one EC beam (empty if missing or inconsistent)."
function ec_power_trace(beam)
    try
        t = collect(beam.power_launched.time)
        p = collect(beam.power_launched.data)
        (isempty(t) || length(t) != length(p)) && return Float64[], Float64[]
        return t, p
    catch
        return Float64[], Float64[]
    end
end

"Total EC launched power [MW] vs time, summed across beams (interpolating if grids differ)."
function ec_total_power_trace(dd)
    try
        traces = [ec_power_trace(b) for b in dd.ec_launchers.beam]
        traces = [(t, p) for (t, p) in traces if !isempty(t)]
        isempty(traces) && return Float64[], Float64[]
        ref_t, p1 = traces[1]
        total = copy(p1)
        for (t, p) in traces[2:end]
            total .+= (t == ref_t ? p : IMAS.interp1d(t, p).(ref_t))
        end
        return ref_t, total ./ 1e6
    catch
        return Float64[], Float64[]
    end
end

"""
R,Z ray path of one EC beam from `waves.coherent_wave[:].beam_tracing` (matched
by antenna name). Empty if no beam tracing is available at time `t0`.
"""
function ec_beam_path(dd, name, t0)
    try
        for cw in dd.waves.coherent_wave
            nm = try cw.identifier.antenna_name catch; "" end
            if nm == name && !isempty(cw.beam_tracing) && t0 >= cw.beam_tracing[1].time
                b1 = cw.beam_tracing[t0].beam[1]
                return collect(b1.position.r), collect(b1.position.z)
            end
        end
    catch
    end
    return Float64[], Float64[]
end

"""
Poloidal (R-Z) cross-section: wall + plasma boundary outline plus each EC beam's
ray path (from coherent_wave beam tracing), launch points marked. `on_names` get
a solid colored line, the rest a faint gray one.
"""
function build_ec_crosssection!(pos, t_obs, dd, names, color_for, on_names)
    ax = Axis(pos; title="cross-section", titlesize=11, titlegap=2,
              xlabel="R [m]", ylabel="Z [m]", aspect=DataAspect())
    tnow = Makie.to_value(t_obs)
    rlim = (1.0, 2.5); zlim = (-1.5, 1.5)
    try
        fw = IMAS.first_wall(dd.wall)
        if !isempty(fw.r)
            lines!(ax, fw.r, fw.z; color=:black, linewidth=1.2)
            rlim = extrema(fw.r); zlim = extrema(fw.z)
        end
    catch
    end
    try
        bnd = dd.equilibrium.time_slice[tnow].boundary.outline
        lines!(ax, collect(bnd.r), collect(bnd.z); color=(:dodgerblue, 0.9), linewidth=1.0)
    catch
    end
    for (i, nm) in enumerate(names)
        r, z = ec_beam_path(dd, nm, tnow)
        isempty(r) && continue
        col = nm in on_names ? color_for(i) : (:gray, 0.4)
        lines!(ax, r, z; color=col, linewidth=2.0)
        scatter!(ax, [r[1]], [z[1]]; color=col, markersize=6)
    end
    limits!(ax, rlim[1] - 0.15, rlim[2] + 0.25, zlim[1] * 1.1, zlim[2] * 1.1)
    return ax
end

"""
EC view — one panel per launcher showing launched power vs time (sim solid,
exp dashed), a total panel, and a poloidal cross-section of the beam paths.
Mirrors the NBI view's shared-axis small-multiples layout.
"""
function build_ec_view!(view_box::GridLayout, t_obs, dd_sim, dd_exp, show_sim, show_exp)
    if dd_sim === nothing || isempty(dd_sim.ec_launchers.beam)
        Label(view_box[1, 1], "No EC data in dd_sim"; halign=:center, fontsize=14)
        return view_box
    end

    palette = Makie.wong_colors()
    color_for(i) = palette[mod1(i, length(palette))]

    names_sim = ec_beam_names(dd_sim)
    names_exp = dd_exp !== nothing ? ec_beam_names(dd_exp) : String[]
    exp_by_name = Dict{String,Any}()
    if dd_exp !== nothing
        for (i, b) in enumerate(dd_exp.ec_launchers.beam)
            exp_by_name[names_exp[i]] = b
        end
    end

    n_beams  = length(dd_sim.ec_launchers.beam)
    n_panels = n_beams + 1           # beams + total, all in the main grid
    ncols    = min(n_panels, 3)
    grid     = GridLayout(view_box[1, 1])

    pos_of(idx) = (div(idx - 1, ncols) + 1, mod(idx - 1, ncols) + 1)
    bottom_row = Dict{Int,Int}()
    for idx in 1:n_panels
        r, c = pos_of(idx)
        bottom_row[c] = max(get(bottom_row, c, 0), r)
    end

    beam_axes = Axis[]
    beam_max  = Dict{String,Float64}()
    gymax = Ref(0.0); gtmin = Ref(Inf); gtmax = Ref(-Inf)
    function track_t!(t)
        isempty(t) && return
        gtmin[] = min(gtmin[], minimum(t)); gtmax[] = max(gtmax[], maximum(t))
    end

    function add_beam!(idx, title, color, beam_sim, beam_exp)
        r, c = pos_of(idx)
        ax = Axis(grid[r, c]; title=title, titlesize=11, titlegap=2,
                  xlabel="t [s]", ylabel="Power [MW]")
        n = 0; pk = 0.0
        t, p = ec_power_trace(beam_sim)
        if !isempty(t)
            pmw = p ./ 1e6
            lines!(ax, t, pmw; color=color, label="sim", visible=show_sim)
            pk = max(pk, maximum(pmw)); track_t!(t); n += 1
        end
        if beam_exp !== nothing
            te, pe = ec_power_trace(beam_exp)
            if !isempty(te)
                pemw = pe ./ 1e6
                lines!(ax, te, pemw; color=color, linestyle=:dash, label="exp", visible=show_exp)
                pk = max(pk, maximum(pemw)); track_t!(te); n += 1
            end
        end
        beam_max[title] = pk; gymax[] = max(gymax[], pk)
        vlines!(ax, t_obs; color=:gray, linewidth=1.0)
        if r != bottom_row[c]
            ax.xticklabelsvisible = false; ax.xlabelvisible = false
        end
        if c != 1
            ax.yticklabelsvisible = false; ax.ylabelvisible = false
        end
        idx == 1 && n > 0 && axislegend(ax; position=:lt, framevisible=false,
                                        labelsize=8, padding=(3, 3, 2, 2))
        push!(beam_axes, ax)
        return ax
    end

    for (i, b) in enumerate(dd_sim.ec_launchers.beam)
        add_beam!(i, names_sim[i], color_for(i), b, get(exp_by_name, names_sim[i], nothing))
    end

    # Total panel (own y-scale) in the main grid, bottom slot
    trow, tcol = pos_of(n_panels)
    ax_tot = Axis(grid[trow, tcol]; title="EC — total", titlesize=11, titlegap=2,
                  xlabel="t [s]", ylabel="Power [MW]")
    ymax_tot = 0.0
    t_tot, p_tot = ec_total_power_trace(dd_sim)
    if !isempty(t_tot)
        lines!(ax_tot, t_tot, p_tot; color=:black, linewidth=2, label="sim", visible=show_sim)
        ymax_tot = maximum(p_tot); track_t!(t_tot)
    end
    vlines!(ax_tot, t_obs; color=:gray, linewidth=1.0)
    ylims!(ax_tot, 0, ymax_tot > 0 ? ymax_tot * 1.15 : 2.0)

    # right column: poloidal cross-section of the beam paths
    on_thresh = max(0.05, 0.02 * gymax[])
    on_names  = Set(nm for (nm, m) in beam_max if m > on_thresh)
    build_ec_crosssection!(view_box[1, 2], t_obs, dd_sim, names_sim, color_for, on_names)

    gylim = gymax[] > 0 ? gymax[] * 1.15 : 1.0
    have_x = gtmin[] < gtmax[]
    for ax in beam_axes
        ylims!(ax, 0, gylim); have_x && xlims!(ax, gtmin[], gtmax[])
    end
    have_x && xlims!(ax_tot, gtmin[], gtmax[])

    rowgap!(grid, 4)
    colgap!(grid, 10)
    colsize!(view_box, 2, Relative(0.24))

    return view_box
end

"""
Sources view — 5-panel breakdown of `dd.core_sources` channels (Pe, Pi, J‖,
Γe, Mφ). When `aggregate_hcd=true` (default), individual HCD sources (NBI,
EC, LH, IC, pellet) are summed per category so there is at most one line per
HCD type instead of one per beam/antenna. Non-HCD sources (ohmic, bootstrap,
radiation, etc.) are always shown individually. A bold black total line comes
from `IMAS.total_sources`; if `dd_exp` is present its total is dashed black.

Per-panel thresholds drop sources/categories below the limit:
- `min_power` — Pe and Pi (W)
- `min_particles` — Γe (s⁻¹)
- `min_momentum` — Mφ (N·m)
J‖ is unfiltered.
"""
function build_sources_view!(view_box::GridLayout,
                             t_obs, dd_sim, dd_exp,
                             show_sim, show_exp;
                             min_power::Real=0.0,
                             min_particles::Real=0.0,
                             min_momentum::Real=0.0,
                             aggregate_hcd::Bool=true)

    if dd_sim === nothing || isempty(dd_sim.core_sources.source)
        Label(view_box[1, 1], "No core_sources data in dd_sim";
              halign=:center, fontsize=14)
        return view_box
    end

    grid = GridLayout(view_box[1, 1])

    times = collect(dd_sim.core_sources.time)
    isempty(times) && (times = collect(dd_sim.equilibrium.time))

    names   = source_names(dd_sim)
    indices = source_indices(dd_sim)
    palette = Makie.wong_colors()
    color_for(i) = palette[mod1(i, length(palette))]

    function build_channel!(pos, accessor, title, ylabel;
                            scale::Real=1.0,
                            integration::Symbol=:volume,
                            threshold::Real=0.0,
                            unit_label::String="",
                            unit_scale::Real=1.0)
        ax = Axis(pos; xlabel="ρ_tor_norm", ylabel=ylabel, title=title)
        n_labeled = 0

        # Individual sources (non-HCD, or all when aggregate_hcd=false)
        for (i, name) in enumerate(names)
            aggregate_hcd && get(indices, i, 0) ∈ _HCD_INDEXES && continue

            integ = source_max_integral(dd_sim.core_sources.source[i],
                                        accessor, times; integration)
            integ > threshold || continue

            pts = lift(t_obs) do t
                source_profile_points(dd_sim.core_sources.source[i], t, accessor; scale)
            end
            label = name
            lines!(ax, pts; color=color_for(i), linewidth=1.0,
                   label=label, visible=show_sim)
            n_labeled += 1
        end

        # One aggregated line per HCD category
        if aggregate_hcd
            for cat in _HCD_CATEGORIES
                js = findall(j -> get(indices, j, 0) == cat.idx, 1:length(indices))
                isempty(js) && continue

                cat_max = sum(source_max_integral(dd_sim.core_sources.source[j],
                                                  accessor, times; integration)
                              for j in js; init=0.0)
                cat_max > threshold || continue

                pts = lift(t_obs) do t
                    hcd_category_points(dd_sim, t, accessor, [cat.idx]; scale)
                end
                label = cat.name
                lines!(ax, pts; color=cat.color, linewidth=1.5,
                       label=label, visible=show_sim)
                n_labeled += 1
            end
        end

        tot_sim = lift(t -> total_source_points(dd_sim, t, accessor; scale), t_obs)
        lines!(ax, tot_sim; color=:black, linewidth=2.0,
               label="total (sim)", visible=show_sim)
        n_labeled += 1

        if dd_exp !== nothing
            tot_exp = lift(t -> total_source_points(dd_exp, t, accessor; scale), t_obs)
            lines!(ax, tot_exp; color=:black, linewidth=2.0,
                   linestyle=:dash, label="total (exp)", visible=show_exp)
            n_labeled += 1
        end

        if n_labeled > 0
            axislegend(ax; position=:rt, framevisible=false,
                       labelsize=9, rowgap=1, padding=(4, 4, 2, 2),
                       merge=true)
        end
        return ax
    end

    build_channel!(grid[1, 1], acc_src_pe,   "Pe — electron heating",   "[MW/m³]";
                   scale=1e-6, integration=:volume,
                   threshold=min_power, unit_label="MW", unit_scale=1e-6)
    build_channel!(grid[1, 2], acc_src_pi,   "Pi — ion heating",        "[MW/m³]";
                   scale=1e-6, integration=:volume,
                   threshold=min_power, unit_label="MW", unit_scale=1e-6)
    build_channel!(grid[1, 3], acc_src_jpar, "J∥ — parallel current",   "[MA/m²]";
                   scale=1e-6, integration=:area,
                   unit_label="MA", unit_scale=1e-6)
    build_channel!(grid[2, 1], acc_src_part, "Γe — electron particles", "[10¹⁹/m³/s]";
                   scale=1e-19, integration=:volume,
                   threshold=min_particles, unit_label="10¹⁹/s", unit_scale=1e-19)
    build_channel!(grid[2, 2], acc_src_mom,  "Mφ — toroidal momentum",  "[N/m²]";
                   integration=:volume,
                   threshold=min_momentum, unit_label="N·m")

    return view_box
end

# ---------------------------------------------------------------------------
# core_transport view
# ---------------------------------------------------------------------------

"""
core_transport view — 4-panel breakdown of transport fluxes (Qe, Qi, Γe, Πφ).

Each panel shows:
  • one faint line per transport model in `dd.core_transport.model`
    (skipping aggregate/placeholder types: combined, transport_solver, etc.)
  • a bold red "total transport" line summing all models via `IMAS.total_fluxes`
  • a bold black "source-implied" line — the steady-state flux you get from
    `power_inside / surface` (resp. particles_inside, torque_tor_inside).
    If `dd_exp` is provided, its source-implied flux is dashed.

If no transport models are populated (e.g., the turbulent transport actor
wasn't run), only the source-implied lines are shown.
"""
function build_transport_view!(view_box::GridLayout,
                                t_obs, dd_sim, dd_exp,
                                show_sim, show_exp)
    if dd_sim === nothing
        Label(view_box[1, 1], "No data loaded"; halign=:center, fontsize=14)
        return view_box
    end

    has_models  = !isempty(dd_sim.core_transport.model)
    has_sources = !isempty(dd_sim.core_sources.source)
    if !has_models && !has_sources
        Label(view_box[1, 1], "No core_transport or core_sources data in dd_sim";
              halign=:center, fontsize=14)
        return view_box
    end

    grid = GridLayout(view_box[1, 1])
    palette = Makie.wong_colors()
    color_for(i) = palette[mod1(i, length(palette))]

    names   = transport_model_names(dd_sim)
    indices = transport_model_indices(dd_sim)

    function build_channel!(pos, accessor_flux, accessor_inside, title, ylabel;
                            scale::Real=1.0)
        ax = Axis(pos; xlabel="ρ_tor_norm", ylabel=ylabel, title=title)
        n = 0

        # Per-model lines (only if any models are loaded)
        for (i, name) in enumerate(names)
            indices[i] ∈ _CT_SKIP_INDEXES && continue
            pts = lift(t_obs) do t
                transport_flux_points(dd_sim.core_transport.model[i], t,
                                      accessor_flux; scale)
            end
            lines!(ax, pts; color=color_for(i), linewidth=1.0,
                   label=name, visible=show_sim)
            n += 1
        end

        # Total transport (sum of models). Only meaningful if ≥ 1 real model.
        n_real_models = count(i -> !(i ∈ _CT_SKIP_INDEXES), indices)
        if n_real_models > 0
            tot_t = lift(t_obs) do t
                total_transport_flux_points(dd_sim, t, accessor_flux; scale)
            end
            lines!(ax, tot_t; color=:red, linewidth=2.0,
                   label="total transport", visible=show_sim)
            n += 1
        end

        # Source-implied flux (sim) — always shown when sources exist
        if has_sources
            tot_s = lift(t_obs) do t
                source_implied_flux_points(dd_sim, t, accessor_inside; scale)
            end
            lines!(ax, tot_s; color=:black, linewidth=2.0,
                   label="source (sim)", visible=show_sim)
            n += 1
        end

        # Source-implied flux (exp), dashed
        if dd_exp !== nothing && !isempty(dd_exp.core_sources.source)
            tot_e = lift(t_obs) do t
                source_implied_flux_points(dd_exp, t, accessor_inside; scale)
            end
            lines!(ax, tot_e; color=:black, linewidth=2.0, linestyle=:dash,
                   label="source (exp)", visible=show_exp)
            n += 1
        end

        n > 0 && axislegend(ax; position=:rt, framevisible=false,
                            labelsize=9, rowgap=1, padding=(4,4,2,2),
                            merge=true)
        return ax
    end

    build_channel!(grid[1, 1], acc_qe_flux,   acc_qe_inside,
                   "Qe — electron heat flux",      "[MW/m²]"; scale=1e-6)
    build_channel!(grid[1, 2], acc_qi_flux,   acc_qi_inside,
                   "Qi — ion heat flux",           "[MW/m²]"; scale=1e-6)
    build_channel!(grid[2, 1], acc_part_flux, acc_part_inside,
                   "Γe — electron particle flux",  "[10¹⁹/m²/s]"; scale=1e-19)
    build_channel!(grid[2, 2], acc_mom_flux,  acc_mom_inside,
                   "Πφ — toroidal momentum flux",  "[N/m]")

    return view_box
end

# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

const VIEWS = ["Summary", "Magnetics", "Sources", "NBI", "EC", "Transport"]

"""
    launch(dd_sim::IMAS.dd, dd_exp::Union{Nothing,IMAS.dd}=nothing; kwargs...)

Open the FuseViewer window. `dd_sim` is the simulation result (drawn in
steel blue); pass `dd_exp` as the experimental/reference dataset (drawn
in black) for side-by-side comparison.

Pick views from the `view:` dropdown in the control row. Currently
implemented: Summary, Magnetics.

`display_window=false` returns the Figure without calling `display` — use it
from a notebook that renders the figure itself. The widgets only respond
under an interactive backend, so the default (a GLMakie window) is still the
one to use for scrubbing through time.
"""
function launch(dd_sim::IMAS.dd, dd_exp::Union{Nothing,IMAS.dd}=nothing;
                time::Union{Nothing,Real}=nothing,
                view::AbstractString="Summary",
                display_window::Bool=true,
                size=(1750, 1100),
                aggregate_hcd::Bool=true,
                min_power::Real=0.0,
                min_particles::Real=0.0,
                min_momentum::Real=0.0)

    times = collect(dd_sim.equilibrium.time)
    isempty(times) && error("dd_sim.equilibrium.time is empty — nothing to show.")
    t0 = something(time, times[length(times) ÷ 2 + 1])

    # Keep resizing usable over a remote display (NoMachine/VNC), where every
    # intermediate size in a drag re-renders and ships the whole window over
    # the wire. scalefactor pins us to one framebuffer pixel per screen pixel,
    # which a HiDPI client would otherwise push to 2x — 4x the pixels to draw
    # and to send; framerate caps how many of those frames go out mid-drag.
    #
    # Do NOT add oit/ssao/fxaa here. They look like they should help, but
    # GLFramebuffer allocates color/objectid/HDR_color/OIT_weight/depth
    # unconditionally (glwindow.jl), so switching the passes off reallocates
    # exactly as much per resize step. Worse, apply_config! reacts to a
    # changed flag by destroying the corresponding postprocessor, and that
    # frees framebuffer textures it only borrowed (:color for fxaa) while
    # leaving them registered in fb.buffers — the next resize! then binds a
    # freed texture and throws. Reusing a screen across launch() calls is
    # what makes the flag "change".
    #
    # Still sluggish? px_per_unit=0.6 renders into a smaller buffer and
    # upscales: ~a third of the fragment work, at the cost of softer text.
    GLMakie.activate!(; scalefactor=1.0f0, framerate=15.0, vsync=false,
                        render_on_demand=true)

    fig = Figure(; size, fontsize=13)

    # ----- top control row ------------------------------------------------
    ctrl = GridLayout(fig[1, 1:2]; tellwidth=false)
    Label(ctrl[1, 1], "FUSE Viewer";
          font=:bold, fontsize=18, halign=:left)

    show_sim = Observable(true)
    show_exp = Observable(dd_exp !== nothing)

    sim_toggle = Toggle(ctrl[1, 2]; active=show_sim[])
    Label(ctrl[1, 3], "sim"; halign=:left)
    on(sim_toggle.active) do v; show_sim[] = v; end

    exp_toggle = Toggle(ctrl[1, 4];
                        active=show_exp[],
                        buttoncolor=(dd_exp === nothing ? :gray70 : :black))
    Label(ctrl[1, 5], "exp + diag"; halign=:left)
    on(exp_toggle.active) do v; show_exp[] = v; end

    Label(ctrl[1, 6], "view:"; halign=:right)
    view_menu = Menu(ctrl[1, 7]; options=VIEWS, default=view, width=150)

    play_btn = Button(ctrl[1, 8]; label="▶ play")

    Label(ctrl[1, 9], "load dir:"; halign=:right)
    load_tb = Textbox(ctrl[1, 10];
                      placeholder="path to shot directory — Enter to load",
                      validator=_is_loadable_dir,
                      width=330)
    on(load_tb.stored_string) do new_dir
        isempty(strip(new_dir)) && return
        @info "FuseViewer: loading from $new_dir"
        try
            open_dir(new_dir)
        catch err
            @error "FuseViewer: load failed" exception=(err, catch_backtrace())
        end
    end
    colgap!(ctrl, 8)

    # ----- time slider ---------------------------------------------------
    sg = SliderGrid(fig[3, 1:2],
                    (label="time [s]", range=times, startvalue=t0,
                     format=x -> string(round(x; digits=4))))
    t_slider = sg.sliders[1]
    t_obs = t_slider.value

    # ----- swappable view box --------------------------------------------
    view_box = GridLayout()
    fig[2, 1:2] = view_box

    # Snapshot all Blocks that exist in the figure *before* any view is built —
    # these are the controls and slider that must persist across view switches.
    # Blocks created afterwards (axes, axislegend's Legend, in-axis Label, etc.)
    # all live in `fig.content`; on switch we delete every Block whose objectid
    # isn't in this set. axislegend's `Legend` is parented to the figure rather
    # than the axis, so simply deleting axes leaves stale legends behind — this
    # walk catches them.
    persistent_blocks = Set(map(objectid, fig.content))

    function rebuild_view(name::AbstractString)
        for b in copy(fig.content)
            if !(objectid(b) in persistent_blocks)
                try; delete!(b); catch err
                    @warn "FuseViewer: could not delete block" exception=err
                end
            end
        end
        _clear_layout!(view_box)  # also drop any leftover sub-GridLayouts

        if name == "Summary"
            build_summary_view!(view_box, t_obs, dd_sim, dd_exp,
                                show_sim, show_exp, times)
        elseif name == "Magnetics"
            build_magnetics_view!(view_box, t_obs, dd_sim, dd_exp,
                                  show_sim, show_exp)
        elseif name == "Sources"
            build_sources_view!(view_box, t_obs, dd_sim, dd_exp,
                                show_sim, show_exp;
                                aggregate_hcd, min_power, min_particles, min_momentum)
        elseif name == "NBI"
            build_nbi_view!(view_box, t_obs, dd_sim, dd_exp, show_sim, show_exp)
        elseif name == "EC"
            build_ec_view!(view_box, t_obs, dd_sim, dd_exp, show_sim, show_exp)
        elseif name == "Transport"
            build_transport_view!(view_box, t_obs, dd_sim, dd_exp, show_sim, show_exp)
        else
            Label(view_box[1, 1], "Unknown view: $name"; halign=:center)
        end
    end

    on(view_menu.selection) do v; rebuild_view(v); end
    rebuild_view(view)

    # ----- play button ---------------------------------------------------
    playing = Ref(false)
    on(play_btn.clicks) do _
        playing[] = !playing[]
        play_btn.label[] = playing[] ? "⏸ pause" : "▶ play"
        if playing[]
            @async begin
                while playing[]
                    idx = findfirst(==(t_obs[]), times)
                    idx = idx === nothing ? 1 : (idx == length(times) ? 1 : idx + 1)
                    set_close_to!(t_slider, times[idx])
                    sleep(0.05)
                end
            end
        end
    end

    if display_window
        _active_screen[] = display(fig)
    end
    _active_fig[] = fig
    return fig
end

"""
    wait_for_close()

Block the calling thread until all FuseViewer windows have been closed.
Use this in non-interactive scripts (without `julia -i`) so Julia does not
exit while the viewer is still open.

    julia --project=@v1.11 fuseviewer/setup.jl   # no -i needed
"""
function wait_for_close(; timeout=30.0)
    # Wait up to `timeout` seconds for the first screen to appear.
    t0 = time()
    while isnothing(_active_screen[])
        time() - t0 > timeout && return
        sleep(0.1)
    end
    # Follow window switches: open_dir calls closeall() then sets a new screen
    # via launch(). wait(screen) blocks until GLMakie's render task exits.
    prev = nothing
    while true
        s = _active_screen[]
        (isnothing(s) || s === prev) && break
        prev = s
        try
            wait(s)      # blocks until this window's render loop exits
        catch
            # GLMakie throws "Context is not alive anymore!" when open_dir
            # calls closeall() to replace the splash with the viewer — ignore.
        end
        sleep(0.5)       # give open_dir time to set _active_screen[] to the new viewer
    end
end

end # module
