"""
Hubbard model on the honeycomb lattice (graphene): AFM transition at half-filling.

Model:
  H = -t Σ_{<ij>,σ} (c†_{iσ}c_{jσ} + h.c.) + U Σ_i n_{i↑}n_{i↓}

Honeycomb unit cell: 2 sublattices (A,B). Mean-field predicts an AFM transition
at U/t ≈ 2.2(3), where opposite sublattice moments form and a gap opens.

This example:
  1) U=0 graphene bands along Γ–K–M–Γ.
  2) U=0 zigzag cylinder bands (edge states).
  3) Interacting Hubbard model: AFM transition and gap opening.

Run:
    julia --project=examples -t 8 examples/SM_AFM/run.jl
"""

using Printf
using LinearAlgebra
using MeanFieldTheories
using Plots

# ── Parameters ───────────────────────────────────────────────────────────────
const t  = 1.0
const Uc = 2.2

U_sweep = collect(range(0.0, 4.0, length=41))
U_bands = [0.0, 1.0, 2.2, 2.3, 3.0, 4.0]

# ── Lattice geometry ───────────────────────────────────────────────────────────
const a  = 1.0
const a1 = [a, 0.0]
const a2 = [0.5*a, √3/2*a]

# Honeycomb unit cell: A at origin, B shifted by δ = (0, a/√3)
# Include a size-1 cell DOF in the unit cell
unitcell = Lattice(
    [Dof(:cell, 1), Dof(:sub, 2, [:A, :B])],
    [QN(cell=1, sub=1), QN(cell=1, sub=2)],
    [[0.0, 0.0], [0.0, a/√3]];
    vectors=[a1, a2]
)

# System DOFs: 1 cell × 2 sublattices × 2 spins → d = 4
dofs = SystemDofs([Dof(:cell, 1), Dof(:sub, 2, [:A, :B]), Dof(:spin, 2, [:up, :dn])])

nn_bonds     = bonds(unitcell, (:p, :p), 1)
onsite_bonds = bonds(unitcell, (:p, :p), 0)

# One-body hopping
onebody_hop = generate_onebody(dofs, nn_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -t : 0.0)

# k-grid for SCF
kpoints = build_kpoints([a1, a2], (100, 100))
Nk      = length(kpoints)

# Half-filling: 2 electrons per unit cell ⇒ 2*Nk total
n_elec = 2 * Nk

# ── DOF index map: (cell, sub, spin_idx) → linear index ──────────────────────
# valid_states ordering for [Dof(:cell,1), Dof(:sub,2), Dof(:spin,2)]:
#   [(cell=1,A,↑),(cell=1,B,↑),(cell=1,A,↓),(cell=1,B,↓)]
idx = Dict((qn[:cell], qn[:sub], qn[:spin]) => i for (i,qn) in enumerate(dofs.valid_states))

# ── Biased initial Green's functions ─────────────────────────────────────────
d = length(dofs.valid_states)

# AFM: A → ↑, B → ↓
function make_G_afm(Nk_local)
    G = zeros(ComplexF64, d, d, Nk_local)
    for ki in 1:Nk_local
        G[idx[(1,1,1)], idx[(1,1,1)], ki] = 1.0   # A, ↑
        G[idx[(1,2,2)], idx[(1,2,2)], ki] = 1.0   # B, ↓
    end
    return G
end

# PM: uniform half-filling (0.5 per spin per sublattice)
function make_G_pm(Nk_local)
    G = zeros(ComplexF64, d, d, Nk_local)
    for ki in 1:Nk_local
        G[idx[(1,1,1)], idx[(1,1,1)], ki] = 0.5
        G[idx[(1,2,1)], idx[(1,2,1)], ki] = 0.5
        G[idx[(1,1,2)], idx[(1,1,2)], ki] = 0.5
        G[idx[(1,2,2)], idx[(1,2,2)], ki] = 0.5
    end
    return G
end

# ── Order parameter and onsite potentials ────────────────────────────────────
function afm_order_and_densities(G_k)
    Nk_local = size(G_k, 3)
    G_loc = dropdims(sum(G_k, dims=3), dims=3) ./ Nk_local

    nA_up = real(G_loc[idx[(1,1,1)], idx[(1,1,1)]])
    nA_dn = real(G_loc[idx[(1,1,2)], idx[(1,1,2)]])
    nB_up = real(G_loc[idx[(1,2,1)], idx[(1,2,1)]])
    nB_dn = real(G_loc[idx[(1,2,2)], idx[(1,2,2)]])

    sA = (nA_up - nA_dn) / 2
    sB = (nB_up - nB_dn) / 2
    m_afm = abs((sA - sB) / 2)

    return m_afm, (nA_up, nA_dn, nB_up, nB_dn)
end

function build_on_site_onebody(nA_up, nA_dn, nB_up, nB_dn, U)
    # Hartree potentials for Hubbard U: ε_{iσ} = U * ⟨n_{i,-σ}⟩
    εA_up = U * nA_dn
    εA_dn = U * nA_up
    εB_up = U * nB_dn
    εB_dn = U * nB_up

    onsite = generate_onebody(dofs, onsite_bonds,
        (delta, qn1, qn2) -> begin
            qn1 == qn2 || return 0.0
            if qn1.sub == 1 && qn1.spin == 1
                return εA_up
            elseif qn1.sub == 1 && qn1.spin == 2
                return εA_dn
            elseif qn1.sub == 2 && qn1.spin == 1
                return εB_up
            elseif qn1.sub == 2 && qn1.spin == 2
                return εB_dn
            end
            return 0.0
        end;
        hc=false
    )

    onebody = (
        ops   = [onebody_hop.ops;   onsite.ops],
        delta = [onebody_hop.delta; onsite.delta],
        irvec = [onebody_hop.irvec; onsite.irvec]
    )
    return onebody
end

# ── Two-body Hubbard interaction ─────────────────────────────────────────────
# Note: for k=1 (onsite) generate_twobody only enumerates one spin ordering.
# The U term is U(n↑n↓ + n↓n↑), so we use the same convention as SDW_CDW.
function build_U_ops(U)
    return generate_twobody(dofs, onsite_bonds,
        (deltas, qn1, qn2, qn3, qn4) ->
            (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1,1,2,2) ? U : 0.0,
        order = (cdag, :i, c, :i, cdag, :i, c, :i))
end

# ── k-path for band structure ────────────────────────────────────────────────
A_mat = hcat(a1, a2)
B_mat = 2π * inv(A_mat)'
const b1 = B_mat[:, 1]
const b2 = B_mat[:, 2]

const Γ = [0.0, 0.0]
const K = (2b1 + b2) / 3
const M = (b1 + b2) / 2

kpath(p1, p2, n) = [p1 .+ t .* (p2 .- p1) for t in range(0.0, 1.0; length=n)]

nk = 120
k_ΓK = kpath(Γ, K, nk)
k_KM = kpath(K, M, nk)
k_MΓ = kpath(M, Γ, nk)
k_path = [k_ΓK; k_KM[2:end]; k_MΓ[2:end]]

# Arc-length x-axis
d_ΓK = norm(K - Γ)
d_KM = norm(M - K)
d_MΓ = norm(Γ - M)
arc = [collect(range(0,         d_ΓK;             length=nk));
       collect(range(d_ΓK,      d_ΓK+d_KM;        length=nk))[2:end];
       collect(range(d_ΓK+d_KM, d_ΓK+d_KM+d_MΓ;  length=nk))[2:end]]
xtick_pos = [0.0, d_ΓK, d_ΓK+d_KM, d_ΓK+d_KM+d_MΓ]
xtick_lab = ["Γ", "K", "M", "Γ"]

# ── Part 1: U=0 graphene band structure (2D) ─────────────────────────────────
println("=" ^ 60)
println("Part 1: U=0 graphene band structure (2D)")
println("=" ^ 60)

dofs_2d = SystemDofs([Dof(:cell, 1), Dof(:sub, 2, [:A, :B])])
nn_bonds_2d = bonds(unitcell, (:p, :p), 1)
onebody_2d = generate_onebody(dofs_2d, nn_bonds_2d, -t)
T_r_2d = build_Tr(dofs_2d, onebody_2d.ops, onebody_2d.irvec)
T_func_2d = build_Tk(T_r_2d)

bands_2d = [eigvals(Hermitian(T_func_2d(k))) for k in k_path]
bands_2d_mat = hcat(bands_2d...)

p0 = plot(; title="Graphene band structure (U=0)",
    xlabel="k-path", ylabel="E / t",
    xticks=(xtick_pos, xtick_lab),
    xlims=(0, d_ΓK+d_KM+d_MΓ), ylims=(-3.5, 3.5),
    legend=false, framestyle=:box)
hline!(p0, [0.0]; color=:gray, linestyle=:dash, linewidth=0.8)
for x in xtick_pos
    vline!(p0, [x]; color=:gray, linestyle=:dot, linewidth=0.8)
end
for n in axes(bands_2d_mat, 1)
    plot!(p0, arc, bands_2d_mat[n, :]; color=:steelblue, linewidth=1.5)
end
savefig(p0, joinpath(@__DIR__, "graphene_bands_2d.png"))
println("Saved: examples/SM_AFM/graphene_bands_2d.png")

# ── Part 2: Zigzag cylinder band structure (U=0) ─────────────────────────────
println("\n" * "=" ^ 60)
println("Part 2: Zigzag cylinder band structure (U=0)")
println("=" ^ 60)

Ny = 50
cylinder = Lattice(unitcell, (1, Ny))
dofs_cyl = SystemDofs([Dof(:cell, Ny), Dof(:sub, 2, [:A, :B])])
nn_bonds_cyl = bonds(cylinder, (:p, :o), 1)
onebody_cyl = generate_onebody(dofs_cyl, nn_bonds_cyl, -t)
T_r_cyl = build_Tr(dofs_cyl, onebody_cyl.ops, onebody_cyl.irvec)
T_func_cyl = build_Tk(T_r_cyl)

nkx = 200
kx_vals = range(-π/a, π/a; length=nkx)
k_cyl = [[kx, 0.0] for kx in kx_vals]

bands_cyl = [eigvals(Hermitian(T_func_cyl(k))) for k in k_cyl]
bands_cyl_mat = hcat(bands_cyl...)

p1 = plot(; title="Graphene zigzag cylinder (Ny=$Ny, U=0)",
    xlabel="kx", ylabel="E / t",
    xticks=([-π, -2π/3, 0, 2π/3, π], ["-π", "-2π/3", "0", "2π/3", "π"]),
    xlims=(-π, π), ylims=(-3.5, 3.5),
    legend=false, framestyle=:box)
hline!(p1, [0.0]; color=:gray, linestyle=:dash, linewidth=0.8)
vline!(p1, [-2π/3, 2π/3]; color=:gray, linestyle=:dot, linewidth=0.8)
for n in axes(bands_cyl_mat, 1)
    plot!(p1, kx_vals, bands_cyl_mat[n, :]; color=:steelblue, linewidth=0.6, alpha=0.7)
end
savefig(p1, joinpath(@__DIR__, "graphene_bands_cylinder.png"))
println("Saved: examples/SM_AFM/graphene_bands_cylinder.png")

# ── Part 3: Interactions, AFM transition, and gap opening ────────────────────
# ── Run U sweep (order parameter part) ───────────────────────────────────────
println("# Hubbard model on honeycomb lattice (half-filling)")
println("# k-grid: 100×100, Nk=$Nk")
println(@sprintf("# Expected mean-field Uc ≈ %.2f", Uc))
println()
println(@sprintf("# %-6s  %-12s  %-10s  %-10s  %s",
                 "U", "E_gs", "m_AF(gs)", "m_AF(afm)", "phase"))

results = Dict{Float64, Any}()

for U in unique(sort([U_sweep; U_bands]))
    U_ops = build_U_ops(U)
    twobody = (ops=U_ops.ops, delta=U_ops.delta, irvec=U_ops.irvec)

    r_afm = solve_hfk(dofs, onebody_hop, twobody, kpoints, n_elec;
        G_init=make_G_afm(Nk), n_restarts=1, tol=1e-12, verbose=false)
    r_pm  = solve_hfk(dofs, onebody_hop, twobody, kpoints, n_elec;
        G_init=make_G_pm(Nk),  n_restarts=1, tol=1e-12, verbose=false)

    E_afm = r_afm.energies.total
    E_pm  = r_pm.energies.total
    r_gs  = E_afm <= E_pm ? r_afm : r_pm
    phase = E_afm <= E_pm ? "AFM" : "PM"

    m_afm_gs, densities = afm_order_and_densities(r_gs.G_k)
    m_afm_afm, _ = afm_order_and_densities(r_afm.G_k)
    m_afm_pm,  _ = afm_order_and_densities(r_pm.G_k)
    results[U] = (r_gs=r_gs, m_afm_gs=m_afm_gs, m_afm_afm=m_afm_afm, m_afm_pm=m_afm_pm, densities=densities, phase=phase)

    println(@sprintf("  %-6.3f  %+12.6f  %-10.6f  %-10.6f  %s",
        U, r_gs.energies.total, m_afm_gs, m_afm_afm, phase))
end

open(joinpath(@__DIR__, "res.dat"), "w") do f
    println(f, "# U  m_AF(gs)  m_AF(afm)  phase")
    for U in U_sweep
        r = results[U]
        println(f, @sprintf("%.4f  %.8f  %.8f  %s",
                            U, r.m_afm_gs, r.m_afm_afm, r.phase))
    end
end

p_ord = plot(
    xlabel = "U / t",
    ylabel = "m_AF",
    title  = "Honeycomb Hubbard: AFM order at half-filling",
    legend = :topleft,
    framestyle = :box,
    size = (600, 400),
)
plot!(p_ord, U_sweep, [results[U].m_afm_afm for U in U_sweep];
    marker = :circle,
    lw = 2,
    color = :darkred,
    label = "AFM branch (m_AF)",
)
vline!(p_ord, [Uc];
    label = "Uc ≈ $(Uc)",
    linestyle = :dash,
    color = :gray,
    lw = 1,
)

ord_out = joinpath(@__DIR__, "afm_order_parameter.png")
savefig(p_ord, ord_out)
println("\nSaved: $ord_out")

# ── Plot: Mean-field bands for selected U ────────────────────────────────────

nb = length(U_bands)
rows = nb <= 3 ? 1 : 2
cols = ceil(Int, nb / rows)

band_plots = Vector{Plots.Plot}(undef, nb)

for (i, U) in enumerate(U_bands)
    r = results[U]
    U_ops = build_U_ops(U)
    twobody = (ops=U_ops.ops, delta=U_ops.delta, irvec=U_ops.irvec)
    bands_mat, _ = energy_bands(dofs, onebody_hop, twobody, kpoints, r.r_gs.G_k, k_path)

    p = plot(; title = "U/t = $(U)",
        xticks = (xtick_pos, xtick_lab),
        xlims = (0, d_ΓK+d_KM+d_MΓ),
        ylims = (-4.0+U/2, 4.0+U/2),
        legend = false, framestyle = :box)

    hline!(p, [0.0+U/2]; color=:gray, linestyle=:dash, linewidth=0.8)
    for x in xtick_pos
        vline!(p, [x]; color=:gray, linestyle=:dot, linewidth=0.8)
    end
    for n in axes(bands_mat, 1)
        plot!(p, arc, bands_mat[n, :]; color=:steelblue, linewidth=1.2)
    end

    band_plots[i] = p
end

p_bands = plot(band_plots...; layout=(rows, cols), size=(350*cols, 300*rows))

bands_out = joinpath(@__DIR__, "afm_bands.png")
savefig(p_bands, bands_out)
println("Saved: $bands_out")
