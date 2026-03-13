using MeanFieldTheories: Lattice as MLattice, bonds as mbonds
using MeanFieldTheories
using LinearAlgebra
using Plots
"""
Kane-Mele-Hubbard model in twisted MoTe2

H = t₁ Σ_{<ij>,σ} (c†_{iσ}c_{jσ} + h.c.)  +  |t₂|Σ_{<<ij>>,σ}(exp(im*ϕ*σ*νᵢⱼ)c†_{iσ}c_{jσ} + h.c.) +  U Σ_i n_{i↑}n_{i↓} 

"""


# ── Parameters ───────────────────────────────────────────────────────────────
# const t1  = -3.49
# const t2  = 0.81-1.32im
# const t3  = 0.72
# const t4  = -0.26
# const t5  = 0.08

const t1  = -2.5
const t2  = -0.4 + 0.65*im
const t3  = 0.4
const t4  = -0.1
const t5  = 0.02

U = 100

# ── Lattice geometry ───────────────────────────────────────────────────────────

const a1 = [0.0, √3]
const a2 = [3/2, √3/2]

# Honeycomb unit cell: A at origin, B shifted by δ = (1.0, 0.0)
# Include a size-1 cell DOF in the unit cell
unitcell = MLattice(
    [Dof(:cell, 1), Dof(:sub, 2, [:A, :B])],
    [QN(cell=1, sub=1), QN(cell=1, sub=2)],
    [[0.0, 0.0], [1.0, 0.0]];
    vectors=[a1, a2]
)

# System DOFs: 1 cell × 2 sublattices × 2 spins → d = 4
dofs = SystemDofs([Dof(:cell, 1), Dof(:sub, 2, [:A, :B]), Dof(:spin, 2, [:up, :dn])])

n1_bonds = mbonds(unitcell, (:p, :p), 1)
n2_bonds = mbonds(unitcell, (:p, :p), 2)
n3_bonds = mbonds(unitcell, (:p, :p), 3)
n4_bonds = mbonds(unitcell, (:p, :p), 4)
n5_bonds = mbonds(unitcell, (:p, :p), 5)

on_bonds = mbonds(unitcell, (:p, :p), 0)

# lattice = Lattice(unitcell, (9, 9))
# n1_bonds     = bonds(lattice, (:p, :p), 1)
# n2_bonds     = bonds(lattice, (:p, :p), 2)
# dofs = SystemDofs([Dof(:cell, 81), Dof(:sub, 2, [:A, :B]), Dof(:spin, 2, [:up, :dn])])

onebody_t1 = generate_onebody(dofs, n1_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? t1 : 0.0)

onebody_t2 = generate_onebody(dofs, n2_bonds,
    (delta, qn1, qn2) -> begin 
    t, phi = abs(t2), angle(t2)
    if (qn1.spin !== qn2.spin)
        return 0.0
    else
        sigma = (qn1.spin == qn2.spin == 1) ? 1.0 : -1.0 
        if (qn1.sub == qn2.sub == 1)
            nu = any(x -> ≈(delta, x; atol=1e-6), [[3/2, √3/2], [-3/2, √3/2], [0, -√3]]) ? 1.0 : any(x -> ≈(delta, x; atol=1e-6), [[-3/2, -√3/2], [3/2, -√3/2], [0, √3]]) ? -1.0 : throw(ArgumentError("Invalid nn Bond with distance=$(delta) on sub=1"))
        elseif (qn1.sub == qn1.sub == 2)
            nu = any(x -> ≈(delta, x; atol=1e-6), [[3/2, -√3/2], [0, √3], [-3/2, -√3/2]]) ? 1.0 : any(x -> ≈(delta, x; atol=1e-6), [[-3/2, √3/2], [0, -√3], [3/2, √3/2]]) ? -1.0 : throw(ArgumentError("Invalid nn Bond with distance=$(delta) on sub=2"))
        else
            return 0.0
        end
        return t*exp(im*sigma*nu*phi)
    end
end)

onebody_t3 = generate_onebody(dofs, n3_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? t3 : 0.0)

onebody_t4 = generate_onebody(dofs, n4_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? t4 : 0.0)

onebody_t5 = generate_onebody(dofs, n5_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? t5 : 0.0)

function build_U_ops(U)
    return generate_twobody(dofs, on_bonds,
        (deltas, qn1, qn2, qn3, qn4) ->
            (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1,1,2,2) ? U : 0.0,
        order = (cdag, :i, c, :i, cdag, :i, c, :i))
end

# ── Merge one-body terms ─────────────────────────────────────────────────────
function merge_onebody(parts...)
    return (
        ops   = vcat((p.ops   for p in parts)...),
        delta = vcat((p.delta for p in parts)...),
        irvec = vcat((p.irvec for p in parts)...),
    )
end

onebody_hop = merge_onebody(onebody_t1, onebody_t2, onebody_t3, onebody_t4, onebody_t5)

# ── U=0 self-consistent solution on uniform k-grid ───────────────────────────
U = 0.0
twobody = build_U_ops(U)

Nk1, Nk2 = 80, 80
kpoints = build_kpoints([a1, a2], (Nk1, Nk2))
n_elec  = 2 * length(kpoints)  # half-filling: 2 electrons per unit cell

r0 = solve_hfk(dofs, onebody_hop, twobody, kpoints, n_elec;
    n_restarts=1, tol=1e-12, verbose=false)

# ── k-path Γ–K1–K2–Γ for band structure ─────────────────────────────────────
A_mat = hcat(a1, a2)
B_mat = 2π * inv(A_mat)'
const b1 = B_mat[:, 1]
const b2 = B_mat[:, 2]

const γ   = [0.0, 0.0]
const k_p = (2b1 + b2) / 3
const k_m = (b1 + 2b2) / 3

kpath(p1, p2, n) = [p1 .+ t .* (p2 .- p1) for t in range(0.0, 1.0; length=n)]

nk = 120
k_ΓK1 = kpath(γ,  k_p, nk)
k_K1K2 = kpath(k_p, k_m, nk)
k_K2Γ = kpath(k_m, γ, nk)
k_path = [k_ΓK1; k_K1K2[2:end]; k_K2Γ[2:end]]

d_ΓK1 = norm(k_p - γ)
d_K1K2 = norm(k_m - k_p)
d_K2Γ = norm(γ - k_m)
arc = [collect(range(0,         d_ΓK1;               length=nk));
       collect(range(d_ΓK1,      d_ΓK1+d_K1K2;       length=nk))[2:end];
       collect(range(d_ΓK1+d_K1K2, d_ΓK1+d_K1K2+d_K2Γ; length=nk))[2:end]]
xtick_pos = [0.0, d_ΓK1, d_ΓK1+d_K1K2, d_ΓK1+d_K1K2+d_K2Γ]
xtick_lab = ["γ", "k₊", "k₋", "γ"]

# ── Bands at U=0 from converged G_k ─────────────────────────────────────────
bands_mat, _ = energy_bands(dofs, onebody_hop, twobody, kpoints, r0.G_k, k_path)

p = plot(; title="Kane-Mele-Hubbard (U=0)",
    xlabel="k-path", ylabel="E",
    xticks=(xtick_pos, xtick_lab),
    xlims=(0, d_ΓK1+d_K1K2+d_K2Γ),
    legend=false, framestyle=:box)

hline!(p, [0.0]; color=:gray, linestyle=:dash, linewidth=0.8)
for x in xtick_pos
    vline!(p, [x]; color=:gray, linestyle=:dot, linewidth=0.8)
end
for n in axes(bands_mat, 1)
    plot!(p, arc, bands_mat[n, :]; color=:steelblue, linewidth=1.5)
end

out = joinpath(@__DIR__, "qahi_bands_U0.png")
savefig(p, out)
println("Saved: $out")
