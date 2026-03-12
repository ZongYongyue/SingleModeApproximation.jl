using MeanFieldTheories
"""
Kane-Mele-Hubbard model in twisted MoTe2

H = t₁ Σ_{<ij>,σ} (c†_{iσ}c_{jσ} + h.c.)  +  |t₂|Σ_{<<ij>>,σ}(exp(im*ϕ*σ*νᵢⱼ)c†_{iσ}c_{jσ} + h.c.) +  U Σ_i n_{i↑}n_{i↓} 

"""


# ── Parameters ───────────────────────────────────────────────────────────────
const t1  = -3.49
const t2  = 0.81-1.32im
const t3  = 0.72
const t4  = -0.26
const t5  = 0.08

U = 100

# ── Lattice geometry ───────────────────────────────────────────────────────────

const a1 = [0.0, √3]
const a2 = [3/2, √3/2]

# Honeycomb unit cell: A at origin, B shifted by δ = (1.0, 0.0)
# Include a size-1 cell DOF in the unit cell
unitcell = Lattice(
    [Dof(:cell, 1), Dof(:sub, 2, [:A, :B])],
    [QN(cell=1, sub=1), QN(cell=1, sub=2)],
    [[0.0, 0.0], [1.0, 0.0]];
    vectors=[a1, a2]
)

# System DOFs: 1 cell × 2 sublattices × 2 spins → d = 4
dofs = SystemDofs([Dof(:cell, 1), Dof(:sub, 2, [:A, :B]), Dof(:spin, 2, [:up, :dn])])

n1_bonds     = bonds(unitcell, (:p, :p), 1)
n2_bonds     = bonds(unitcell, (:p, :p), 2)

onsite_bonds = bonds(unitcell, (:p, :p), 0)
