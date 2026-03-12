# Quantum System

The `quantumsystem` module provides a flexible, general-purpose framework for constructing quantum many-body Hamiltonians. It is organized around three levels that follow a natural physical hierarchy:

**Degrees of freedom → Lattice → Operators**

Each level is deliberately general: no restrictions are imposed on the number or nature of quantum numbers, the geometry of the lattice, or the form of operator terms. Physical constraints are enforced entirely through user-supplied functions, without modifying the library.

---

## 1. Degrees of Freedom

### Quantum numbers

A single-particle state is labelled by a `QuantumNumber` — a named-tuple wrapper that provides clean field access and display:

```julia
using MeanFieldTheories

qn = QN(site=2, spin=1, orbital=2)
qn.site      # 2
qn.spin      # 1
```

The `QN` alias is equivalent to `QuantumNumber`.

### Defining the Hilbert space

`SystemDofs` is built from a list of `Dof` objects, one per quantum-number component:

```julia
dofs = SystemDofs([
    Dof(:site, 4),
    Dof(:spin, 2),
    Dof(:orbital, 3)
])
```

`Dof(:name, n)` declares that this component takes integer values `1:n`. An optional label vector can be supplied for readability:

```julia
Dof(:spin, 2, [:up, :down])
```

The full set of valid single-particle states is the Cartesian product of all components:

```julia
dofs.valid_states   # Vector of QuantumNumber, length = 4 × 2 × 3 = 24
```

Any combination of components is supported — site, spin, orbital, sublattice, valley, layer, and so on — with no upper limit.

A `constraint` function can restrict the state space to a physically relevant subspace.
Without a constraint the total dimension is the full Cartesian product of all DOF sizes.
A constraint eliminates combinations that violate a physical rule, reducing the dimension:

```julia
# No constraint: 4 sites × 2 spins × 2 valleys = 16 states
dofs_full = SystemDofs([
    Dof(:site, 4),
    Dof(:spin,   2, [:up, :down]),
    Dof(:valley, 2, [:K, :Kprime])
])
length(dofs_full.valid_states)   # 16

# Spin-valley locking: ↑ locks to K, ↓ locks to K′ (spin index == valley index)
# Eliminates mixed combinations → 4 sites × 2 locked (spin, valley) pairs = 8 states
dofs_locked = SystemDofs([
    Dof(:site, 4),
    Dof(:spin,   2, [:up, :down]),
    Dof(:valley, 2, [:K, :Kprime])
], constraint = qn -> qn.spin == qn.valley)
length(dofs_locked.valid_states)  # 8
```

### Sorting and block structure

The order in which states appear in `valid_states` (and thus the row/column ordering of matrices) is controlled by `sortrule`, a permutation of `1:ndofs`:

```julia
# Default: reverse order, i.e. last DOF varies slowest
# [Dof 2 slow, Dof 1 fast] for a 2-DOF system
dofs = SystemDofs([Dof(:site, 3), Dof(:spin, 2)])
# valid_states: QN(site=1,spin=1), QN(site=2,spin=1), QN(site=3,spin=1),
#               QN(site=1,spin=2), QN(site=2,spin=2), QN(site=3,spin=2)
```

When a conserved quantum number makes the Hamiltonian block-diagonal, wrapping its index in a nested vector declares it as the block label:

```julia
# Spin is conserved → two blocks: spin=1 (indices 1:3) and spin=2 (indices 4:6)
dofs = SystemDofs([Dof(:site, 3), Dof(:spin, 2)], sortrule = [[2], 1])
dofs.blocks   # [1:3, 4:6]
```

Blocks are exploited automatically by the mean-field solvers to reduce memory and computation: the interaction matrix only needs to store elements within each block, giving a factor of $B^2$ reduction for $B$ symmetry blocks (e.g., 75% reduction for two spin blocks).

Multiple nested levels can be specified for simultaneous conservation of several quantum numbers:

```julia
# Block by valley first, then by spin within each valley block
dofs = SystemDofs([
    Dof(:site, 4), Dof(:orbital, 2), Dof(:spin, 2), Dof(:valley, 2)
], sortrule = [[4, 3], 2, 1])
# → 4 blocks: (K,↑), (K,↓), (K′,↑), (K′,↓)
```

### Index conversion

Linear indices map quantum numbers to matrix row/column positions. The lookup is O(1) via an internal hash table:

```julia
idx = qn2linear(dofs, QN(site=2, spin=1))   # Int
qn  = linear2qn(dofs, idx)                  # QuantumNumber
```

---

## 2. Lattice

### Unit cell

A `Lattice` maps position states to real-space coordinates. It is constructed from position degrees of freedom, the list of position states (one per site in the unit cell), and their Cartesian coordinates:

```julia
# Single-site square-lattice unit cell
unitcell = Lattice(
    [Dof(:cell, 1)],    # first DOF must have size=1 for tiling
    [QN(cell=1)],
    [[0.0, 0.0]];
    vectors=[[1.0, 0.0], [0.0, 1.0]]
)

# Two-site (honeycomb) unit cell
unitcell = Lattice(
    [Dof(:cell, 1), Dof(:sublattice, 2, [:A, :B])],
    [QN(cell=1, sublattice=1), QN(cell=1, sublattice=2)],
    [[0.0, 0.0], [0.5, 0.289]];
    vectors=[[1.0, 0.0], [0.5, 0.866]]
)
```

The position DOFs of a `Lattice` can be a subset of the full system DOFs: `SystemDofs` may include additional internal DOFs (spin, orbital, ...) that the lattice does not know about. This separation means spatial and internal structure can be defined independently.

### Supercell construction

The physical simulation cell is built by tiling the unit cell along primitive lattice vectors. The first DOF of the unit cell (which must have `size=1`) is expanded to the total number of unit cells:

```julia
a1, a2 = [1.0, 0.0], [0.0, 1.0]
lattice = Lattice(unitcell, (4, 4))
# First DOF now has size = 16 (= 4×4 unit cells)
# vectors = [4*a1, 4*a2] are set automatically
```

The `vectors` field must be set in the unit cell; it is used for periodic boundary conditions in bond generation and for the reciprocal-space grid in momentum-space HF.

Site coordinates can be queried by quantum number, even when the QN contains extra components beyond the position DOFs:

```julia
get_coordinate(lattice, QN(cell=3, sublattice=2))
get_coordinate(lattice, QN(cell=3, sublattice=2, spin=1))  # extra keys are ignored
```

### Bonds

A `Bond` connects one or more sites and stores, for each site:

- `states`: the position QN (site label within the unit cell),
- `coordinates`: absolute (unwrapped) Cartesian coordinate,
- `icoordinates`: lattice vector of the unit cell containing that site.

The `icoordinates` field is essential for momentum-space calculations: the difference `icoordinates[2] - icoordinates[1]` gives the integer lattice vector crossed by the bond, which becomes the displacement $\mathbf{R}$ in $T_{\mathbf{R}}$. For bonds that cross a periodic boundary, the periodic image site carries a non-zero `icoordinate` shift (e.g., `[-1.0, 0.0]`), while the originating site has `[0.0, 0.0]`.

`bonds(lattice, boundary, neighbors)` generates all bonds of a given neighbor shell. The boundary tuple specifies periodic (`:p`) or open (`:o`) conditions per spatial direction:

```julia
# Onsite bonds (one site per bond)
onsite = bonds(lattice, (:p, :p), 0)

# Nearest-neighbor bonds (two sites per bond), periodic in both directions
nn = bonds(lattice, (:p, :p), 1)

# Next-nearest-neighbor, open in y
nnn = bonds(lattice, (:p, :o), 2)

# Multiple shells at once
nn_and_nnn = bonds(lattice, (:p, :p), [1, 2])

# Specific distances
specific = bonds(lattice, (:p, :p), [1.0, 1.732])
```

Bonds are not restricted to two sites. Three-site and four-site bonds can be constructed manually (`Bond(states, coordinates, icoordinates)`) and passed to `generate_twobody` with the appropriate `order` to model ring-exchange and other multi-site interactions.

---

## 3. Operators

### Primitive operators and operator terms

Single creation and annihilation operators are constructed from quantum numbers:

```julia
i = QN(site=1, spin=1)
j = QN(site=2, spin=1)

cdag(i)   # creation  c†_{site=1, spin=1}
c(j)      # annihilation c_{site=2, spin=1}
```

An `Operators` object is a product of primitive operators with a (possibly complex) coefficient:

```julia
# Hopping term: -t c†_i c_j
Operators(-1.0, [cdag(i), c(j)])

# Complex hopping
Operators(0.81 - 1.32im, [cdag(i), c(j)])

# Four-operator interaction: U c†_↑ c_↑ c†_↓ c_↓
Operators(U, [cdag(QN(site=1,spin=1)), c(QN(site=1,spin=1)),
              cdag(QN(site=1,spin=2)), c(QN(site=1,spin=2))])
```

Operators are stored in the user-supplied order; reordering to creation-annihilation alternating order (c†c c†c …) happens automatically during matrix construction.

### Generating one-body terms

`generate_onebody` iterates over all bonds and all combinations of internal DOFs, calling a user-supplied `value` function to assign the coefficient for each term:

```julia
onebody = generate_onebody(dofs, bonds, value; order=(cdag, 1, c, 2), hc=true)
```

**Arguments:**
- `value`: a `Number` (uniform coefficient) or a function `(delta, qn1, qn2) -> Number`.
  `delta` is generated as `bond.coordinates[1] - bond.coordinates[2]` — the vector
  from site 2 to site 1, encoding the full bond direction (magnitude and angle).
- `order`: `(op_type1, site1, op_type2, site2)` — which site of the bond each operator
  acts on. Default `(cdag, 1, c, 2)` places the creation operator on site 1 and the
  annihilation operator on site 2, giving the standard hopping $c^\dagger_1 c_2$.
- `hc=true`: automatically appends the Hermitian conjugate of each term with the
  displacement sign flipped. Together with a `value` function that returns real
  coefficients, this guarantees a Hermitian Hamiltonian without manual bookkeeping.

**Return value:** a NamedTuple `(ops, delta, irvec)` — three parallel vectors.
`irvec[n] = icoordinates[2] - icoordinates[1]` for term `n` is the unit-cell
displacement, which is passed directly to `build_Tr` for momentum-space HF.

**Example — spin-conserving nearest-neighbor hopping:**

```julia
t_onebody = generate_onebody(dofs, nn_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -1.0 : 0.0)
# t_onebody is a NamedTuple with three parallel vectors:
#   t_onebody.ops   → Vector{Operators}        (one entry per generated term; hc=true doubles the count)
#   t_onebody.delta → Vector{SVector{D,T}}     (physical bond vector for each term)
#   t_onebody.irvec → Vector{SVector{D,T}}     (unit-cell displacement R for each term; used by build_Tr)
t_ops = t_onebody.ops
```

**Example — direction-dependent complex hopping:**

The `delta` vector carries the full bond direction, enabling valley- or direction-dependent phases without any separate bookkeeping:

```julia
t2_onebody = generate_onebody(dofs, nnn_bonds,
    (delta, qn1, qn2) -> begin
        qn1.spin != qn2.spin && return 0.0          # spin-conserving
        tau = qn1.spin == 1 ? +1 : -1               # valley sign
        nu  = sign(delta[1] + 0.5*delta[2]) |> Int  # circulation direction
        return (0.81 + 1.32im) * (tau * nu)
    end)
# t2_onebody.ops   → Vector{Operators} with complex coefficients
# t2_onebody.delta → Vector{SVector{D,T}} encoding bond direction (used above for nu)
# t2_onebody.irvec → Vector{SVector{D,T}} encoding lattice vector R
t2_ops = t2_onebody.ops
```

**Example — orbital-selective intra- and inter-orbital hopping:**

```julia
# Intra-orbital: different hopping per orbital, same spin
t_intra = generate_onebody(dofs, nn_bonds,
    (delta, qn1, qn2) -> begin
        qn1.spin != qn2.spin && return 0.0
        qn1.orbital != qn2.orbital && return 0.0
        qn1.orbital == 1 ? -1.0 : -0.8
    end).ops
# .ops → Vector{Operators}: one term per (bond, spin, orbital) combination passing the filter

# Inter-orbital: mixing between orbitals
t_inter = generate_onebody(dofs, nn_bonds,
    (delta, qn1, qn2) ->
        qn1.spin == qn2.spin && qn1.orbital != qn2.orbital ? -0.2 : 0.0).ops
```

### Generating two-body terms

`generate_twobody` follows the same pattern for four-operator interaction terms, but offers additional generality:

```julia
twobody = generate_twobody(dofs, bonds, value; order=(cdag, :i, c, :i, cdag, :j, c, :j))
```

**Arguments:**
- `value`: a `Number` or a function `(deltas, qn1, qn2, qn3, qn4) -> Number`.
  `deltas[m] = coord(op_m) - coord(op_4)` for `m = 1,2,3` (always 3 entries).
  Values depend on the current site assignment and can be used to filter by direction.
- `order`: `(type1, label1, type2, label2, type3, label3, type4, label4)` where labels are **Symbols**
  acting as free Einstein summation indices. Each unique Symbol is an independent position index
  summed over all **injective assignments** to bond sites (different Symbols → different sites).
  The number of unique Symbols must not exceed the number of sites in each bond.

**Return value:** a NamedTuple `(ops, delta, irvec)`.
- `ops`: operators already reordered to creation-annihilation alternating order (c†c c†c) with the fermionic sign absorbed into the coefficient.
- `delta`: a `Vector{NTuple{3, SVector{D,T}}}`, where each tuple `(δ1, δ2, δ3)`
  gives the three physical (Cartesian) displacements after reordering:
  $\delta_n = \mathrm{coord}(\mathrm{op}_n) - \mathrm{coord}(\mathrm{op}_4)$.
- `irvec`: a `Vector{NTuple{3, SVector{D,T}}}`, where each tuple `(τ1, τ2, τ3)`
  gives the three unit-cell displacements: $\tau_n = \mathrm{icoord}(\mathrm{op}_n) - \mathrm{icoord}(\mathrm{op}_4)$.

**Example — nearest-neighbor Coulomb $V \sum_{i \neq j, \sigma\sigma'} n_{i\sigma} n_{j\sigma'}$:**

```julia
# Default order=(cdag,:i,c,:i,cdag,:j,c,:j): :i and :j independently loop over bond sites,
# generating both n_A n_B and n_B n_A for each bond (full ordered-pair sum).
nn_coulomb = generate_twobody(dofs, nn_bonds,
    (deltas, qn1, qn2, qn3, qn4) ->
        qn1.spin == qn2.spin && qn3.spin == qn4.spin ? V : 0.0)
```

**Example — onsite Hubbard $U \sum_i n_{i\uparrow} n_{i\downarrow}$:**

```julia
# order=(cdag,:i,c,:i,cdag,:i,c,:i): single label :i, all operators at the same site.
# For 1-site bonds (onsite_bonds), only one assignment exists.
hubbard = generate_twobody(dofs, onsite_bonds,
    (deltas, qn1, qn2, qn3, qn4) ->
        (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1, 1, 2, 2) ? U : 0.0,
    order = (cdag, :i, c, :i, cdag, :i, c, :i))
```

**Example — Hund's pair hopping $J \sum_{i \neq j} c^\dagger_{i\uparrow} c^\dagger_{i\downarrow} c_{j\downarrow} c_{j\uparrow}$:**

```julia
# order=(cdag,:i,cdag,:i,c,:j,c,:j): pair created at :i, annihilated at :j.
hund = generate_twobody(dofs, nn_bonds,
    (deltas, qn1, qn2, qn3, qn4) ->
        (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1, 2, 2, 1) ? J : 0.0,
    order = (cdag, :i, cdag, :i, c, :j, c, :j))
```

**Example — NN Coulomb along x only (filter by direction):**

```julia
# deltas[1] = coord(:i site) - coord(:j site); filter keeps only |Δx|=1, Δy=0.
nn_x = generate_twobody(dofs, nn_bonds,
    (deltas, qn1, qn2, qn3, qn4) ->
        abs(deltas[1][1]) ≈ 1.0 && iszero(deltas[1][2]) ? V : 0.0)
```

The `value` function mechanism makes all these variants expressible without any modifications to the library. Physical constraints (spin conservation, orbital selection, direction filtering) are encoded entirely in user-supplied anonymous functions.

### Building matrices

Once operator lists are assembled, matrix representations are built by:

```julia
# Full one-body Hamiltonian H[i,j]
H = build_onebody_matrix(dofs, onebody.ops)

# Full two-body interaction tensor V[i,j,k,l] in creation-annihilation alternating order
V = build_interaction_tensor(dofs, twobody.ops)
```

For mean-field calculations the interaction is more efficiently accessed through `build_U` (real-space HF) or `build_Vr` / `build_Vk` (momentum-space HF), which bypass the $O(N^4)$ dense tensor entirely. See the respective HF solver documentation for details.

---

## Summary

| Component | Key design choices |
|-----------|-------------------|
| `SystemDofs` | Arbitrary DOF names and sizes; `constraint` prunes the state space; nested `sortrule` creates block-diagonal structure for conserved quantum numbers |
| `Lattice` / `Bond` | Any unit-cell geometry; tiling constructor handles PBC automatically; `icoordinates` carries unit-cell origin information needed for Bloch-space calculations |
| `generate_onebody` | `delta`-aware value function for direction-dependent coefficients; `hc=true` for automatic Hermitian conjugate; returns `(ops, delta, irvec)` |
| `generate_twobody` | 1–4 site bonds; `deltas` vector for full cluster geometry; arbitrary operator ordering; returns `(ops, delta, irvec)` where operators are in creation-annihilation alternating order, `delta` holds physical position differences `(δ1,δ2,δ3)`, and `irvec` holds unit-cell displacements `(τ1,τ2,τ3)`; all vectors are `SVector{D,T}` |
| `build_*` | Sparse matrix/tensor construction; specialized paths in HF solvers avoid dense $N^4$ tensors |
