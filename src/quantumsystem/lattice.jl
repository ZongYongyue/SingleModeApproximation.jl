"""
Lattice Structures for Quantum Systems

Provides types and functions for defining lattice structures,
coordinate mappings, and bond generation for many-body systems.
"""

#==================== Lattice Specification ====================#

"""
    Lattice{D, Q, T}

Represents a lattice structure mapping position states to real-space coordinates.

# Type Parameters
- `D`: Spatial dimension (e.g., 2 for 2D, 3 for 3D)
- `Q`: QuantumNumber type for position states
- `T`: Coordinate element type (subtype of Real)

# Fields
- `position_dofs::Vector{Dof}`: Degrees of freedom that define positions
- `position_states::Vector{Q}`: All valid position states
- `coordinates::Vector{SVector{D,T}}`: Coordinate for each position state (length D)
- `vectors::Union{Nothing, Vector{SVector{D,T}}}`: Supercell vectors for PBC (optional)

# Notes
- `position_states[i]` corresponds to `coordinates[i]`
- `vectors` is set automatically when using tiling constructor

# Examples
```julia
# 2D honeycomb lattice with 16 unit cells and 2 sublattices
position_dofs = [Dof(:moire_cell, 16), Dof(:sublattice, 2, [:A, :B])]
position_states = [QN(moire_cell=1, sublattice=1), QN(moire_cell=1, sublattice=2), ...]
coordinates = [[0.0, 0.0], [1.0, 0.0], ...]

lattice = Lattice(position_dofs, position_states, coordinates)
```
"""
struct Lattice{D, Q<:QuantumNumber, T<:Real}
    position_dofs::Vector{Dof}
    position_states::Vector{Q}
    coordinates::Vector{SVector{D,T}}
    vectors::Union{Nothing, Vector{SVector{D,T}}}
end

# Public constructor: accepts Vector{Vector{T}} and converts to SVector (backward compatible)
function Lattice(
    position_dofs::Vector{Dof},
    position_states::Vector{Q},
    coordinates::Vector{Vector{T}};
    vectors::Union{Nothing, Vector{Vector{T}}} = nothing
) where {Q<:QuantumNumber, T<:Real}
    D = length(coordinates[1])
    sv_coords    = SVector{D,T}.(coordinates)
    sv_supercell = isnothing(vectors) ? nothing : SVector{D,T}.(vectors)
    Lattice{D,Q,T}(position_dofs, position_states, sv_coords, sv_supercell)
end

# Internal constructor: already-converted SVectors (used by tiling constructor)
function Lattice(
    position_dofs::Vector{Dof},
    position_states::Vector{Q},
    coordinates::Vector{SVector{D,T}};
    vectors::Union{Nothing, Vector{SVector{D,T}}} = nothing
) where {D, Q<:QuantumNumber, T<:Real}
    Lattice{D,Q,T}(position_dofs, position_states, coordinates, vectors)
end

"""
    get_coordinate(lattice::Lattice, state) -> Vector{T}

Return the coordinate for a given quantum state. The input can be QuantumNumber
or NamedTuple, and can contain extra keys beyond the position DOFs.

# Examples
```julia
lattice = Lattice(...)  # with position_dofs = [:moire_cell, :sublattice]

# Exact position state
coord = get_coordinate(lattice, QN(moire_cell=1, sublattice=2))

# Full quantum state with extra DOFs
coord = get_coordinate(lattice, QN(moire_cell=1, sublattice=2, spin=1, valley=1))
```
"""
function get_coordinate(lattice::Lattice, state::QuantumNumber)
    # Extract position-related keys from the input state
    pos_keys = Tuple(dof.name for dof in lattice.position_dofs)
    pos_nt = NamedTuple{pos_keys}(Tuple(state[k] for k in pos_keys))
    pos_state = QuantumNumber(pos_nt)

    idx = findfirst(==(pos_state), lattice.position_states)
    if isnothing(idx)
        error("Position state $pos_state not found in lattice")
    end
    return lattice.coordinates[idx]
end

# Also accept NamedTuple
function get_coordinate(lattice::Lattice, state::NamedTuple)
    return get_coordinate(lattice, QuantumNumber(state))
end

"""
    Base.show(io::IO, lattice::Lattice)

Custom display for Lattice in REPL.
"""
function Base.show(io::IO, lattice::Lattice{D, Q, T}) where {D, Q, T}
    n = length(lattice.position_dofs)
    nstates = length(lattice.position_states)

    println(io, "Lattice{$D} with $n position degree(s) of freedom:")
    for (i, dof) in enumerate(lattice.position_dofs)
        print(io, "  $i. $(dof.name): $(dof.size) states")
        if dof.labels !== nothing
            print(io, " [", join(dof.labels, ", "), "]")
        end
        println(io)
    end
    println(io, "Number of position states: $nstates")
    if isnothing(lattice.vectors)
        print(io, "Supercell vectors: not set")
    else
        print(io, "Supercell vectors: set")
    end
end

"""
    Lattice(unitcell, box_size)

Create a lattice by tiling a unit cell.

The first DOF of unitcell defines the unit cell type (with size=1 in unitcell).
After tiling, its size becomes Nx*Ny*... (total number of unit cells).
The unitcell must define `vectors` which are used as the unit-cell
translation vectors; the new lattice uses `box_size .* unitcell.vectors`.

# Arguments
- `unitcell::Lattice`: The unit cell definition (first DOF has size=1)
- `box_size::NTuple{D, Int}`: Number of unit cells in each direction (Nx, Ny, ...)

# Returns
A new `Lattice` with the first DOF's size expanded to prod(box_size).
Coordinates: R = (i1-1)*a1 + (i2-1)*a2 + ... + r_unitcell

# Examples
```julia
# Define honeycomb unit cell (first DOF = cell type with size 1)
unitcell = Lattice(
    [Dof(:honeycomb_cell, 1), Dof(:sublattice, 2, [:A, :B])],
    [QN(honeycomb_cell=1, sublattice=1), QN(honeycomb_cell=1, sublattice=2)],
    [[0.0, 0.0], [0.5, 0.289]]
)

# Tile with unit-cell vectors already stored in unitcell.vectors
lattice = Lattice(unitcell, (4, 4))
# Result: DOFs [honeycomb_cell (16), sublattice (2)], 32 sites total
# vectors = [4*a1, 4*a2] automatically set
```
"""
function Lattice(
    unitcell::Lattice{D, Q, T},
    box_size::NTuple{M, Int}
) where {D, Q, T, M}
    @assert M == D "box_size dimension ($M) must match lattice dimension ($D)"
    @assert unitcell.position_dofs[1].size == 1 "First DOF of unitcell must have size=1"
    isnothing(unitcell.vectors) && error("unitcell.vectors is not set. Please set it in unitcell before tiling.")

    unitcell_vectors = unitcell.vectors

    # Number of unit cells
    n_cells = prod(box_size)

    # Compute supercell vectors (convert user-provided Vector{Vector{T}} to SVector)
    vectors = SVector{D,T}[rd(SVector{D,T}(box_size[i] .* unitcell_vectors[i])) for i in 1:D]

    # Create new DOFs: expand first DOF's size, keep others
    first_dof = unitcell.position_dofs[1]
    new_dofs = [Dof(first_dof.name, n_cells, first_dof.labels); unitcell.position_dofs[2:end]]

    # Get names from unitcell (same structure, just different values for first DOF)
    all_names = Tuple(dof.name for dof in unitcell.position_dofs)

    # Generate position_states and coordinates
    N = QuantumNumber{all_names, NTuple{length(all_names), Int}}
    position_states = N[]
    sv_coordinates  = SVector{D,T}[]

    # Iterate over all unit cells (linear index)
    cell_idx = 1
    for indices in Iterators.product([1:s for s in box_size]...)
        # Calculate unit cell origin: R0 = (i1-1)*a1 + (i2-1)*a2 + ...
        R0 = zero(SVector{D,T})
        for (dim, i) in enumerate(indices)
            R0 = R0 + (i - 1) * unitcell_vectors[dim]
        end
        R0 = rd(R0)

        # Add each site in the unit cell
        for (cell_state, cell_coord) in zip(unitcell.position_states, unitcell.coordinates)
            # Replace first DOF value with cell_idx, keep rest
            state_values = values(cell_state)
            new_nt = NamedTuple{all_names}((cell_idx, state_values[2:end]...))
            push!(position_states, QuantumNumber(new_nt))

            # Calculate coordinate
            new_coord = rd(R0 + cell_coord)
            push!(sv_coordinates, new_coord)
        end

        cell_idx += 1
    end

    return Lattice(new_dofs, position_states, sv_coordinates; vectors=vectors)
end

#==================== Bond Specification ====================#

"""
    Bond{Q, T}

Represents a bond connecting one or more sites.

# Fields
- `states::Vector{Q}`: Position states (length determines bond order: 1=onsite, 2=two-body, etc.)
- `coordinates::Vector{SVector{D,T}}`: Physical (unwrapped) coordinates for each site.
  For periodic bonds these may lie outside the simulation cell.
- `icoordinates::Vector{SVector{D,T}}`: Unit-cell origin for each site.
  `icoordinates[n]` is the lattice vector of the unit cell that contains `coordinates[n]`.
  For intra-cell bonds all entries are zero vectors. For bonds that cross a periodic
  boundary, the entry for the image site carries the lattice shift, e.g. `[-2.0, 0.0]`.

# Examples
```julia
# Onsite bond — icoordinates default to zero vectors
bond = Bond([state1], [coord1])

# Two-body intra-cell bond
bond = Bond([state1, state2], [coord1, coord2])

# Two-body bond crossing a periodic boundary (site 2 lives in cell [-2, 0])
bond = Bond([state1, state2], [coord1, coord2], [[0.0, 0.0], [-2.0, 0.0]])
```
"""
struct Bond{Q<:QuantumNumber, T<:Real, D}
    states::Vector{Q}
    coordinates::Vector{SVector{D,T}}
    icoordinates::Vector{SVector{D,T}}
end

# Public constructors: accept Vector{Vector{T}} and convert to SVector (backward compatible)
function Bond(
    states::Vector{Q},
    coordinates::Vector{Vector{T}},
    icoordinates::Vector{Vector{T}}
) where {Q<:QuantumNumber, T<:Real}
    D = length(coordinates[1])
    Bond{Q,T,D}(states, SVector{D,T}.(coordinates), SVector{D,T}.(icoordinates))
end

function Bond(states::Vector{Q}, coordinates::Vector{Vector{T}}) where {Q<:QuantumNumber, T<:Real}
    D = length(coordinates[1])
    sv = SVector{D,T}.(coordinates)
    Bond{Q,T,D}(states, sv, [zero(SVector{D,T}) for _ in states])
end

function Base.show(io::IO, bond::Bond)
    _show_bond_plain(io, bond)
end

function Base.show(io::IO, ::MIME"text/plain", bond::Bond)
    _show_bond_plain(io, bond)
end

function _show_bond_plain(io::IO, bond::Bond)
    states_str = "[" * join(string.(bond.states), ", ") * "]"
    coords = [collect(v) for v in bond.coordinates]
    icoords = [collect(v) for v in bond.icoordinates]
    print(io, "Bond(", states_str, ", ", coords, ", ", icoords, ")")
end

"""
    bonds(lattice, boundary, neighbors)

Generate bonds based on neighbor specification.

# Arguments
- `lattice::Lattice`: The lattice structure (must have vectors set)
- `boundary::NTuple{D, Symbol}`: Boundary conditions, :p (periodic) or :o (open)
- `neighbors`: Neighbor specification:
  - `Int`: n-th nearest neighbor (0=onsite, 1=nearest, 2=next-nearest, ...)
  - `Vector{Int}`: multiple neighbor orders [n1, n2, ...]
  - `Vector{<:Real}`: specific distances [r1, r2, ...]

# Returns
`Vector{Bond}` with appropriate number of sites per bond

# Examples
```julia
# Create lattice with tiling (vectors auto-set)
lattice = Lattice(unitcell, (4, 4))

# Onsite bonds - returns Vector{Bond} with 1 site per bond
onsite = bonds(lattice, (:p, :p), 0)
# onsite[1] displays as: Bond{1}(QN(cell=1, sub=1) @ [0.000, 0.000])

# Nearest neighbor bonds - returns Vector{Bond} with 2 sites per bond
nn = bonds(lattice, (:p, :o), 1)
# nn[1] displays as: Bond{2}(QN(cell=1, sub=1) @ [0.000, 0.000], QN(cell=1, sub=2) @ [0.500, 0.289])

# Multiple neighbor orders
nn_and_nnn = bonds(lattice, (:p, :p), [1, 2])

# Specific distances
specific = bonds(lattice, (:p, :p), [1.0, 1.732])
```
"""
function bonds(
    lattice::Lattice{D, Q, T},
    boundary::NTuple{D, Symbol},
    neighbors::Int
) where {D, Q, T}
    if neighbors == 0
        return _generate_onsite_bonds(lattice)
    else
        vectors = _get_vectors(lattice)
        return _generate_neighbor_bonds(lattice, vectors, boundary, neighbors)
    end
end

# Multiple neighbor orders
function bonds(
    lattice::Lattice{D, Q, T},
    boundary::NTuple{D, Symbol},
    neighbors::Vector{Int}
) where {D, Q, T}
    result = Bond{Q, T, D}[]
    for n in neighbors
        append!(result, bonds(lattice, boundary, n))
    end
    return result
end

# Specific distances
function bonds(
    lattice::Lattice{D, Q, T},
    boundary::NTuple{D, Symbol},
    distances::Vector{<:Real};
    tolerance::Real = 1e-6
) where {D, Q, T}
    vectors = _get_vectors(lattice)
    result = Bond{Q, T, D}[]
    for d in distances
        append!(result, _generate_distance_bonds(lattice, vectors, boundary, d, tolerance))
    end
    return result
end

# Helper: get vectors or error
function _get_vectors(lattice::Lattice)
    if isnothing(lattice.vectors)
        error("lattice.vectors is not set. Use Lattice(unitcell, box_size) with unitcell.vectors set, or set vectors manually.")
    end
    return lattice.vectors
end

# Helper: generate onsite bonds
function _generate_onsite_bonds(lattice::Lattice{D, Q, T}) where {D, Q, T}
    result = Bond{Q, T, D}[]
    zero_sv = zero(SVector{D,T})
    for (state, coord) in zip(lattice.position_states, lattice.coordinates)
        push!(result, Bond{Q,T,D}([state], [coord], [zero_sv]))
    end
    return result
end

# Helper: integer-shift positive direction (lexicographic).
function _is_positive_shift(shift)
    for s in shift
        if s > 0
            return true
        elseif s < 0
            return false
        end
    end
    return false
end

"""
    is_positive_direction(delta::Vector) -> Bool

Check if delta vector is in "positive" direction using lexicographic ordering.
The first non-zero component must be positive.
"""
function is_positive_direction(delta::AbstractVector{T}) where T
    for d in delta
        rd_d = rd(d)
        if rd_d > 0
            return true
        elseif rd_d < 0
            return false
        end
    end
    return true  # zero vector
end

# Helper: build translations for periodic directions
function _translations(boundary::NTuple{D, Symbol}, n::Int) where {D}
    n <= 0 && return Tuple{}[]
    ranges = [boundary[d] == :p ? (-n:n) : (0:0) for d in 1:D]
    shifts = vec(collect(Iterators.product(ranges...)))
    # drop zero shift and keep only one of ±shift
    filter!(s -> any(!=(0), s) && _is_positive_shift(s), shifts)
    return shifts
end

# Helper: convert integer shift to real-space vector
function _shift_vector(
    vectors::Vector{SVector{D,T}},
    shift
) where {D, T}
    v = zero(SVector{D,T})
    for d in 1:D
        v = v + shift[d] * vectors[d]
    end
    return rd(v)
end

# Helper: estimate how many translations are needed for a target distance
function _estimate_shift_range(
    distance::Real,
    vectors::Vector{SVector{D,T}},
    boundary::NTuple{D, Symbol}
) where {D, T}
    lens = [norm(vectors[d]) for d in 1:D if boundary[d] == :p]
    isempty(lens) && return 0
    min_len = minimum(lens)
    return max(1, ceil(Int, distance / min_len) + 1)
end

# Helper: collect distance shells using explicit tiling (QuantumLattices style)
function _distance_levels(
    lattice::Lattice{D, Q, T},
    vectors::Vector{SVector{D,T}},
    boundary::NTuple{D, Symbol},
    max_n::Int,
    tol::Real
) where {D, Q, T}
    coords = lattice.coordinates
    n_sites = length(coords)
    dists = Float64[]

    # Intra-cell distances
    for i in 1:n_sites-1
        for j in (i+1):n_sites
            dist = rd(norm(coords[j] - coords[i]))
            dist > 1e-10 && push!(dists, dist)
        end
    end

    # Inter-cell distances via periodic images
    for shift in _translations(boundary, max_n)
        shift_vec = _shift_vector(vectors, shift)
        for i in 1:n_sites
            c_img = coords[i] + shift_vec
            for j in 1:n_sites
                dist = rd(norm(c_img - coords[j]))
                dist > 1e-10 && push!(dists, dist)
            end
        end
    end

    sort!(dists)
    levels = Float64[]
    for d in dists
        if isempty(levels) || abs(d - levels[end]) > tol
            push!(levels, d)
        end
    end
    return levels
end

# Helper: generate neighbor bonds by order (QuantumLattices style)
function _generate_neighbor_bonds(
    lattice::Lattice{D, Q, T},
    vectors::Vector{SVector{D,T}},
    boundary::NTuple{D, Symbol},
    neighbor_order::Int
) where {D, Q, T}
    tol = 1e-6
    levels = _distance_levels(lattice, vectors, boundary, neighbor_order, tol)
    if neighbor_order > length(levels)
        error("neighbor_order=$neighbor_order exceeds available distance levels ($(length(levels)))")
    end
    return _generate_distance_bonds(lattice, vectors, boundary, levels[neighbor_order], tol)
end

# Helper: generate bonds at specific distance (unidirectional).
#
# Design: explicitly enumerate ALL periodic images of each site (analogous to
# tiling the lattice into a supercluster), then collect every image whose
# distance from site i equals the target distance AND whose displacement is in
# a "positive" direction (to keep exactly one of each i↔j pair).
function _generate_distance_bonds(
    lattice::Lattice{D, Q, T},
    vectors::Vector{SVector{D,T}},
    boundary::NTuple{D, Symbol},
    distance::Real,
    tolerance::Real
) where {D, Q, T}
    result   = Bond{Q, T, D}[]
    coords   = lattice.coordinates
    states   = lattice.position_states
    n_sites  = length(states)
    zero_sv  = zero(SVector{D,T})

    # Intra-cell bonds
    for i in 1:n_sites-1
        for j in (i+1):n_sites
            dist = rd(norm(coords[j] - coords[i]))
            abs(dist - distance) < tolerance || continue
            push!(result, Bond{Q,T,D}(
                [states[i], states[j]],
                [coords[i], coords[j]],
                [zero_sv, zero_sv]
            ))
        end
    end

    # Inter-cell bonds (periodic images)
    nshift = _estimate_shift_range(distance, vectors, boundary)
    for shift in _translations(boundary, nshift)
        shift_vec = _shift_vector(vectors, shift)
        for i in 1:n_sites
            c_img = coords[i] + shift_vec
            for j in 1:n_sites
                dist = rd(norm(c_img - coords[j]))
                abs(dist - distance) < tolerance || continue
                push!(result, Bond{Q,T,D}(
                    [states[j], states[i]],
                    [coords[j], rd(c_img)],
                    [zero_sv, rd(shift_vec)]
                ))
            end
        end
    end

    return result
end
