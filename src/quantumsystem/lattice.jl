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
- `coordinates::Vector{Vector{T}}`: Coordinate for each position state (length D)
- `supercell_vectors::Union{Nothing, Vector{Vector{T}}}`: Supercell vectors for PBC (optional)

# Notes
- `position_states[i]` corresponds to `coordinates[i]`
- `supercell_vectors` is set automatically when using tiling constructor

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
    coordinates::Vector{Vector{T}}
    supercell_vectors::Union{Nothing, Vector{Vector{T}}}

    # Inner constructor that infers type parameters
    function Lattice(
        position_dofs::Vector{Dof},
        position_states::Vector{Q},
        coordinates::Vector{Vector{T}};
        supercell_vectors::Union{Nothing, Vector{Vector{T}}} = nothing
    ) where {Q<:QuantumNumber, T<:Real}
        D = length(coordinates[1])
        new{D, Q, T}(position_dofs, position_states, coordinates, supercell_vectors)
    end
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
    if isnothing(lattice.supercell_vectors)
        print(io, "Supercell vectors: not set")
    else
        print(io, "Supercell vectors: set")
    end
end

"""
    Lattice(unitcell, unitcell_vectors, box_size)

Create a lattice by tiling a unit cell.

The first DOF of unitcell defines the unit cell type (with size=1 in unitcell).
After tiling, its size becomes Nx*Ny*... (total number of unit cells).
The supercell_vectors are automatically computed as box_size .* unitcell_vectors.

# Arguments
- `unitcell::Lattice`: The unit cell definition (first DOF has size=1)
- `unitcell_vectors::Vector{Vector{T}}`: Unit cell translation vectors [a1, a2, ...]
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

# Tile with unit cell vectors
a1, a2 = [1.0, 0.0], [0.5, 0.866]
lattice = Lattice(unitcell, [a1, a2], (4, 4))
# Result: DOFs [honeycomb_cell (16), sublattice (2)], 32 sites total
# supercell_vectors = [4*a1, 4*a2] automatically set
```
"""
function Lattice(
    unitcell::Lattice{D, Q, T},
    unitcell_vectors::Vector{Vector{T}},
    box_size::NTuple{M, Int}
) where {D, Q, T, M}
    @assert M == D "box_size dimension ($M) must match lattice dimension ($D)"
    @assert length(unitcell_vectors) == D "Need $D unit cell vectors for $D-dimensional lattice"
    @assert unitcell.position_dofs[1].size == 1 "First DOF of unitcell must have size=1"

    # Number of unit cells
    n_cells = prod(box_size)

    # Compute supercell vectors
    supercell_vectors = [rd(box_size[i] .* unitcell_vectors[i]) for i in 1:D]

    # Create new DOFs: expand first DOF's size, keep others
    first_dof = unitcell.position_dofs[1]
    new_dofs = [Dof(first_dof.name, n_cells, first_dof.labels); unitcell.position_dofs[2:end]]

    # Get names from unitcell (same structure, just different values for first DOF)
    all_names = Tuple(dof.name for dof in unitcell.position_dofs)

    # Generate position_states and coordinates
    N = QuantumNumber{all_names, NTuple{length(all_names), Int}}
    position_states = N[]
    coordinates = Vector{T}[]

    # Iterate over all unit cells (linear index)
    cell_idx = 1
    for indices in Iterators.product([1:s for s in box_size]...)
        # Calculate unit cell origin: R0 = (i1-1)*a1 + (i2-1)*a2 + ...
        R0 = zeros(T, D)
        for (dim, i) in enumerate(indices)
            R0 .+= (i - 1) .* unitcell_vectors[dim]
        end
        R0 = rd(R0)

        # Add each site in the unit cell
        for (cell_state, cell_coord) in zip(unitcell.position_states, unitcell.coordinates)
            # Replace first DOF value with cell_idx, keep rest
            state_values = values(cell_state)
            new_nt = NamedTuple{all_names}((cell_idx, state_values[2:end]...))
            push!(position_states, QuantumNumber(new_nt))

            # Calculate coordinate
            new_coord = rd(R0 .+ cell_coord)
            push!(coordinates, new_coord)
        end

        cell_idx += 1
    end

    return Lattice(new_dofs, position_states, coordinates; supercell_vectors=supercell_vectors)
end

#==================== Bond Specification ====================#

"""
    Bond{Q, T}

Represents a bond connecting one or more sites.

# Fields
- `states::Vector{Q}`: Position states (length determines bond order: 1=onsite, 2=two-body, etc.)
- `coordinates::Vector{Vector{T}}`: Coordinates for each site

# Examples
```julia
# Onsite bond (1 site)
bond = Bond([state1], [coord1])

# Two-body bond (2 sites)
bond = Bond([state1, state2], [coord1, coord2])
```
"""
struct Bond{Q<:QuantumNumber, T<:Real}
    states::Vector{Q}
    coordinates::Vector{Vector{T}}
end

function Base.show(io::IO, bond::Bond)
    print(io, "Bond(", bond.states, ", ", bond.coordinates, ")")
end

"""
    bonds(lattice, boundary, neighbors)

Generate bonds based on neighbor specification.

# Arguments
- `lattice::Lattice`: The lattice structure (must have supercell_vectors set)
- `boundary::NTuple{D, Symbol}`: Boundary conditions, :p (periodic) or :o (open)
- `neighbors`: Neighbor specification:
  - `Int`: n-th nearest neighbor (0=onsite, 1=nearest, 2=next-nearest, ...)
  - `Vector{Int}`: multiple neighbor orders [n1, n2, ...]
  - `Vector{<:Real}`: specific distances [r1, r2, ...]

# Returns
`Vector{Bond}` with appropriate number of sites per bond

# Examples
```julia
# Create lattice with tiling (supercell_vectors auto-set)
lattice = Lattice(unitcell, [a1, a2], (4, 4))

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
        supercell_vectors = _get_supercell_vectors(lattice)
        return _generate_neighbor_bonds(lattice, supercell_vectors, boundary, neighbors)
    end
end

# Multiple neighbor orders
function bonds(
    lattice::Lattice{D, Q, T},
    boundary::NTuple{D, Symbol},
    neighbors::Vector{Int}
) where {D, Q, T}
    result = Bond{Q, T}[]
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
    supercell_vectors = _get_supercell_vectors(lattice)
    result = Bond{Q, T}[]
    for d in distances
        append!(result, _generate_distance_bonds(lattice, supercell_vectors, boundary, d, tolerance))
    end
    return result
end

# Helper: get supercell_vectors or error
function _get_supercell_vectors(lattice::Lattice)
    if isnothing(lattice.supercell_vectors)
        error("lattice.supercell_vectors is not set. Use Lattice(unitcell, vectors, box_size) or set supercell_vectors manually.")
    end
    return lattice.supercell_vectors
end

# Helper: generate onsite bonds
function _generate_onsite_bonds(lattice::Lattice{D, Q, T}) where {D, Q, T}
    result = Bond{Q, T}[]
    for (state, coord) in zip(lattice.position_states, lattice.coordinates)
        push!(result, Bond([state], [coord]))
    end
    return result
end

# Helper: compute minimum image distance considering PBC
function _min_image_distance(
    coord1::Vector{T},
    coord2::Vector{T},
    supercell_vectors::Vector{Vector{T}},
    boundary::NTuple{D, Symbol}
) where {D, T}
    # Try all periodic images and find minimum distance
    min_dist = Inf
    delta = coord2 .- coord1

    # For each periodic direction, try shifts -1, 0, +1
    for shifts in Iterators.product([boundary[i] == :p ? (-1, 0, 1) : (0,) for i in 1:D]...)
        shifted_delta = copy(delta)
        for (dim, shift) in enumerate(shifts)
            shifted_delta .+= shift .* supercell_vectors[dim]
        end
        dist = sqrt(sum(shifted_delta .^ 2))
        if dist < min_dist
            min_dist = dist
        end
    end

    return rd(min_dist)
end

# Helper: get actual delta vector with minimum image convention
function _min_image_delta(
    coord1::Vector{T},
    coord2::Vector{T},
    supercell_vectors::Vector{Vector{T}},
    boundary::NTuple{D, Symbol}
) where {D, T}
    min_dist = Inf
    best_delta = coord2 .- coord1

    for shifts in Iterators.product([boundary[i] == :p ? (-1, 0, 1) : (0,) for i in 1:D]...)
        delta = coord2 .- coord1
        for (dim, shift) in enumerate(shifts)
            delta .+= shift .* supercell_vectors[dim]
        end
        dist = sqrt(sum(delta .^ 2))
        if dist < min_dist
            min_dist = dist
            best_delta = delta
        end
    end

    return rd(best_delta), rd(min_dist)
end

# Helper: generate neighbor bonds by order
function _generate_neighbor_bonds(
    lattice::Lattice{D, Q, T},
    supercell_vectors::Vector{Vector{T}},
    boundary::NTuple{D, Symbol},
    neighbor_order::Int
) where {D, Q, T}
    n_sites = length(lattice.position_states)

    # Compute all unique distances
    distances = Set{T}()
    for i in 1:n_sites
        for j in (i+1):n_sites
            d = _min_image_distance(
                lattice.coordinates[i],
                lattice.coordinates[j],
                supercell_vectors, boundary
            )
            if d > 1e-10  # exclude self
                push!(distances, d)  # already rounded by _min_image_distance
            end
        end
    end

    # Sort distances and get the n-th smallest
    sorted_distances = sort(collect(distances))
    if neighbor_order > length(sorted_distances)
        error("neighbor_order=$neighbor_order exceeds available distance levels ($(length(sorted_distances)))")
    end

    target_distance = sorted_distances[neighbor_order]
    tolerance = 1e-6

    return _generate_distance_bonds(lattice, supercell_vectors, boundary, target_distance, tolerance)
end

"""
    is_positive_direction(delta::Vector) -> Bool

Check if delta vector is in "positive" direction using lexicographic ordering.
The first non-zero component must be positive.

Used to generate unidirectional bonds, avoiding duplicates like (i→j) and (j→i).
"""
function is_positive_direction(delta::Vector{T}) where T
    for d in delta
        rd_d = rd(d)
        if rd_d > 0
            return true
        elseif rd_d < 0
            return false
        end
    end
    return true  # zero vector (shouldn't happen for valid bonds)
end

# Helper: generate bonds at specific distance (unidirectional)
function _generate_distance_bonds(
    lattice::Lattice{D, Q, T},
    supercell_vectors::Vector{Vector{T}},
    boundary::NTuple{D, Symbol},
    distance::Real,
    tolerance::Real
) where {D, Q, T}
    result = Bond{Q, T}[]
    n_sites = length(lattice.position_states)

    for i in 1:n_sites
        for j in 1:n_sites
            if i == j
                continue
            end

            delta, d = _min_image_delta(
                lattice.coordinates[i],
                lattice.coordinates[j],
                supercell_vectors, boundary
            )

            if abs(d - distance) < tolerance
                # Only keep positive direction bonds (unidirectional)
                if !is_positive_direction(delta)
                    continue
                end

                # Create bond with actual coordinates (not wrapped)
                coord1 = lattice.coordinates[i]
                coord2 = rd(coord1 .+ delta)

                bond = Bond(
                    [lattice.position_states[i], lattice.position_states[j]],
                    [coord1, coord2]
                )
                push!(result, bond)
            end
        end
    end

    return result
end
