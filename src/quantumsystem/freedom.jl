"""
Degree of Freedom Definitions for Quantum Systems

Provides types and functions for defining arbitrary quantum numbers
(site, spin, orbital, valley, layer, etc.) and their valid combinations.
"""

using LinearAlgebra

#==================== Global Precision Setting ====================#

"""
Global precision for numerical rounding (number of decimal places).
All floating-point values (coordinates, hopping amplitudes, interactions, phases, etc.)
are rounded to this precision to avoid floating-point comparison issues.
"""
const PRECISION = 6

"""
    rd(x) -> Number

Round a number to PRECISION decimal places. Works with Real and Complex.
"""
rd(x::Real) = round(x, digits=PRECISION)
rd(x::Complex) = Complex(rd(real(x)), rd(imag(x)))
rd(v::AbstractVector) = rd.(v)

#==================== Quantum Number Type ====================#

"""
    QuantumNumber{Names, Types}

A wrapper around NamedTuple for representing quantum numbers with cleaner display.

# Type Parameters
- `Names`: Tuple of field names (Symbols)
- `Types`: Tuple type for values

# Examples
```julia
# Create from keyword arguments
qn = QuantumNumber(site=1, spin=2)
qn.site  # 1
qn.spin  # 2

# Display: QN(site=1, spin=2)
```
"""
struct QuantumNumber{Names, Types<:Tuple}
    data::NamedTuple{Names, Types}

    # Inner constructor from NamedTuple
    function QuantumNumber{Names, Types}(nt::NamedTuple{Names, Types}) where {Names, Types<:Tuple}
        new{Names, Types}(nt)
    end
end

# Outer constructor from keyword arguments
QuantumNumber(; kwargs...) = QuantumNumber(NamedTuple(kwargs))

# Outer constructor from NamedTuple (generic)
QuantumNumber(nt::NamedTuple{Names, Types}) where {Names, Types<:Tuple} =
    QuantumNumber{Names, Types}(nt)

# Property access: qn.field -> qn.data.field
function Base.getproperty(qn::QuantumNumber, name::Symbol)
    if name === :data
        return getfield(qn, :data)
    else
        return getproperty(getfield(qn, :data), name)
    end
end

# Support tab completion
Base.propertynames(qn::QuantumNumber) = propertynames(getfield(qn, :data))

# Equality
Base.:(==)(a::QuantumNumber, b::QuantumNumber) = getfield(a, :data) == getfield(b, :data)
Base.hash(qn::QuantumNumber, h::UInt) = hash(getfield(qn, :data), h)

# Iteration and indexing (for compatibility)
Base.keys(qn::QuantumNumber) = keys(getfield(qn, :data))
Base.values(qn::QuantumNumber) = values(getfield(qn, :data))
Base.getindex(qn::QuantumNumber, key::Symbol) = getfield(qn, :data)[key]
Base.length(qn::QuantumNumber) = length(getfield(qn, :data))

# Display: QN(field1=val1, field2=val2)
function Base.show(io::IO, qn::QuantumNumber)
    print(io, "QN(")
    data = getfield(qn, :data)
    ks = keys(data)
    for (i, k) in enumerate(ks)
        i > 1 && print(io, ", ")
        print(io, k, "=", data[k])
    end
    print(io, ")")
end

# Convert to NamedTuple
Base.convert(::Type{NamedTuple}, qn::QuantumNumber) = getfield(qn, :data)
Base.NamedTuple(qn::QuantumNumber) = getfield(qn, :data)

# Type alias for convenience
const QN = QuantumNumber

#==================== Quantum Number Specification ====================#

"""
    Dof

Specification for a single degree of freedom (DOF).

# Fields
- `name::Symbol`: Name of the DOF (e.g., :site, :spin, :orbital)
- `size::Int`: Number of states for this DOF
- `labels::Union{Nothing, Vector}`: Optional labels for states (e.g., [:up, :down])
"""
struct Dof
    name::Symbol
    size::Int
    labels::Union{Nothing, Vector}

    function Dof(name::Symbol, size::Int, labels=nothing)
        if labels !== nothing && length(labels) != size
            error("Number of labels must match size")
        end
        new(name, size, labels)
    end
end

"""
    SystemDofs{Q}

Defines all degrees of freedom in the system, with optional constraints and custom ordering.

# Type Parameter
- `Q`: QuantumNumber type with DOF names as keys

# Fields
- `dofs::Vector{Dof}`: List of degrees of freedom
- `valid_states::Vector{Q}`: All valid quantum number combinations (sorted)
- `qn_to_idx::Dict{Q, Int}`: Fast O(1) lookup from quantum numbers to linear indices

# Constructor
```julia
SystemDofs(dofs::Vector{Dof}; constraint=nothing, sortrule=collect(length(dofs):-1:1))
```

# Arguments
- `dofs`: Vector of Dof specifications
- `constraint`: Optional function `qn -> Bool` to filter valid states (not saved after construction)
- `sortrule`: Vector of integers specifying sort priority (not saved after construction)
  - Default: `[N, N-1, ..., 1]` (reverse order: first dof varies fastest)
  - Elements are indices into `dofs` array
  - Sorts by `dofs[sortrule[1]]` first, then `dofs[sortrule[2]]`, etc.

# Design Rationale
The default reverse order `[N, N-1, ..., 1]` ensures that the first DOF (typically site)
varies fastest, making adjacent sites have adjacent linear indices, which is cache-friendly.

# Examples
```julia
# Example 1: Default ordering (reverse)
dofs = SystemDofs([
    Dof(:site, 3),
    Dof(:spin, 2)
])
# sortrule = [2, 1] by default: spin varies slowest, site varies fastest
# valid_states = [QN(site=1,spin=1), QN(site=2,spin=1), QN(site=3,spin=1),
#                 QN(site=1,spin=2), QN(site=2,spin=2), QN(site=3,spin=2)]

# Example 2: With constraint (spin-valley locking)
dofs = SystemDofs([
    Dof(:site, 4),
    Dof(:spin, 2),
    Dof(:valley, 2)
], constraint = qn -> qn.spin == qn.valley)

# Example 3: Custom sort order
dofs = SystemDofs([
    Dof(:moire_unitcell, 4),   # index 1
    Dof(:sublattice, 2),        # index 2
    Dof(:spin, 2),              # index 3
    Dof(:valley, 2)             # index 4
], sortrule = [4, 3, 2, 1])
# Sorts by: valley → spin → sublattice → moire_unitcell (slow to fast)

# Example 4: Forward ordering (first dof slowest)
dofs = SystemDofs([
    Dof(:site, 3),
    Dof(:spin, 2)
], sortrule = [1, 2])
# valid_states = [QN(site=1,spin=1), QN(site=1,spin=2),
#                 QN(site=2,spin=1), QN(site=2,spin=2), ...]
```
"""
struct SystemDofs{Q<:QuantumNumber}
    dofs::Vector{Dof}
    valid_states::Vector{Q}
    qn_to_idx::Dict{Q, Int}
end

# Constructor with keyword arguments
function SystemDofs(
    dofs::Vector{Dof};
    constraint::Union{Nothing, Function} = nothing,
    sortrule::Vector{Int} = collect(length(dofs):-1:1)
)
    # Validate sortrule
    @assert length(sortrule) == length(dofs) "sortrule must have same length as dofs"
    @assert sort(sortrule) == collect(1:length(dofs)) "sortrule must be a permutation of 1:$(length(dofs))"

    # Generate all valid states
    names = Tuple(dof.name for dof in dofs)
    Q = QuantumNumber{names, NTuple{length(dofs), Int}}
    valid_states = Q[]
    _enumerate_states!(valid_states, dofs, names, Int[], 1, constraint)

    # Sort according to sortrule
    sort!(valid_states, by = qn -> tuple([qn[dofs[i].name] for i in sortrule]...))

    # Build index lookup dictionary
    qn_to_idx = Dict(qn => i for (i, qn) in enumerate(valid_states))

    SystemDofs{Q}(dofs, valid_states, qn_to_idx)
end

# Helper: recursively enumerate states (unified for with/without constraint)
function _enumerate_states!(result::Vector{Q}, dofs::Vector{Dof}, names::Tuple,
                            current::Vector{Int}, depth::Int,
                            constraint::Union{Nothing, Function}) where {Q<:QuantumNumber}
    if depth > length(dofs)
        nt = NamedTuple{names}(Tuple(current))
        qn = QuantumNumber(nt)
        if constraint === nothing || constraint(qn)
            push!(result, qn)
        end
        return
    end

    for i in 1:dofs[depth].size
        push!(current, i)
        _enumerate_states!(result, dofs, names, current, depth + 1, constraint)
        pop!(current)
    end
end

# Accessor functions
"""
    ndofs(dofs::SystemDofs) -> Int

Return number of degrees of freedom.
"""
ndofs(d::SystemDofs) = length(d.dofs)

"""
    total_dim(dofs::SystemDofs) -> Int

Return total dimension (number of valid states).
"""
total_dim(d::SystemDofs) = length(d.valid_states)

# Convenience constructors
"""
    site_spin_system(Nsite::Int) -> SystemDofs

Create a simple site+spin system (standard case).
"""
function site_spin_system(Nsite::Int)
    return SystemDofs([
        Dof(:site, Nsite),
        Dof(:spin, 2, [:up, :down])
    ])
end

"""
    site_orbital_spin_system(Nsite::Int, Norbital::Int) -> SystemDofs

Create a site+orbital+spin system (multi-orbital case).
"""
function site_orbital_spin_system(Nsite::Int, Norbital::Int)
    return SystemDofs([
        Dof(:site, Nsite),
        Dof(:orbital, Norbital),
        Dof(:spin, 2, [:up, :down])
    ])
end

#==================== SystemDofs Display ====================#

"""
    Base.show(io::IO, dofs::SystemDofs)

Custom display for SystemDofs in REPL.
"""
function Base.show(io::IO, dofs::SystemDofs)
    n = length(dofs.dofs)
    dim = length(dofs.valid_states)

    println(io, "SystemDofs with $n degree(s) of freedom:")
    for (i, dof) in enumerate(dofs.dofs)
        print(io, "  $i. $(dof.name): $(dof.size) states")
        if dof.labels !== nothing
            print(io, " [", join(dof.labels, ", "), "]")
        end
        println(io)
    end
    print(io, "Total dimension: $dim")
end

#==================== Quantum Number Index Functions ====================#

"""
    qn2linear(dofs::SystemDofs, qn) -> Int

Convert quantum numbers to linear index (position in valid_states).
Uses O(1) hash table lookup for fast performance.
Accepts QuantumNumber or NamedTuple.

# Examples
```julia
dofs = SystemDofs([Dof(:site, 4), Dof(:spin, 2)])
idx = qn2linear(dofs, QN(site=1, spin=2))  # O(1) lookup
idx = qn2linear(dofs, (site=1, spin=2))    # also works
```
"""
function qn2linear(dofs::SystemDofs, qn::QuantumNumber)
    idx = get(dofs.qn_to_idx, qn, nothing)
    if isnothing(idx)
        error("Quantum numbers $qn not in valid states")
    end
    return idx
end

# Also accept NamedTuple for convenience
function qn2linear(dofs::SystemDofs, qn::NamedTuple)
    return qn2linear(dofs, QuantumNumber(qn))
end

"""
    linear2qn(dofs::SystemDofs, idx::Int) -> QuantumNumber

Convert linear index back to quantum numbers.

# Examples
```julia
dofs = site_spin_system(4)
qn = linear2qn(dofs, 2)  # returns QN(site=1, spin=2)
```
"""
function linear2qn(dofs::SystemDofs, idx::Int)
    return dofs.valid_states[idx]
end
