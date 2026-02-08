"""
    SingleModeApproximation

A Julia package for quantum many-body systems using single-mode approximation methods.

This package provides tools for:
- Defining quantum systems with arbitrary degrees of freedom
- Constructing lattice structures
- Generating operator terms with flexible orderings
- Building matrix and tensor representations
- Hartree-Fock calculations (ground state and excited states)

# Modules
- `quantumsystem`: Core quantum system definitions and operator algebra
- `groundstate`: Ground state calculations (Hartree-Fock, mean-field methods)
- `excitations`: Excited state calculations (RPA, TDHF, etc.)

# Example
```julia
using SingleModeApproximation

# Define system
dofs = SystemDofs(site=1:4, spin=[:up, :dn])
lattice = Lattice(:Square, 2, 2, pbc=true)

# Generate hopping terms
bonds = bonds(lattice, 1)
hopping = generate_onebody(dofs, bonds, -1.0)

# Build Hamiltonian
H = build_onebody_matrix(dofs, hopping)
```
"""
module SingleModeApproximation

# Export quantum system functionality
include("quantumsystem/freedom.jl")
include("quantumsystem/lattice.jl")
include("quantumsystem/operators.jl")

# Export core types
export QuantumNumber, QN, Dof, SystemDofs
export Lattice, Bond
export Operator, FermionOp, Operators
export c, cdag

# Export lattice functions
export bonds, get_coordinate, is_positive_direction

# Export quantum number functions
export qn2linear, linear2qn
export ndofs, total_dim, site_spin_system

# Export operator generation functions
export generate_onebody, generate_twobody
export generate_coulomb_intra, generate_coulomb_inter
export generate_hund, generate_ising, generate_exchange, generate_pair_hop

# Export matrix/tensor builders
export build_onebody_matrix, build_interaction_tensor

# Export utility constants
export PRECISION, rd

end # module SingleModeApproximation
