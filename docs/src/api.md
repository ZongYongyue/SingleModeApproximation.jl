# API Reference

## Quantum System

### Types

```@docs
QuantumNumber
Dof
SystemDofs
Lattice
Bond
FermionOp
Operators
```

### Quantum Number Functions

```@docs
qn2linear
linear2qn
ndofs
total_dim
site_spin_system
```

### Lattice Functions

```@docs
get_coordinate
bonds
is_positive_direction
```

### Operator Constructors

```@docs
c
cdag
```

### Operator Generators

```@docs
generate_onebody
generate_twobody
build_onebody_matrix
build_interaction_tensor
```

---

## Real-Space Hartree-Fock

```@docs
build_T
build_U
solve_hf
```

---

## Momentum-Space Hartree-Fock

### Preprocessing

```@docs
build_Tr
build_Tk
build_Vr
build_Vk
build_Uk
```

### k-point Utilities

```@docs
build_kpoints
initialize_green_k
green_k_to_tau
```

### SCF Solver

```@docs
solve_hfk
```

---

## Utilities

```@docs
rd
PRECISION
```
