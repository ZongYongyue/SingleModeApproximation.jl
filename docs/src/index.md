```@raw html
<img src="https://raw.githubusercontent.com/Quantum-Many-Body/MeanFieldTheories.jl/main/docs/src/logo.png" width="300">
```

# MeanFieldTheories.jl

**MeanFieldTheories.jl** is a Julia package for studying quantum many-body systems using mean-field theory and related methods. It provides a complete workflow from constructing many-body Hamiltonians to obtaining self-consistent ground states and calculating collective excitation spectra, covering methods such as Hartree-Fock (HF), Single-Mode Approximation (SMA), and Random Phase Approximation (RPA).

## Overview

The package implements a systematic approach to quantum many-body problems:

1. **Quantum System Construction**: Define degrees of freedom (spin, orbital, sublattice, valley, etc.) on arbitrary lattice geometries, generate one-body and two-body operator terms in creation-annihilation alternating order.
2. **Real-Space Hartree-Fock** (`solve_hfr`): Self-consistent field (SCF) iteration in the full $N \times N$ basis. Supports symmetry blocks, DIIS acceleration, multiple random restarts, and finite-temperature Fermi-Dirac occupation.
3. **Momentum-Space Hartree-Fock** (`solve_hfk`): Exploits translational symmetry to reduce the problem to $N_k$ independent $d \times d$ eigenvalue problems. FFT-accelerated self-energy evaluation for density-density (Case A), exchange-type (Case B), and pair-hopping (Case C) interactions.
4. **Single-Mode Approximation**: Calculate momentum-dependent collective excitation spectra on top of the HF ground state.

## Installation

```julia
using Pkg
Pkg.develop(url="https://github.com/Quantum-Many-Body/MeanFieldTheories.jl")
```

## Quick Start

### Real-Space Hartree-Fock (`solve_hfr`)

Hubbard model ($t = 1$, $U = 8$) on an 8-site $\sqrt{8}\times\sqrt{8}$ square-lattice cluster (45°-rotated supercell with PBC), half-filling ($N_e = 8$, $S_z = 0$):

$$H = -t \sum_{\langle ij \rangle, \sigma} \left(c^\dagger_{i\sigma} c_{j\sigma} + \mathrm{h.c.}\right) + U \sum_i n_{i\uparrow} n_{i\downarrow}$$

```julia
using MeanFieldTheories

t = 1.0;  U = 8.0

dofs = SystemDofs([Dof(:site, 8), Dof(:spin, 2)], sortrule = [[2], 1])

# 8 sites in the √8×√8 fundamental domain; supercell_vectors encodes the tilted PBC
cluster = Lattice([Dof(:site, 8)],
                  [QN(site=i) for i in 1:8],
                  [[0.0,0.0], [-1.0,1.0], [0.0,1.0], [1.0,1.0],
                   [-1.0,2.0], [0.0,2.0], [1.0,2.0], [0.0,3.0]];
                  supercell_vectors=[[2.0,2.0],[-2.0,2.0]])

t_ops = generate_onebody(dofs, bonds(cluster, (:p, :p), 1),
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -t : 0.0,
    hc = true).ops

U_ops = generate_twobody(dofs, bonds(cluster, (:p, :p), 0),
    (deltas, qn1, qn2, qn3, qn4) ->
        (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1, 1, 2, 2) ? U : 0.0,
    order = (cdag, :i, c, :i, cdag, :i, c, :i)).ops

result = solve_hfr(dofs, vcat(t_ops, U_ops), [4, 4];
    n_restarts = 20, tol = 1e-8, mix_alpha = 0.5)
```

```
============================================================
Best result from 20 restarts:
SCF CONVERGED  (15 iterations)
  Band energy:        -0.9265415224
  Interaction energy: -2.4816524388
  Total energy:       -3.4081939612
  NCond:              8.000000
  Sz:                 -0.000000
```

### Momentum-Space Hartree-Fock (`solve_hfk`)

V-model ($t = 1$, $V = 4$) on the 2D square lattice at half-filling, with a $2\times 2$ magnetic unit cell that allows $Q = (\pi,\pi)$ antiferromagnetic order:

$$H = -t \sum_{\langle ij \rangle, \sigma} \left(c^\dagger_{i\sigma} c_{j\sigma} + \mathrm{h.c.}\right) + V \sum_{i \neq j, \sigma\sigma'} n_{i\sigma} n_{j\sigma'}$$

```julia
using MeanFieldTheories

t = 1.0;  V = 4.0

# 2×2 magnetic unit cell: 4 sites × 2 spins → d = 8 per k-point
dofs = SystemDofs([Dof(:site, 4), Dof(:spin, 2, [:up, :dn])])

unitcell = Lattice([Dof(:site, 4)],
                   [QN(site=i) for i in 1:4],
                   [[0.0,0.0], [1.0,0.0], [0.0,1.0], [1.0,1.0]];
                   supercell_vectors=[[2.0,0.0],[0.0,2.0]])

nn_bonds = bonds(unitcell, (:p, :p), 1)  # intra- and inter-cell NN bonds

onebody = generate_onebody(dofs, nn_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -t : 0.0)

# :i and :j are independent position indices → full ordered-pair sum Σ_{i≠j}
twobody = generate_twobody(dofs, nn_bonds,
    (deltas, qn1, qn2, qn3, qn4) ->
        qn1.spin == qn2.spin && qn3.spin == qn4.spin ? V : 0.0)

# 2×2 k-grid in the magnetic BZ; total 4×4 = 16 sites, 16 electrons
kpoints = build_kpoints([[2.0,0.0],[0.0,2.0]], (2, 2))

result = solve_hfk(dofs, onebody, twobody, kpoints, 16; n_restarts = 20, tol = 1e-8)
```

```
============================================================
Best result from 20 restarts:
SCF CONVERGED  (9 iterations)
  Band energy:        -0.0095509951
  Interaction energy: -0.5570138562
  Total energy:       -0.5665648513   ← per magnetic unit cell
  NCond:              4.000000
```

## Documentation Contents

```@contents
Pages = ["quantumsystem.md", "hartreefock_real.md", "hartreefock_momentum.md", "singlemode.md", "api.md"]
Depth = 2
```

## Package Structure

- **quantumsystem**: Core quantum system definitions, lattice structures, and operator algebra
- **groundstate**: Ground state calculations — real-space HF (`solve_hfr`) and momentum-space HF (`solve_hfk`)
- **excitations**: Excited state calculations — Single-Mode Approximation and Random Phase Approximation

## Citation

If you use this package in your research, please cite:

```bibtex
@software{meanfieldtheories,
  author = {Yong-Yue Zong},
  title  = {MeanFieldTheories.jl: A Julia package for quantum many-body systems using mean-field theory and related methods},
  year   = {2026},
  url    = {https://github.com/Quantum-Many-Body/MeanFieldTheories.jl}
}
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests on [GitHub](https://github.com/Quantum-Many-Body/MeanFieldTheories.jl).

## License

This project is licensed under the MIT License — see the LICENSE file for details.
