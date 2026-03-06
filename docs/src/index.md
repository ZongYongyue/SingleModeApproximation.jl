```@raw html
<img src="https://raw.githubusercontent.com/ZongYongyue/MeanFieldTheories.jl/main/docs/src/logo.png" width="300">
```

# MeanFieldTheories.jl

**MeanFieldTheories.jl** is a Julia package for studying quantum many-body systems using mean-field theory and related methods. It provides a complete workflow from constructing many-body Hamiltonians to obtaining self-consistent ground states and calculating collective excitation spectra, covering methods such as Hartree-Fock (HF), Single-Mode Approximation (SMA), and Random Phase Approximation (RPA).

## Overview

The package implements a systematic approach to quantum many-body problems:

1. **Quantum System Construction**: Define degrees of freedom (spin, orbital, sublattice, valley, etc.) on arbitrary lattice geometries, generate one-body and two-body operator terms in creation-annihilation alternating order.
2. **Real-Space Hartree-Fock** (`solve_hf`): Self-consistent field (SCF) iteration in the full $N \times N$ basis. Supports symmetry blocks, DIIS acceleration, multiple random restarts, and finite-temperature Fermi-Dirac occupation.
3. **Momentum-Space Hartree-Fock** (`solve_hfk`): Exploits translational symmetry to reduce the problem to $N_k$ independent $d \times d$ eigenvalue problems. FFT-accelerated self-energy evaluation for density-density (Case A), exchange-type (Case B), and pair-hopping (Case C) interactions.
4. **Single-Mode Approximation**: Calculate momentum-dependent collective excitation spectra on top of the HF ground state.

## Installation

```julia
using Pkg
Pkg.develop(url="https://github.com/ZongYongyue/MeanFieldTheories.jl")
```

## Quick Start

Square-lattice Hubbard model ($t = 1$, $U = 4$, half-filling on a 4×4 lattice) solved by real-space Hartree-Fock. The four cases below demonstrate the trade-off between spin-block structure, memory efficiency, and convergence to the global minimum.

```julia
using MeanFieldTheories

# Square lattice with periodic boundary conditions
unitcell = Lattice([Dof(:site, 1)], [QN(site=1)], [[0.0, 0.0]])
lattice  = Lattice(unitcell, [[1.0, 0.0], [0.0, 1.0]], (4, 4))
t = 1.0;  U = 4.0

# Operators (shared by all four cases below)
function make_ops(dofs)
    t_ops = generate_onebody(dofs, bonds(lattice, (:p, :p), 1),
        (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -t : 0.0).ops
    U_ops = generate_twobody(dofs, bonds(lattice, (:p, :p), 0),
        (deltas, qn1, qn2, qn3, qn4) ->
            (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1, 1, 2, 2) ? U : 0.0,
        order = (cdag, 1, c, 1, cdag, 1, c, 1)).ops
    vcat(t_ops, U_ops)
end
```

**Case 1 — spin blocks, single run:**

```julia
dofs_b = SystemDofs([Dof(:site, 16), Dof(:spin, 2, [:up, :down])], sortrule = [[2], 1])
result = solve_hf(dofs_b, make_ops(dofs_b), [8, 8])
```

```
============================================================
Hartree-Fock SCF Solver
============================================================
Building Hamiltonian  (144 operators)
           t matrix: (32, 32), nnz = 128      243.021ms
           U matrix: (1024, 1024), nnz = 32        28.811ms
System: N = 32, blocks = 2, particles = [8, 8] (total = 16)
T = 0,  mixing = DIIS(m=8),  tol = 1e-08,  max_iter = 1000
============================================================
G initialized   139.263ms
Iter    1  res = 7.805e-03  E = -24.155896  NCond = 1.4441
...
Iter   53  res = 9.972e-09 < 1.000e-08  CONVERGED
============================================================
SCF CONVERGED  (53 iterations)
  Band energy:        +2.3594869312
  Interaction energy: -12.4724545253
  Total energy:       -10.1129675941
  NCond:              16.000000
  Sz:                 -0.000000
  μ (block 1):       +2.0000000007
  μ (block 2):       +1.9999999993
```

**Case 2 — spin blocks, 10 random restarts:**

```julia
result = solve_hf(dofs_b, make_ops(dofs_b), [8, 8], n_restarts = 10)
```

```
============================================================
  Restart  1: E = -10.7471610309  (CONVERGED, 20 iters)
  Restart  2: E = -10.1129683059  (CONVERGED, 44 iters)
  ...
  Restart  7: E = -12.5665533115  (CONVERGED, 16 iters)
  ...
============================================================
Best result from 10 restarts:
SCF CONVERGED  (16 iterations)
  Total energy:       -12.5665533115
  NCond:              16.000000
  Sz:                 +0.000000
```

**Case 3 — no block structure, single run:**

```julia
dofs_nb = SystemDofs([Dof(:site, 16), Dof(:spin, 2, [:up, :down])])
result = solve_hf(dofs_nb, make_ops(dofs_nb), [16])
```

```
SCF CONVERGED  (19 iterations)
  Total energy:       -12.5665518227
  NCond:              16.000000
  Sz:                 -0.000005
```

**Case 4 — no block structure, 10 random restarts:**

```julia
result = solve_hf(dofs_nb, make_ops(dofs_nb), [16], n_restarts = 10)
```

```
  Restart  2: E = -12.5665594104  (CONVERGED, 17 iters)
  ...
Best result from 10 restarts:
SCF CONVERGED  (22 iterations)
  Total energy:       -12.5665661434
  NCond:              16.000000
```

> **Block structure trade-off.** Declaring spin as a symmetry block (Cases 1–2 vs 3–4) halves the number of non-zeros in the interaction matrix (nnz 64 → 32), reducing memory and speeding up each SCF iteration. However, enforcing spin conservation at the density-matrix level means the Green's function is constrained to be strictly block-diagonal throughout the iteration — off-diagonal spin correlations are permanently forbidden. As a result, the blocked solver cannot use spin-mixed intermediate configurations to escape local minima. The unblocked solver, by contrast, has full freedom in the density matrix during the SCF and can reach the global minimum via spin-mixed intermediates. In 10 runs the unblocked solver (Case 4) finds the global minimum ($E \approx -12.57$) 7 out of 10 times, while the blocked solver (Case 2) finds it only 1 out of 10 times. The practical recommendation: use block structure for efficiency, but always pair it with `n_restarts` to compensate for the reduced basin of attraction.

## Documentation Contents

```@contents
Pages = ["quantumsystem.md", "hartreefock_real.md", "hartreefock_momentum.md", "singlemode.md"]
Depth = 2
```

## Package Structure

- **quantumsystem**: Core quantum system definitions, lattice structures, and operator algebra
- **groundstate**: Ground state calculations — real-space HF (`solve_hf`) and momentum-space HF (`solve_hfk`)
- **excitations**: Excited state calculations — Single-Mode Approximation and Random Phase Approximation

## Citation

If you use this package in your research, please cite:

```bibtex
@software{meanfieldtheories,
  author = {Yong-Yue Zong},
  title  = {MeanFieldTheories.jl: A Julia package for quantum many-body systems using mean-field theory and related methods},
  year   = {2026},
  url    = {https://github.com/ZongYongyue/MeanFieldTheories.jl}
}
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests on [GitHub](https://github.com/ZongYongyue/MeanFieldTheories.jl).

## License

This project is licensed under the MIT License — see the LICENSE file for details.
