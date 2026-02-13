# Hartree-Fock Approximation

This section introduces the Hartree-Fock (HF) approximation method for obtaining self-consistent mean-field solutions of quantum many-body systems.

## Physical Hamiltonian

Consider a general quantum many-body Hamiltonian consisting of one-body and two-body interaction terms:

$$H = H_0 + H_{\text{int}}$$

where the one-body term (kinetic energy, external potential, etc.) is:

$$H_0 = \sum_{ij} t_{ij} \, c^\dagger_i c_j$$

and the two-body interaction term is:

$$H_{\text{int}} = \sum_{ijkl} V_{ijkl} \, c^\dagger_i c_j c^\dagger_k c_l$$

Here:
- $c^\dagger_i, c_j$ are fermionic creation/annihilation operators
- $i = (r_i, \alpha_i, \sigma_i)$ represents a complete set of quantum numbers (position, orbital, spin, etc.)
- $t_{ij}$ is the one-body matrix element (hopping amplitude)
- $V_{ijkl}$ is the two-body interaction matrix element in **InterAll format** ($c^\dagger c c^\dagger c$)

## Mean-Field Decomposition

The Hartree-Fock approximation replaces the four-fermion interaction with a quadratic form by applying **Wick's theorem**. For the interaction term, we perform the mean-field decomposition:

$$c^\dagger_i c_j c^\dagger_k c_l \approx \langle c^\dagger_i c_j \rangle c^\dagger_k c_l + c^\dagger_i c_j \langle c^\dagger_k c_l \rangle - \langle c^\dagger_i c_l \rangle c^\dagger_k c_j - c^\dagger_i c_l \langle c^\dagger_k c_j \rangle$$

Here we omit the constant term $\langle c^\dagger_i c_j \rangle \langle c^\dagger_k c_l \rangle$ as it contributes only to the total energy.

### Green's Function Notation

Define the one-particle Green's function (density matrix):

$$G_{ij} = \langle c^\dagger_i c_j \rangle$$

The mean-field Hamiltonian becomes:

$$H_{\text{int}}^{\text{MF}} = \sum_{ijkl} V_{ijkl} \left[ G_{ij} c^\dagger_k c_l + c^\dagger_i c_j G_{kl} - G_{il} c^\dagger_k c_j - c^\dagger_i c_l G_{kj} \right]$$

This consists of **four terms**:
1. **Hartree direct term 1**: $\sum_{ijkl} V_{ijkl} G_{ij} c^\dagger_k c_l$
2. **Hartree direct term 2**: $\sum_{ijkl} V_{ijkl} c^\dagger_i c_j G_{kl}$
3. **Fock exchange term 1**: $-\sum_{ijkl} V_{ijkl} G_{il} c^\dagger_k c_j$
4. **Fock exchange term 2**: $-\sum_{ijkl} V_{ijkl} c^\dagger_i c_l G_{kj}$

## Effective Single-Particle Hamiltonian

We now express the mean-field Hamiltonian as an effective single-particle problem:

$$H^{\text{eff}} = \sum_{ij} h_{ij}^{\text{eff}} c^\dagger_i c_j$$

where the effective single-particle matrix is:

$$h_{ij}^{\text{eff}} = t_{ij} + \sum_{kl} U_{ijkl} \, G_{kl}$$

### Derivation of $U_{ijkl}$

To obtain a unified form, we rewrite each of the four mean-field terms to isolate the operator part $c^\dagger_i c_j$:

**Term 1 (Hartree)**: $\sum_{ijkl} V_{ijkl} G_{ij} c^\dagger_k c_l$

Relabel summation indices so that the operator becomes $c^\dagger_i c_j$:
- Operator: $c^\dagger_k c_l \to c^\dagger_i c_j$ implies $(k,l) \to (i,j)$
- Green's function: $G_{ij} \to G_{kl}$ implies $(i,j) \to (k,l)$

Result: $\sum_{ij} c^\dagger_i c_j \sum_{kl} V_{klij} G_{kl}$

**Term 2 (Hartree)**: $\sum_{ijkl} V_{ijkl} c^\dagger_i c_j G_{kl}$

Already in standard form:

Result: $\sum_{ij} c^\dagger_i c_j \sum_{kl} V_{ijkl} G_{kl}$

**Term 3 (Fock)**: $-\sum_{ijkl} V_{ijkl} G_{il} c^\dagger_k c_j$

Relabel to match $c^\dagger_i c_j$ and $G_{kl}$:
- Operator: $c^\dagger_k c_j \to c^\dagger_i c_j$ implies $k \to i$, $j$ unchanged
- Green's function: $G_{il} \to G_{kl}$ implies $i \to k$, $l$ unchanged

Result: $-\sum_{ij} c^\dagger_i c_j \sum_{kl} V_{kjil} G_{kl}$

**Term 4 (Fock)**: $-\sum_{ijkl} V_{ijkl} c^\dagger_i c_l G_{kj}$

Relabel to match $c^\dagger_i c_j$ and $G_{kl}$:
- Operator: $c^\dagger_i c_l \to c^\dagger_i c_j$ implies $l \to j$, $i$ unchanged
- Green's function: $G_{kj} \to G_{kl}$ implies $j \to l$, $k$ unchanged

Result: $-\sum_{ij} c^\dagger_i c_j \sum_{kl} V_{ilkj} G_{kl}$

### Final Formula

Combining all four terms, we obtain the **effective interaction tensor**:

$$\boxed{U_{ijkl} = V_{ijkl} + V_{klij} - V_{kjil} - V_{ilkj}}$$

For each set of indices $(i,j,k,l)$, $U_{ijkl}$ collects contributions from four different elements of the original interaction tensor $V$:

**Physical interpretation**:
- **Hartree terms** ($V_{ijkl} + V_{klij}$): Direct Coulomb repulsion, treating other electrons as a static charge density
- **Fock terms** ($-V_{kjil} - V_{ilkj}$): Exchange interaction arising from quantum statistics (Pauli principle)

## Matrix Representation

### 4D Tensor to 2D Matrix

The effective interaction tensor $U_{ijkl}$ is a rank-4 tensor. For efficient computation, we reshape it into a 2D matrix:

$$U_{(i,j),(k,l)} = U_{ijkl}$$

where the composite index $(i,j)$ runs over all pairs of single-particle states. For a system with $N$ single-particle states:
- 4D tensor $U_{ijkl}$: shape $(N, N, N, N)$
- 2D matrix $U_{(i,j),(k,l)}$: shape $(N^2, N^2)$

### Effective Hamiltonian Construction

In the self-consistent field iteration, we compute the effective Hamiltonian at each step:

$$H^{\text{eff}} = H_0 + U \cdot G$$

where:
- $H_0$: one-body Hamiltonian (shape $N \times N$)
- $U$: effective interaction matrix $U_{(i,j),(k,l)}$ (shape $N^2 \times N^2$)
- $G$: Green's function reshaped to a vector (length $N^2$)
- Result: effective Hamiltonian matrix (shape $N \times N$)

## Self-Consistent Field (SCF) Iteration

The Hartree-Fock equations are solved iteratively:

1. **Initialize**: Start with an initial guess for $G_{ij}$ (random or from previous calculation)

2. **Build effective Hamiltonian**:
   $$h^{\text{eff}} = t + U \cdot G$$

3. **Diagonalize**: Find eigenvalues $\{\varepsilon_n\}$ and eigenvectors $\{|n\rangle\}$ of $h^{\text{eff}}$

4. **Update Green's function**:
   - **Zero temperature**: Occupy lowest $N_e$ states
     $$G_{ij} = \sum_{n=1}^{N_e} \langle i | n \rangle \langle n | j \rangle$$

   - **Finite temperature**: Use Fermi-Dirac distribution with chemical potential $\mu$
     $$G_{ij} = \sum_n f(\varepsilon_n - \mu) \langle i | n \rangle \langle n | j \rangle$$
     where $f(\varepsilon) = 1/(e^{\varepsilon/T} + 1)$ and $\mu$ is determined by particle number conservation

5. **Check convergence**:
   $$\text{Rest} = \frac{1}{N^2} \|G_{\text{new}} - G_{\text{old}}\|$$

   If $\text{Rest} < \epsilon$, stop. Otherwise, mix and return to step 2:
   $$G \leftarrow (1-\alpha) G_{\text{old}} + \alpha G_{\text{new}}$$

6. **Calculate physical quantities**:
   - **Band energy**: $E_{\text{band}} = \sum_i H^{\text{eff}}_{ii} G_{ii}$ (zero T) or use grand potential (finite T)
   - **Interaction energy**: $E_{\text{int}} = -\frac{1}{2} G^T \cdot U \cdot G$
   - **Total energy**: $E_{\text{total}} = E_{\text{band}} + E_{\text{int}}$
   - **Particle number**: $N = \sum_i G_{ii}$
   - **Spin polarization**: $S_z = \frac{1}{2} \sum_i (G_{i\uparrow,i\uparrow} - G_{i\downarrow,i\downarrow})$

## Implementation Notes

### Recommended: Direct U Matrix Construction

**For efficient Hartree-Fock calculations, use `build_U_matrix` to directly construct the sparse U matrix from operators:**

```julia
U_matrix = build_U_matrix(dofs, interaction_ops, blocks)
```

This approach:
- **Skips the intermediate V tensor** (saves N⁴ memory, e.g., 1.3 TB for 30×30 system)
- **Directly generates sparse matrix** representation
- **Applies the 4-term formula** during construction
- **Exploits block structure** when provided

For each interaction operator $V_{ijkl}$, the function contributes to **four U matrix elements**:

```julia
U[(i-1)*N+j, (k-1)*N+l] += V[i,j,k,l]  # Hartree term 2
U[(k-1)*N+l, (i-1)*N+j] += V[i,j,k,l]  # Hartree term 1  (transpose)
U[(k-1)*N+j, (i-1)*N+l] -= V[i,j,k,l]  # Fock term 3     (exchange)
U[(i-1)*N+l, (k-1)*N+j] -= V[i,j,k,l]  # Fock term 4     (exchange)
```

### Alternative: Two-Step Construction via V Tensor

For debugging or small systems, you can use the two-step approach:

```julia
V = build_interaction_tensor(dofs, interaction_ops)  # Dense N×N×N×N tensor
# Then manually apply 4-term formula to build U
```

**Note**: This approach requires O(N⁴) memory and is **not recommended for large systems**.

### Symmetry and Blocking

When the system has conserved quantum numbers (e.g., spin $S_z$, particle number), the effective Hamiltonian and Green's function can be block-diagonalized:

- Use `SystemDofs` with `sortrule = [[conserved_dof_indices], other_dofs]`
- `blocks` field provides index ranges for each symmetry sector
- Diagonalize each block independently for efficiency
- Drastically reduces computational cost for large systems

### Block Structure and U Matrix Sparsification

**Important**: The `blocks` parameter optimizes the U matrix construction by exploiting the block-diagonal structure of G, **not** by making U itself block-diagonal.

**Key insight**: If G is block-diagonal (e.g., spin-up and spin-down blocks), then $G_{kl} = 0$ when $k$ and $l$ are in different blocks. When computing:

$$h^{\text{eff}}_{ij} = t_{ij} + \sum_{kl} U_{ijkl} G_{kl}$$

any $U_{ijkl}$ term multiplied by $G_{kl} = 0$ contributes nothing. Therefore, we can **skip storing** U matrix elements where $(k,l)$ are in different blocks.

**Storage reduction**: For B equal-sized blocks, this saves approximately $(1 - 1/B^2)$ of memory and computation:
- 2 spin blocks: **75% reduction**
- 4 blocks: **94% reduction**

**Implementation**: `build_U_matrix(dofs, ops, blocks)` only stores $U_{(i,j),(k,l)}$ when $(k,l)$ are in the same block, creating a column-sparse (not block-diagonal) matrix optimized for multiplication with block-diagonal G.

## Example: Hubbard Model

For the single-band Hubbard model:

$$H = -t \sum_{\langle ij \rangle, \sigma} c^\dagger_{i\sigma} c_{j\sigma} + U \sum_i n_{i\uparrow} n_{i\downarrow}$$

The on-site interaction $U n_{i\uparrow} n_{i\downarrow} = U c^\dagger_{i\uparrow} c_{i\uparrow} c^\dagger_{i\downarrow} c_{i\downarrow}$ is implemented as:

```julia
# Setup system with spin blocking
dofs = SystemDofs([Dof(:site, N), Dof(:spin, 2)], sortrule = [[2], 1])
onsite = bonds(lattice, (:o, :o), 0)  # On-site bonds

# Generate Hubbard U operators
U_ops = generate_twobody(dofs, onsite,
    (delta, qn1, qn2, qn3, qn4) ->
        (qn1.site == qn2.site == qn3.site == qn4.site) &&
        (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1, 1, 2, 2) ? U : 0.0,
    order = (cdag, 1, c, 1, cdag, 1, c, 1)
)

# Build U matrix directly (recommended for large systems)
U_matrix = build_U_matrix(dofs, U_ops, dofs.blocks)

# In SCF iteration:
# H_eff = H_0 + reshape(U_matrix * vec(G), N, N)
```

The `build_U_matrix` function automatically applies the mean-field 4-term formula and exploits the spin block structure for efficiency.

## References

1. Hartree, D. R. (1928). "The Wave Mechanics of an Atom with a Non-Coulomb Central Field"
2. Fock, V. (1930). "Näherungsmethode zur Lösung des quantenmechanischen Mehrkörperproblems"