# Momentum-Space Hartree-Fock: Derivation

## Step 1: General Two-Body Hamiltonian

The most general two-body interaction Hamiltonian in InterAll format is:

$$H_{\text{int}} = \frac{1}{2}\sum_{ijkl}\sum_{abcd} V_{ijkl}^{abcd} \, c^\dagger_{ia} c_{jb} c^\dagger_{kc} c_{ld}$$

Here $i,j,k,l$ are lattice site indices, and $a,b,c,d$ are **intra-site degrees of freedom** (a collective label for spin, orbital, or any other internal quantum numbers).

---

## Step 2: Reduction via Translational Invariance

### Vertex Locality

For physical Coulomb-type interactions, charge is conserved at each interaction vertex — the interaction does not directly mediate inter-site hopping. This means each $c^\dagger c$ pair acts on the **same site**:

$$j = i, \qquad l = k$$

The first operator pair $c^\dagger_{ia} c_{ib}$ sits at site $i$; the second $c^\dagger_{kc} c_{kd}$ sits at site $k$. This encompasses all common interaction types:

| Interaction | Orbital structure | Site structure |
|-------------|-------------------|----------------|
| Hubbard $U n_{a\uparrow}n_{a\downarrow}$ | $a=b=c=d$ | $i=k$ (on-site) |
| Inter-site Coulomb $V n_i n_j$ | $a=b$, $c=d$ | $i \neq k$ |
| Hund, Pair-Lift, Exchange | general $a,b,c,d$ | $i=k$ or $i \neq k$ |

Note that $a,b,c,d$ are **four independent internal indices** — the structure can be $aaaa$, $aabb$, $abab$, or fully general $abcd$ depending on the interaction type.

### Translational Invariance

For a translationally invariant system, $V_{ii,kk}^{abcd}$ (after imposing vertex locality $j=i$, $l=k$) depends only on the **relative displacement** $\mathbf{r} = \mathbf{R}_k - \mathbf{R}_i$:

$$V_{ii,kk}^{abcd} = W^{abcd}(\mathbf{R}_k - \mathbf{R}_i) \equiv W^{abcd}(\mathbf{r})$$

Renaming $k \to j$ for clarity, the reduced Hamiltonian is:

$$\boxed{H_{\text{int}} = \frac{1}{2}\sum_{i,j}\sum_{abcd} W^{abcd}(\mathbf{R}_j - \mathbf{R}_i) \, c^\dagger_{ia} c_{ib} c^\dagger_{jc} c_{jd}}$$

where $W^{abcd}(\mathbf{r})$ depends only on the inter-site distance. For $i=j$, $W^{abcd}(\mathbf{0})$ reduces to the local multi-orbital interaction matrix.

---

## Step 3: Fourier Transform to Momentum Space

For a lattice with $N$ sites and periodic boundary conditions, define Bloch creation/annihilation operators:

$$c^\dagger_{ia} = \frac{1}{\sqrt{N}}\sum_{\mathbf{k}} e^{-i\mathbf{k}\cdot\mathbf{R}_i} c^\dagger_{\mathbf{k}a}, \qquad c_{ia} = \frac{1}{\sqrt{N}}\sum_{\mathbf{k}} e^{+i\mathbf{k}\cdot\mathbf{R}_i} c_{\mathbf{k}a}$$

Substituting all four operators into $H_{\text{int}}$ and introducing momentum labels $\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3,\mathbf{k}_4$:

$$H_{\text{int}} = \frac{1}{2N^2}\sum_{ij,abcd}\sum_{\mathbf{k}_1\mathbf{k}_2\mathbf{k}_3\mathbf{k}_4} W^{abcd}(\mathbf{R}_j - \mathbf{R}_i)\, e^{-i\mathbf{k}_1\cdot\mathbf{R}_i} e^{+i\mathbf{k}_2\cdot\mathbf{R}_i} e^{-i\mathbf{k}_3\cdot\mathbf{R}_j} e^{+i\mathbf{k}_4\cdot\mathbf{R}_j} \, c^\dagger_{\mathbf{k}_1 a} c_{\mathbf{k}_2 b} c^\dagger_{\mathbf{k}_3 c} c_{\mathbf{k}_4 d}$$

**Summing over site $i$**, with $\mathbf{r} = \mathbf{R}_j - \mathbf{R}_i$:

$$\sum_i e^{i(-\mathbf{k}_1+\mathbf{k}_2-\mathbf{k}_3+\mathbf{k}_4)\cdot\mathbf{R}_i} = N\,\delta_{\mathbf{k}_1-\mathbf{k}_2+\mathbf{k}_3-\mathbf{k}_4,\,\mathbf{0}}$$

This enforces **momentum conservation**:

$$\mathbf{k}_1 - \mathbf{k}_2 + \mathbf{k}_3 - \mathbf{k}_4 = \mathbf{0} \pmod{\text{reciprocal lattice vector}}$$

**Summing over the relative displacement $\mathbf{r}$** defines the Fourier transform of the interaction:

$$\tilde{W}^{abcd}(\mathbf{q}) \equiv \sum_\mathbf{r} W^{abcd}(\mathbf{r})\,e^{i\mathbf{q}\cdot\mathbf{r}}$$

The overall prefactor becomes $\frac{1}{2N^2} \cdot N = \frac{1}{2N}$.

### Introducing the Momentum Transfer

Define the momentum transfer $\mathbf{q} = \mathbf{k}_1 - \mathbf{k}_2 = \mathbf{k}_4 - \mathbf{k}_3$ (momentum conservation guarantees equality). Taking $\mathbf{k}_2$ and $\mathbf{k}_4$ as the two independent free momenta:

$$\mathbf{k}_1 = \mathbf{k}_2 + \mathbf{q}, \qquad \mathbf{k}_3 = \mathbf{k}_4 - \mathbf{q}$$

The interaction Hamiltonian in momentum space becomes:

$$\boxed{H_{\text{int}} = \frac{1}{2N}\sum_{abcd}\sum_{\mathbf{k}_2,\mathbf{k}_4,\mathbf{q}} \tilde{W}^{abcd}(\mathbf{q})\, c^\dagger_{\mathbf{k}_2+\mathbf{q},a}\; c_{\mathbf{k}_2,b}\; c^\dagger_{\mathbf{k}_4-\mathbf{q},c}\; c_{\mathbf{k}_4,d}}$$

The four operators carry momenta $\mathbf{k}_2+\mathbf{q}$, $\mathbf{k}_2$, $\mathbf{k}_4-\mathbf{q}$, $\mathbf{k}_4$ respectively. The Hamiltonian still depends on **three independent momenta** and cannot yet be cast as a single-particle problem.

---

## Step 4: Hartree-Fock Mean-Field Decoupling

### Wick's Theorem

Applying Wick's theorem to the four-operator product (dropping the constant fully-contracted term):

$$c^\dagger_{\mathbf{k}_2+\mathbf{q},a}\; c_{\mathbf{k}_2,b}\; c^\dagger_{\mathbf{k}_4-\mathbf{q},c}\; c_{\mathbf{k}_4,d}$$
$$\approx \underbrace{\langle c^\dagger_{\mathbf{k}_2+\mathbf{q},a} c_{\mathbf{k}_2,b}\rangle}_{\mathcal{G}_1} c^\dagger_{\mathbf{k}_4-\mathbf{q},c} c_{\mathbf{k}_4,d}
  + c^\dagger_{\mathbf{k}_2+\mathbf{q},a} c_{\mathbf{k}_2,b}\underbrace{\langle c^\dagger_{\mathbf{k}_4-\mathbf{q},c} c_{\mathbf{k}_4,d}\rangle}_{\mathcal{G}_2}$$
$$- \underbrace{\langle c^\dagger_{\mathbf{k}_2+\mathbf{q},a} c_{\mathbf{k}_4,d}\rangle}_{\mathcal{G}_3} c^\dagger_{\mathbf{k}_4-\mathbf{q},c} c_{\mathbf{k}_2,b}
  - c^\dagger_{\mathbf{k}_2+\mathbf{q},a} c_{\mathbf{k}_4,d}\underbrace{\langle c^\dagger_{\mathbf{k}_4-\mathbf{q},c} c_{\mathbf{k}_2,b}\rangle}_{\mathcal{G}_4}$$

The first two terms (direct contractions) give the **Hartree** contributions; the last two (cross contractions) give the **Fock** contributions, with minus signs from fermionic statistics.

### Translational Invariance of the Green's Function

For a ground state with discrete translational symmetry, the single-particle Green's function is **diagonal in $\mathbf{k}$**:

$$\langle c^\dagger_{\mathbf{k},a} c_{\mathbf{k}',b}\rangle = G_{ab}(\mathbf{k})\,\delta_{\mathbf{k},\mathbf{k}'}$$

States at different momenta are uncorrelated, a direct consequence of the commutation of the translation operator with the Hamiltonian.

**Important note on symmetry breaking**: Although the Hamiltonian has translational symmetry, the HF ground state may spontaneously break this symmetry (e.g., antiferromagnetic order, charge density waves). In such cases, one must redefine the **unit cell** to the enlarged magnetic/modulated cell that restores periodicity. From this point forward, $\mathbf{k}$ refers to momenta in the reciprocal lattice of this (possibly enlarged) cell, and indices $a,b,c,d$ label all internal degrees of freedom within the chosen cell (dimension $d \to nd$ for an $n$-fold enlargement, with correspondingly fewer $\mathbf{k}$ points).

Substituting into each contraction:

$$\mathcal{G}_1 = G_{ab}(\mathbf{k}_2)\,\delta_{\mathbf{q},\mathbf{0}}, \qquad \mathcal{G}_2 = G_{cd}(\mathbf{k}_4)\,\delta_{\mathbf{q},\mathbf{0}}$$

$$\mathcal{G}_3 = G_{ad}(\mathbf{k}_4)\,\delta_{\mathbf{k}_2+\mathbf{q},\mathbf{k}_4}, \qquad \mathcal{G}_4 = G_{cb}(\mathbf{k}_2)\,\delta_{\mathbf{k}_4-\mathbf{q},\mathbf{k}_2}$$

---

## Step 5: Hartree Self-Energy — Only $\mathbf{q}=0$ Contributes

The factors $\delta_{\mathbf{q},\mathbf{0}}$ in $\mathcal{G}_1$ and $\mathcal{G}_2$ show that **all $\mathbf{q} \neq 0$ Hartree contributions vanish**.

**From $\mathcal{G}_1$** (setting $\mathbf{q}=0$; remaining operator: $c^\dagger_{\mathbf{k}_4,c}c_{\mathbf{k}_4,d}$):

$$H^{H1} = \frac{1}{2N}\sum_{abcd}\sum_{\mathbf{k}_2,\mathbf{k}_4} \tilde{W}^{abcd}(0)\, G_{ab}(\mathbf{k}_2)\; c^\dagger_{\mathbf{k}_4,c} c_{\mathbf{k}_4,d}$$

Summing over $\mathbf{k}_2$ yields the **on-site density matrix**:

$$G_{ab}^{(0)} \equiv G_{ab}(\mathbf{r}=0) = \frac{1}{N}\sum_{\mathbf{k}_2} G_{ab}(\mathbf{k}_2)$$

so:

$$H^{H1} = \frac{1}{2}\sum_{\mathbf{k}}\sum_{cd} \left(\sum_{ab}\tilde{W}^{abcd}(0)\,G_{ab}^{(0)}\right) c^\dagger_{\mathbf{k},c} c_{\mathbf{k},d}$$

**From $\mathcal{G}_2$** (setting $\mathbf{q}=0$; remaining operator: $c^\dagger_{\mathbf{k}_2,a}c_{\mathbf{k}_2,b}$):

$$H^{H2} = \frac{1}{2}\sum_{\mathbf{k}}\sum_{ab} \left(\sum_{cd}\tilde{W}^{abcd}(0)\,G_{cd}^{(0)}\right) c^\dagger_{\mathbf{k},a} c_{\mathbf{k},b}$$

Both contributions are $\mathbf{k}$-independent. Collecting into the $c^\dagger_{\mathbf{k},a}c_{\mathbf{k},b}$ channel, the **Hartree self-energy matrix** is:

$$\boxed{\Sigma^H_{ab} = \frac{1}{2}\sum_{cd}\left[\tilde{W}^{ab,cd}(0) + \tilde{W}^{cd,ab}(0)\right] G_{cd}^{(0)}}$$

**Physical picture**: $\Sigma^H$ is a constant matrix, identical for all $\mathbf{k}$ points, determined entirely by the local density matrix $G^{(0)}$. Here $\tilde{W}^{abcd}(0) = \sum_\mathbf{r} W^{abcd}(\mathbf{r})$ is the interaction summed over all inter-site distances.

---

## Step 6: Fock Self-Energy — Convolution in $\mathbf{k}$ Space

The constraint $\delta_{\mathbf{k}_2+\mathbf{q},\mathbf{k}_4}$ from $\mathcal{G}_3$ fixes $\mathbf{k}_4 = \mathbf{k}_2+\mathbf{q}$, so $\mathbf{k}_4-\mathbf{q} = \mathbf{k}_2$: Op 3 and Op 2 acquire the same momentum, and the remaining operator product becomes $-c^\dagger_{\mathbf{k}_2,c}\,c_{\mathbf{k}_2,b}$ (the minus sign from the Fock cross-contraction).

**From $\mathcal{G}_3$**:

$$H^{F1} = -\frac{1}{2N}\sum_{abcd}\sum_{\mathbf{k}_2,\mathbf{q}} \tilde{W}^{abcd}(\mathbf{q})\,G_{ad}(\mathbf{k}_2+\mathbf{q})\; c^\dagger_{\mathbf{k}_2,c}\, c_{\mathbf{k}_2,b}$$

**From $\mathcal{G}_4$** (same constraint $\mathbf{k}_4 = \mathbf{k}_2+\mathbf{q}$; remaining operator: $-c^\dagger_{\mathbf{k}_2+\mathbf{q},a}\,c_{\mathbf{k}_2+\mathbf{q},d}$; relabeling $\mathbf{k}_2+\mathbf{q} \to \mathbf{k}$):

$$H^{F2} = -\frac{1}{2N}\sum_{abcd}\sum_{\mathbf{k},\mathbf{q}} \tilde{W}^{abcd}(-\mathbf{q})\,G_{cb}(\mathbf{k}-\mathbf{q})\; c^\dagger_{\mathbf{k},a}\, c_{\mathbf{k},d}$$

Collecting both Fock terms into the $c^\dagger_{\mathbf{k},a}c_{\mathbf{k},b}$ channel:

**From $H^{F1}$**: the remaining operator is $c^\dagger_{\mathbf{k},c}\,c_{\mathbf{k},b}$, so the channel indices are $a\leftarrow c$ (3rd position in $\tilde{W}$) and $b\leftarrow b$ (2nd position), while $a,d$ are summed via $G_{ad}$. Relabelling the internal summation $a\to c,\,d\to d$:

$$\Sigma^{F1}_{ab}(\mathbf{k}) = -\frac{1}{2N}\sum_\mathbf{q}\sum_{cd}\tilde{W}^{cbad}(\mathbf{q})\,G_{cd}(\mathbf{k}+\mathbf{q})$$

**From $H^{F2}$**: the remaining operator is $c^\dagger_{\mathbf{k},a}\,c_{\mathbf{k},d}$, so $a\leftarrow a$ (1st position) and $b\leftarrow d$ (4th position), while $b,c$ are summed via $G_{cb}$. Substituting $\mathbf{q}\to-\mathbf{q}$ to bring $G$ to $G(\mathbf{k}+\mathbf{q})$ form and swapping the dummy pair $b\leftrightarrow c$:

$$\Sigma^{F2}_{ab}(\mathbf{k}) = -\frac{1}{2N}\sum_\mathbf{q}\sum_{cd}\tilde{W}^{adcb}(\mathbf{q})\,G_{cd}(\mathbf{k}+\mathbf{q})$$

The **Fock self-energy** is:

$$\boxed{\Sigma^F_{ab}(\mathbf{k}) = -\frac{1}{N}\sum_\mathbf{q}\sum_{cd}\,\frac{1}{2}\!\left[\tilde{W}^{cbad}(\mathbf{q}) + \tilde{W}^{adcb}(\mathbf{q})\right] G_{cd}(\mathbf{k}+\mathbf{q})}$$

This is the exact analogue of the Hartree formula (Step 5): two symmetry-related permutations of $\tilde{W}$ appear, one from each Wick contraction.

### Convolution Structure and FFT Acceleration

The Fock self-energy is a **convolution in $\mathbf{k}$ space**: $\Sigma^F(\mathbf{k}) \propto \sum_\mathbf{q} \tilde{W}(\mathbf{q})\,G(\mathbf{k}+\mathbf{q})$. By the convolution theorem:

$$\text{convolution in } \mathbf{k}\text{-space} \;\Longleftrightarrow\; \text{pointwise product in } \mathbf{r}\text{-space}$$

$$[\Sigma^F(\mathbf{k})]_{ab} = -\mathcal{F}_{\mathbf{r}\to\mathbf{k}}\!\left[\sum_{cd}\frac{1}{2}\!\left[W^{cbad}(\mathbf{r})+W^{adcb}(\mathbf{r})\right] G_{cd}(\mathbf{r})\right]$$

Concretely:
1. Compute the combined interaction kernel $\frac{1}{2}[W^{cbad}(\mathbf{r})+W^{adcb}(\mathbf{r})]$ and contract with $G_{cd}(\mathbf{r})$ pointwise in real space (summing over $c,d$)
2. Apply a single FFT to obtain $[\Sigma^F(\mathbf{k})]_{ab}$

This reduces the naive $O(N^2)$ convolution to $O(N\log N)$.

**Physical picture**: $\Sigma^F(\mathbf{k})$ depends on the full real-space Green's function $G(\mathbf{r})$ and captures non-local quantum exchange effects.

---

## Step 7: Effective Single-Particle Hamiltonian

Combining the one-body term (block-diagonal in $\mathbf{k}$ after Fourier transform) with the mean-field interaction, the effective Hamiltonian at each $\mathbf{k}$ point is:

$$\boxed{H^{\text{eff}}(\mathbf{k}) = T(\mathbf{k}) + \Sigma^H + \Sigma^F(\mathbf{k})}$$

| Term | Expression | $\mathbf{k}$-dependence | Physical meaning |
|------|-----------|-------------------------|------------------|
| $T(\mathbf{k})$ | $\sum_\mathbf{r} T^{ab}(\mathbf{r})\,e^{i\mathbf{k}\cdot\mathbf{r}}$ | $\mathbf{k}$-dependent | Single-particle hopping (band structure) |
| $\Sigma^H$ | $\sum_{cd}\bar{W}^{ab,cd}(0)\,G_{cd}^{(0)}$ | **$\mathbf{k}$-independent** | Hartree: mean-field Coulomb repulsion |
| $\Sigma^F(\mathbf{k})$ | $-\mathcal{F}[W(\mathbf{r})\cdot G(\mathbf{r})](\mathbf{k})$ | $\mathbf{k}$-dependent | Fock: quantum exchange interaction |

The original four-body interaction depending on three independent momenta has been reduced, via HF decoupling and translational symmetry, to $N$ independent $d\times d$ eigenvalue problems (where $d$ is the number of intra-site degrees of freedom). This is the fundamental efficiency gain of momentum-space Hartree-Fock.

---

## Self-Consistent Field (SCF) Iteration

1. **Initialize** $G_{ab}(\mathbf{k})$ (random, zero, or from a previous calculation)
2. **Build** $H^{\text{eff}}(\mathbf{k}) = T(\mathbf{k}) + \Sigma^H[G^{(0)}] + \Sigma^F[G(\mathbf{r})]$
3. **Diagonalize** $H^{\text{eff}}(\mathbf{k})$ at each $\mathbf{k}$ point: $H^{\text{eff}}(\mathbf{k})\,v_{n\mathbf{k}} = \varepsilon_{n\mathbf{k}}\,v_{n\mathbf{k}}$
4. **Determine occupation** $f_{n\mathbf{k}}$: step function at $T=0$, or Fermi-Dirac at $T>0$
5. **Update Green's function**: $G_{ab}(\mathbf{k}) = \sum_n v^*_{an}(\mathbf{k})\,f_{n\mathbf{k}}\,v_{bn}(\mathbf{k})$, then FFT to real space
6. **Check convergence**: $\|G_{\text{new}} - G_{\text{old}}\| < \varepsilon$
7. **Mix**: $G \leftarrow (1-\alpha)G_{\text{old}} + \alpha G_{\text{new}}$, go to step 2

---
---

# Momentum-Space Hartree-Fock: Implementation

This section describes how to implement the derivation above as a concrete algorithm. We follow the same logic as H-wave's `UHFk` solver, but abstract away specific data formats so that the procedure applies to any system with translational symmetry and arbitrary internal degrees of freedom.

## Input Specification

A momentum-space HF calculation requires three groups of inputs:

### 1. Lattice and Unit Cell

| Parameter | Meaning |
|-----------|---------|
| **Cell shape** $(L_x, L_y, L_z)$ | Total lattice size (number of unit cells along each direction) |
| **Subcell shape** $(B_x, B_y, B_z)$ | Magnetic/modulated unit cell size (for symmetry breaking; $B=1$ if no enlargement) |
| **Geometry** | Lattice vectors and atomic positions within the unit cell |

The **effective k-grid** has $(N_x, N_y, N_z) = (L_x/B_x,\, L_y/B_y,\, L_z/B_z)$ points with $N_k = N_x N_y N_z$. The number of internal degrees of freedom per cell is $d = n_{\mathrm{orb}} \times B_x B_y B_z \times n_{\mathrm{spin}}$ (orbital count is multiplied by the subcell volume because the enlarged cell absorbs extra sites as internal degrees of freedom).

**Example** (4×4 Hubbard model with 2×2 antiferromagnetic subcell):
- Cell shape: $(4, 4, 1)$, Subcell shape: $(2, 2, 1)$
- k-grid: $(2, 2, 1) \Rightarrow N_k = 4$
- Original: 1 orbital, 2 spins → $d = 1 \times 4 \times 2 = 8$ internal DOFs per magnetic cell

### 2. One-Body Terms: Hopping $T_{ab}(\mathbf{r})$

The hopping integral between internal DOFs $a$ and $b$ separated by lattice vector $\mathbf{r}$. Stored as a real-space array:

$$T_{ab}(\mathbf{r}), \qquad \mathbf{r} \in \{0, \ldots, N_x-1\} \times \{0, \ldots, N_y-1\} \times \{0, \ldots, N_z-1\}$$

**Example** (square lattice nearest-neighbor hopping $t=1$):

| $\mathbf{r}$ | $a$ | $b$ | $T_{ab}(\mathbf{r})$ |
|:---:|:---:|:---:|:---:|
| $(1,0,0)$ | 1 | 1 | $-1.0$ |
| $(0,1,0)$ | 1 | 1 | $-1.0$ |
| $(-1,0,0)$ | 1 | 1 | $-1.0$ |
| $(0,-1,0)$ | 1 | 1 | $-1.0$ |

### 3. Two-Body Terms: Interaction $W^{abcd}(\mathbf{r})$

In practice, interactions are specified in terms of physically motivated types (Coulomb, Hund, Exchange, etc.), each of which maps to a specific pattern of $W^{abcd}$. To handle this efficiently, we decompose:

$$W^{abcd}(\mathbf{r}) = \sum_\lambda J^\lambda_{ab}(\mathbf{r}) \cdot S^\lambda_{stuv}$$

where $\lambda$ labels the interaction type, $J^\lambda_{ab}(\mathbf{r})$ encodes the orbital-spatial part (a 2-index interaction matrix at each $\mathbf{r}$), and $S^\lambda_{stuv}$ is a **spin table** that specifies which spin combinations contribute.

| Type | Physical form | Non-zero $S_{stuv}$ |
|------|-------------|---------------------|
| CoulombIntra ($U$) | $U\,n_{a\uparrow}n_{a\downarrow}$ | $(0,1,1,0)$, $(1,0,0,1)$ |
| CoulombInter ($V$) | $V\,n_{ia}n_{jb}$ | $(0,0,0,0)$, $(1,1,1,1)$, $(0,1,1,0)$, $(1,0,0,1)$ |
| Hund ($J_H$) | $-J_H \mathbf{S}_a\cdot\mathbf{S}_b$ | $(0,0,0,0)$, $(1,1,1,1)$ |
| Exchange ($J_{\mathrm{Ex}}$) | pair exchange | $(0,1,0,1)$, $(1,0,1,0)$ |

**Symmetrization**: Each $J^\lambda_{ab}(\mathbf{r})$ must be symmetrized as $\tilde{J}^\lambda_{ab}(\mathbf{r}) = \frac{1}{2}\!\left[J^\lambda_{ab}(\mathbf{r}) + J^{\lambda*}_{ba}(-\mathbf{r})\right]$ before use, which absorbs both Hartree terms (or both Fock terms) into a single computation.

---

## Preprocessing (Performed Once)

### Step A: Build $T(\mathbf{k})$ via FFT

*→ Implements the kinetic term in the Step 7 table: $T(\mathbf{k}) = \sum_\mathbf{r} T^{ab}(\mathbf{r})\,e^{i\mathbf{k}\cdot\mathbf{r}}$.*

Fourier-transform the real-space hopping to momentum space:

$$T_{ab}(\mathbf{k}) = \sum_\mathbf{r} T_{ab}(\mathbf{r})\,e^{i\mathbf{k}\cdot\mathbf{r}}$$

Computationally, this is a 3D FFT over the spatial indices $(r_x, r_y, r_z)$, applied independently for each pair $(a, b)$:

```
Input:  T_r[rx, ry, rz, a, b]     — real-space hopping
Output: T_k[kx, ky, kz, a, b]     — momentum-space hopping
Method: T_k = IFFT_3D(T_r, axes=(0,1,2), norm="forward")
```

If spin does not couple to the hopping (no spin-orbital interaction), the spin structure is trivially diagonal: $T(\mathbf{k})$ is block-diagonal with identical blocks for spin-up and spin-down.

### Step B: Prepare Symmetrized Interactions

*→ Absorbs the two-term structure of both the Hartree self-energy (Step 5: $[\tilde{W}^{ab,cd}(0) + \tilde{W}^{cd,ab}(0)]$) and the Fock self-energy (Step 6: $[\tilde{W}^{cbad}(\mathbf{q}) + \tilde{W}^{adcb}(\mathbf{q})]$) into a single symmetrized coefficient $\tilde{J}_{ab}(\mathbf{r})$, so each self-energy requires only one contraction instead of two.*

For each interaction type $\lambda$:

1. Load the real-space interaction matrix $J^\lambda_{ab}(\mathbf{r})$
2. Construct the Hermitian-conjugated reverse: $J^\lambda_{ba}(-\mathbf{r})^*$
3. Symmetrize: $\tilde{J}^\lambda_{ab}(\mathbf{r}) = \frac{1}{2}[J^\lambda_{ab}(\mathbf{r}) + J^{\lambda*}_{ba}(-\mathbf{r})]$
4. Store the spin table $S^\lambda_{stuv}$

After this step, we have a list of $\{(\tilde{J}^\lambda, S^\lambda)\}$ ready for the SCF loop.

---

## SCF Iteration (Repeated Until Convergence)

At the start of each iteration, we have the current Green's function $G_{ab}(\mathbf{r})$ in real space (shape: $[N_k, d, d]$ where the first axis is the r-space index). The following steps update it.

### Step 1: Construct $H^{\mathrm{eff}}(\mathbf{k})$

*→ Implements Step 7: $H^{\text{eff}}(\mathbf{k}) = T(\mathbf{k}) + \Sigma^H + \Sigma^F(\mathbf{k})$.*

Start from the pre-computed kinetic term:

$$H(\mathbf{k}) = T(\mathbf{k})$$

Then add the self-energy from each interaction type $\lambda$:

**Hartree contribution** (k-independent):

*→ Implements Step 5: $\Sigma^H_{ab} = \frac{1}{2}\sum_{cd}\left[\tilde{W}^{ab,cd}(0) + \tilde{W}^{cd,ab}(0)\right] G_{cd}^{(0)}$. The two-term sum $[\tilde{W}^{ab,cd}+\tilde{W}^{cd,ab}]$ is already encoded in the symmetrized $\tilde{J}$ (Step B). The density $G_{cd}^{(0)} = G_{cd}(\mathbf{r}=0)$ (Step 5) reduces to the on-site orbital-diagonal $G_{cc}(\mathbf{0})$ for density-density interactions.*

The Hartree self-energy uses only the on-site ($\mathbf{r}=0$) diagonal of G:

$$\Sigma^{H,\lambda}_{(s,a),(t,b)} = \delta_{ab} \sum_c \left[\sum_\mathbf{r} \tilde{J}^\lambda_{ac}(\mathbf{r})\right] \sum_{u,v} G_{cc}^{uv}(\mathbf{0})\, S^\lambda_{s,u,v,t}$$

where $G_{cc}^{uv}(\mathbf{0})$ is the on-site orbital-diagonal density matrix with spin indices $(u,v)$.

Concretely:
```
For each interaction type λ:
    gbb[u,v,c]   = G[r=0, u, c, v, c]          # on-site, orbital-diagonal
    hh0[s,t,c]   = sum_{u,v} gbb[u,v,c] * S[s,u,v,t]   # spin contraction
    hh1[s,t,a]   = sum_c J_tilde[r,a,c] * hh0[s,t,c]    # interaction contraction, sum over r and c
    Σ_H[s,a,t,b] = δ_{ab} * hh1[s,t,a]          # orbital-diagonal in output
    H(k) += Σ_H   for all k                     # broadcast to all k-points
```

**Fock contribution** (k-dependent, via FFT):

*→ Implements Step 6's convolution-theorem form: $[\Sigma^F(\mathbf{k})]_{ab} = -\mathcal{F}_{\mathbf{r}\to\mathbf{k}}\!\left[\sum_{cd}\frac{1}{2}[W^{cbad}(\mathbf{r})+W^{adcb}(\mathbf{r})] G_{cd}(\mathbf{r})\right]$. The two index permutations of $W$ are absorbed into $\tilde{J}$ (Step B). Note that $G_{cd}(\mathbf{r})$ appears as $G_{ba}(\mathbf{r})$ in the code because the orbital indices are transposed relative to the Hartree case (Step 6: cross-contraction $\mathcal{G}_3, \mathcal{G}_4$).*

The Fock self-energy is computed as a pointwise product in real space followed by FFT:

$$\Sigma^{F,\lambda}_{(s,a),(t,b)}(\mathbf{r}) = -\tilde{J}^\lambda_{ab}(\mathbf{r}) \sum_{u,v} G^{uv}_{ba}(\mathbf{r})\, S^\lambda_{s,u,t,v}$$

$$\Sigma^{F,\lambda}(\mathbf{k}) = \mathrm{FFT}_{\mathbf{r}\to\mathbf{k}}\!\left[\Sigma^{F,\lambda}(\mathbf{r})\right]$$

Concretely:
```
For each interaction type λ:
    hh[r,s,a,t,b] = J_tilde[r,a,b] * sum_{u,v} G[r,u,b,v,a] * S[s,u,t,v]
    Σ_F(k)        = FFT_3D(hh, axes=(0,1,2))
    H(k) -= Σ_F(k)
```

Note the key difference: Hartree uses $G_{cc}(\mathbf{0})$ (on-site, orbital-diagonal) while Fock uses $G_{ba}(\mathbf{r})$ (off-site, orbital-transposed). Hartree is k-independent; Fock creates k-dependence through the FFT.

### Step 2: Diagonalize

*→ Implements SCF step 3: diagonalize $H^{\text{eff}}(\mathbf{k})\,v_{n\mathbf{k}} = \varepsilon_{n\mathbf{k}}\,v_{n\mathbf{k}}$ at each $\mathbf{k}$ independently. This is the central efficiency gain of Step 7: the four-body problem has been reduced to $N_k$ independent $d\times d$ eigenproblems.*

Diagonalize $H^{\mathrm{eff}}(\mathbf{k})$ at each k-point to obtain eigenvalues and eigenvectors:

$$H^{\mathrm{eff}}(\mathbf{k})\,|\psi_{n\mathbf{k}}\rangle = \varepsilon_{n\mathbf{k}}\,|\psi_{n\mathbf{k}}\rangle$$

If total $S_z$ is conserved (no spin-orbit coupling or spin-flip interactions), $H(\mathbf{k})$ is block-diagonal in spin and can be diagonalized as two smaller $n_{\mathrm{orb}} \times n_{\mathrm{orb}}$ problems per k-point, rather than one $d \times d$ problem.

### Step 3: Determine Occupation

*→ Implements SCF step 4: determine $f_{n\mathbf{k}}$ — step function at $T=0$, Fermi-Dirac at $T>0$.*

**At $T=0$**: fill the lowest $N_e$ states (summed over all k-points and spin blocks). The chemical potential $\mu$ lies in the gap between the $N_e$-th and $(N_e+1)$-th eigenvalue.

**At $T>0$**: solve for the chemical potential $\mu$ such that:

$$\sum_{n,\mathbf{k}} f(\varepsilon_{n\mathbf{k}}, \mu) = N_e, \qquad f(\varepsilon, \mu) = \frac{1}{1 + e^{(\varepsilon - \mu)/T}}$$

This is a root-finding problem (bisection or Newton's method).

### Step 4: Update Green's Function

*→ Implements SCF step 5: $G_{ab}(\mathbf{k}) = \sum_n v^*_{an}(\mathbf{k})\,f_{n\mathbf{k}}\,v_{bn}(\mathbf{k})$. The subsequent IFFT gives the real-space $G_{ab}(\mathbf{r})$ needed by the Fock term (Step 6's convolution structure) in the next iteration.*

Construct the new Green's function in k-space:

$$G_{ab}^{\mathrm{new}}(\mathbf{k}) = \sum_n \psi^*_{na}(\mathbf{k})\,f_{n\mathbf{k}}\,\psi_{nb}(\mathbf{k})$$

Then inverse-FFT to real space:

$$G_{ab}(\mathbf{r}) = \mathrm{FFT}_{\mathbf{k}\to\mathbf{r}}\!\left[G_{ab}(\mathbf{k})\right]$$

This real-space $G(\mathbf{r})$ is needed for the Fock term in the next iteration.

### Step 5: Calculate Total Energy

*→ The band energy $E_{\mathrm{band}} = \sum_{n\mathbf{k}} f_{n\mathbf{k}}\,\varepsilon_{n\mathbf{k}}$ counts each occupied level once. However, $\varepsilon_{n\mathbf{k}}$ already includes the mean-field self-energy $\Sigma^H + \Sigma^F$, so summing eigenvalues double-counts the interaction energy. The interaction energy correction $E_{\mathrm{int}}$ removes this double-count. This double-counting issue arises because the HF Hamiltonian (Step 7) treats the interaction at the mean-field level, not at the full two-body level.*

**Band energy** (at $T=0$):

$$E_{\mathrm{band}} = \sum_{n,\mathbf{k}} f_{n\mathbf{k}}\,\varepsilon_{n\mathbf{k}}$$

At $T>0$, the free energy is:

$$\Omega_{\mathrm{band}} = \mu N_e - T \sum_{n,\mathbf{k}} \ln\!\left(1 + e^{-(\varepsilon_{n\mathbf{k}} - \mu)/T}\right)$$

**Interaction energy** (subtract double-counted mean-field):

For each interaction type $\lambda$:

$$E^{\lambda}_{\mathrm{int}} = -\frac{N_k}{2} \sum_\mathbf{r} \sum_{ab} \tilde{J}^\lambda_{ab}(\mathbf{r}) \left[\sum_{stuv} S^\lambda_{stuv}\,G^{va}_{sa}(\mathbf{0})\,G^{ub}_{tb}(\mathbf{0}) \;-\; \sum_{stuv} S^\lambda_{stuv}\,G^{ub}_{sa}(\mathbf{r})\,G^{va}_{tb}(\mathbf{r})\right]$$

The first term is the Hartree double-count (on-site), the second is the Fock double-count (all r).

**Total energy**: $E_{\mathrm{total}} = E_{\mathrm{band}} + \sum_\lambda E^\lambda_{\mathrm{int}}$

### Step 6: Check Convergence and Mix

*→ Implements SCF steps 6 and 7: check $\|G_{\mathrm{new}} - G_{\mathrm{old}}\| < \varepsilon$, then mix $G \leftarrow (1-\alpha)G_{\mathrm{old}} + \alpha G_{\mathrm{new}}$ and return to Step 1.*

Compute the residual:

$$\mathrm{rest} = \frac{\|G_{\mathrm{new}} - G_{\mathrm{old}}\|}{\mathrm{size}(G)}$$

If $\mathrm{rest} < \varepsilon$, the calculation has converged. Otherwise, mix and return to Step 1:

$$G \leftarrow (1 - \alpha)\,G_{\mathrm{old}} + \alpha\,G_{\mathrm{new}}$$

More advanced mixing strategies (DIIS, Broyden, scheduled mixing) can be used to improve convergence stability.

---

## Summary: Complete Algorithm

```
PREPROCESSING (once):
  A. T(k) = FFT[ T(r) ]                         — kinetic energy in k-space
  B. {J̃_λ(r), S_λ} = symmetrize interactions    — interaction coefficients

INITIALIZE:
  G(r) = random / zero / from previous calculation

SCF LOOP:
  ┌─ 1. Build H(k):
  │     H(k) = T(k)
  │     For each interaction λ:
  │       H(k) += Σ_H[λ]           ← k-independent, from G(r=0) diagonal
  │       H(k) -= FFT[ J̃(r)·G(r) ] ← k-dependent Fock term
  │
  ├─ 2. Diagonalize H(k) → {ε_nk, ψ_nk}
  │
  ├─ 3. Find μ, determine occupation f_nk
  │
  ├─ 4. G_new(k) = Σ_n ψ*ψ f_nk ;  G_new(r) = IFFT[ G_new(k) ]
  │
  ├─ 5. E_total = E_band + E_int
  │
  └─ 6. Converged? → done
         Not converged? → G = mix(G_old, G_new), go to 1
```

### Computational Complexity per SCF Iteration

| Operation | Cost | Note |
|-----------|------|------|
| Hartree self-energy | $O(N_k \cdot d)$ | Only r=0, orbital-diagonal |
| Fock self-energy | $O(N_k \log N_k \cdot d^2)$ | FFT over k-grid |
| Diagonalization | $O(N_k \cdot d^3)$ | Dominates for large d |
| Green's function update | $O(N_k \cdot d^2)$ | Matrix multiply + FFT |

The total cost per iteration scales as $O(N_k \cdot d^3)$, compared with $O(N^3)$ for real-space HF where $N = N_k \cdot d$. For large systems with small unit cells ($N_k \gg d$), the momentum-space method is dramatically cheaper.
