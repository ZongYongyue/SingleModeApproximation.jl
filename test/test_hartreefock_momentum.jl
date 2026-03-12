"""
Tests for Hartree-Fock module (groundstate/hartreefock_momentum.jl)

Model: V-model on 2D square lattice with 2×2 AFM magnetic unit cell.
  H = -t Σ_{<ij>,σ} c†_{iσ}c_{jσ} + V Σ_{i≠j,σσ'} n_{iσ}n_{jσ'}
  t=1, V=4, half-filling (16 electrons on 4×4 sites).
"""

using Test
using MeanFieldTheories

const _t = 1.0
const _V = 4.0

# 2×2 magnetic unit cell: 4 sites × 2 spins → d=8
_dofs = SystemDofs([Dof(:site, 4), Dof(:spin, 2, [:up, :dn])])

_unitcell = Lattice([Dof(:site, 4)],
                    [QN(site=1), QN(site=2), QN(site=3), QN(site=4)],
                    [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
                    vectors=[[2.0, 0.0], [0.0, 2.0]])

_nn_bonds  = bonds(_unitcell, (:p, :p), 1)
_kpoints   = build_kpoints([[2.0, 0.0], [0.0, 2.0]], (2, 2))
_n_elec    = 16

_onebody = generate_onebody(_dofs, _nn_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -_t : 0.0)

# include 1/2 prefactor in the operator definition
_twobody = generate_twobody(_dofs, _nn_bonds,
    (deltas, qn1, qn2, qn3, qn4) ->
        qn1.spin == qn2.spin && qn3.spin == qn4.spin ? 0.5 * _V : 0.0)

_result = solve_hfk(_dofs, _onebody, _twobody, _kpoints, _n_elec;
    n_restarts=20, tol=1e-8, verbose=false)

@testset "solve_hfk — V-model 2×2 AFM cell" begin

    @testset "Convergence" begin
        @test _result.converged
    end

    @testset "Particle number" begin
        # ncond is electrons per magnetic unit cell; 16 electrons / 4 k-points = 4
        @test abs(_result.ncond - _n_elec / length(_kpoints)) < 1e-6
    end

    @testset "Energy" begin
        # Reference value from a converged SCF calculation.
        @test abs(_result.energies.total - (-0.5656600250)) < 1e-2
    end

    @testset "Output shapes" begin
        d  = 8   # 4 sites × 2 spins
        Nk = 4   # 2×2 k-grid
        @test size(_result.G_k)          == (d, d, Nk)
        @test size(_result.eigenvalues)  == (d, Nk)
        @test size(_result.eigenvectors) == (d, d, Nk)
        @test length(_result.kpoints)    == Nk
    end

    @testset "Green's function is Hermitian at each k" begin
        for ki in 1:size(_result.G_k, 3)
            G = _result.G_k[:, :, ki]
            @test G ≈ G'  atol=1e-10
        end
    end

    @testset "Eigenvalues are real" begin
        @test all(isreal.(_result.eigenvalues))
    end

end
