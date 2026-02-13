"""
Tests for Hartree-Fock module (groundstate/hartreefock.jl)
"""

using Test
using LinearAlgebra
using SparseArrays
using SingleModeApproximation

@testset "build_U_matrix" begin
    # Setup a simple 2-site, 2-spin system
    dofs = SystemDofs([Dof(:site, 2), Dof(:spin, 2)], sortrule = [[2], 1])
    N = total_dim(dofs)  # Should be 4

    # Create a simple 2D lattice
    unitcell = Lattice([Dof(:site, 1)], [QN(site=1)], [[0.0, 0.0]])
    lattice = Lattice(unitcell, [[1.0, 0.0], [0.0, 1.0]], (2, 1))
    onsite_bonds = bonds(lattice, (:o, :o), 0)

    # Simple Hubbard U interaction: n_↑ n_↓
    U_ops = generate_twobody(dofs, onsite_bonds,
        (delta, qn1, qn2, qn3, qn4) ->
            (qn1.site == qn2.site == qn3.site == qn4.site) &&
            (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1, 1, 2, 2) ? 1.0 : 0.0,
        order = (cdag, 1, c, 1, cdag, 1, c, 1))

    @testset "Basic functionality" begin
        U_matrix = build_U_matrix(dofs, U_ops, dofs.blocks)

        @test size(U_matrix) == (N^2, N^2)
        @test U_matrix isa SparseMatrixCSC
        @test nnz(U_matrix) > 0
        @test all(isreal, nonzeros(U_matrix))
    end

    @testset "Block optimization" begin
        U_with_blocks = build_U_matrix(dofs, U_ops, dofs.blocks)
        U_without_blocks = build_U_matrix(dofs, U_ops, nothing)

        # With blocks should have fewer or equal non-zeros
        @test nnz(U_with_blocks) <= nnz(U_without_blocks)
    end

    @testset "Consistency with V tensor method" begin
        # Compare with reference implementation
        V = build_interaction_tensor(dofs, U_ops)

        # Apply 4-term formula manually
        U_ref = zeros(ComplexF64, N, N, N, N)
        for i in 1:N, j in 1:N, k in 1:N, l in 1:N
            U_ref[i,j,k,l] = V[i,j,k,l] + V[k,l,i,j] - V[k,j,i,l] - V[i,l,k,j]
        end

        U_matrix = build_U_matrix(dofs, U_ops, nothing)

        @test Matrix(U_matrix) ≈ reshape(U_ref, N^2, N^2)
    end

    @testset "Error handling" begin
        @test_throws ErrorException build_U_matrix(dofs, Operators[], dofs.blocks)
    end
end
