"""
Hartree-Fock U*G computation using sparse matrix multiplication

Test system: Multi-orbital Hubbard model on square lattice
"""

using BenchmarkTools
using SparseArrays
using LinearAlgebra
using Random
using Printf

# Add parent directory to load path
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using SingleModeApproximation

#==================== Direct U Matrix Construction Test ====================#

function test_build_U_matrix()
    println("="^80)
    println("Direct U Matrix Construction Test (build_U_matrix)")
    println("="^80)

    # Test different sizes
    test_sizes = [
        (Lx=5, Ly=5, Norb=2),
        (Lx=10, Ly=10, Norb=2),
        (Lx=20, Ly=20, Norb=2),
        (Lx=30, Ly=30, Norb=2),
    ]

    println("\n| Lx×Ly | N    | Gen Ops (s) | Build U (s) | nnz(U) | Memory (MB) |")
    println("|-------|------|-------------|-------------|--------|-------------|")

    for (Lx, Ly, Norb) in test_sizes
        # Setup system
        Nsite = Lx * Ly
        unitcell = Lattice([Dof(:site, 1)], [QN(site=1)], [[0.0, 0.0]])
        lattice = Lattice(unitcell, [[1.0, 0.0], [0.0, 1.0]], (Lx, Ly))

        dofs = SystemDofs([
            Dof(:site, Nsite),
            Dof(:orbital, Norb),
            Dof(:spin, 2)
        ], sortrule = [[3], 2, 1])

        N = total_dim(dofs)
        blocks = dofs.blocks

        # Prepare bonds
        onsite_bonds = bonds(lattice, (:o, :o), 0)

        # Time operator generation
        t_ops = @elapsed begin
            U_ops = generate_twobody(dofs, onsite_bonds,
                (d, q1, q2, q3, q4) -> (q1.orbital == q2.orbital == q3.orbital == q4.orbital) &&
                                       (q1.spin, q2.spin, q3.spin, q4.spin) == (1, 1, 2, 2) ? 6.0 : 0.0,
                order = (cdag, 1, c, 1, cdag, 1, c, 1))
            Up_is = generate_twobody(dofs, onsite_bonds,
                (d, q1, q2, q3, q4) -> (q1.orbital, q3.orbital) == (q2.orbital, q4.orbital) &&
                                       (q1.orbital < q3.orbital) && (q1.spin, q3.spin) == (q2.spin, q4.spin) &&
                                       (q1.spin !== q3.spin) ? 3.6 : 0.0,
                order = (cdag, 1, c, 1, cdag, 1, c, 1))
            Up_ia = generate_twobody(dofs, onsite_bonds,
                (d, q1, q2, q3, q4) -> (q1.orbital, q3.orbital) == (q2.orbital, q4.orbital) &&
                                       (q1.orbital < q3.orbital) && (q1.spin, q3.spin) == (q2.spin, q4.spin) &&
                                       (q1.spin == q3.spin) ? 2.4 : 0.0,
                order = (cdag, 1, c, 1, cdag, 1, c, 1))
            J_sf = generate_twobody(dofs, onsite_bonds,
                (d, q1, q2, q3, q4) -> (q1.orbital, q2.orbital) == (q3.orbital, q4.orbital) &&
                                       (q1.orbital !== q2.orbital) &&
                                       (q1.spin, q2.spin, q3.spin, q4.spin) == (1, 2, 2, 1) ? 1.2 : 0.0,
                order = (cdag, 1, cdag, 1, c, 1, c, 1))
            J_ph = generate_twobody(dofs, onsite_bonds,
                (d, q1, q2, q3, q4) -> (q1.orbital, q3.orbital) == (q2.orbital, q4.orbital) &&
                                       (q1.orbital !== q3.orbital) &&
                                       (q1.spin, q2.spin, q3.spin, q4.spin) == (1, 2, 2, 1) ? 1.2 : 0.0,
                order = (cdag, 1, cdag, 1, c, 1, c, 1))
            all_ops = vcat(U_ops, Up_is, Up_ia, J_sf, J_ph)
        end

        # Time U matrix construction
        t_U = @elapsed U_matrix = build_U_matrix(dofs, all_ops, blocks)
        nnz_U = nnz(U_matrix)

        # Estimate memory usage (sparse matrix)
        mem_mb = (nnz_U * (2 * sizeof(Int) + sizeof(eltype(U_matrix)))) / 1024^2

        @printf("| %2d×%-2d | %-4d | %11.3f | %11.3f | %-6d | %11.2f |\n",
                Lx, Ly, N, t_ops, t_U, nnz_U, mem_mb)

        # Force garbage collection
        GC.gc()
    end

    println("="^80)
end
