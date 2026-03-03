"""
Test suite entry point for MeanFieldTheories.jl

This file aggregates all test files in the test/ directory.
Each source file should have a corresponding test_*.jl file.
"""

using Test
using SafeTestsets

@testset "MeanFieldTheories.jl" begin
    @safetestset "Quantum System" begin
        include("test_quantum_system.jl")
    end

    @safetestset "Hartree-Fock (real space)" begin
        include("test_hartreefock_real.jl")
    end
end
