"""
Test suite entry point for SingleModeApproximation.jl

This file aggregates all test files in the test/ directory.
Each source file should have a corresponding test_*.jl file.
"""

using Test
using SafeTestsets

@testset "SingleModeApproximation.jl" begin
    @safetestset "Quantum System" begin
        include("test_quantum_system.jl")
    end

    # Future test files can be added here:
    # @safetestset "Hartree-Fock Solver" begin
    #     include("test_hartreefork.jl")
    # end
end
