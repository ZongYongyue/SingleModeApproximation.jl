#!/usr/bin/env julia

"""
Script to run Hartree-Fock benchmarks with proper environment setup
"""

using Pkg

# Activate benchmark environment
println("Activating benchmark environment...")
Pkg.activate(@__DIR__)

# Install dependencies if needed
println("Installing/updating dependencies...")
Pkg.instantiate()

# Add parent package in development mode
println("Adding SingleModeApproximation package...")
parent_dir = dirname(@__DIR__)
Pkg.develop(path=parent_dir)

println("\n" * "=" ^ 80)
println("Running benchmarks...")
println("=" ^ 80 * "\n")

# Include and run benchmark
include("benchmark_hartreefock.jl")

# Test direct U matrix construction
test_build_U_matrix()

println("\n")


