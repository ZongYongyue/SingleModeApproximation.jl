"""
Hartree-Fock mean-field approximation for quantum many-body systems.
"""

using SparseArrays
using LinearAlgebra

"""
    build_U_matrix(dofs::SystemDofs, ops::AbstractVector{<:Operators}, blocks=nothing)

Build sparse Hartree-Fock U matrix (N²×N²) directly from interaction operators.

This function efficiently constructs the U matrix without the intermediate V tensor,
applying the 4-term formula:

    U_{ijkl} = V_{ijkl} + V_{klij} - V_{kjil} - V_{ilkj}

where V is implicitly defined by the interaction operators.

# Arguments
- `dofs::SystemDofs`: System degrees of freedom
- `ops::AbstractVector{<:Operators}`: Two-body interaction operators (each with 4 FermionOp)
- `blocks`: Optional block structure. If provided, only stores elements where (k,l) are in the same block.

# Returns
- Sparse matrix U of size (N²×N²) where N = total_dim(dofs)

# Implementation details
For each operator corresponding to V_{ijkl}, it contributes to 4 U matrix elements:
1. U[(i-1)*N+j, (k-1)*N+l] += V_{ijkl}
2. U[(k-1)*N+l, (i-1)*N+j] += V_{ijkl}
3. U[(k-1)*N+j, (i-1)*N+l] -= V_{ijkl}
4. U[(i-1)*N+l, (k-1)*N+j] -= V_{ijkl}

# Performance
- Avoids allocating N⁴ dense tensor for V
- Uses dictionary accumulation to handle overlapping contributions
- Pre-allocates sparse matrix arrays for final construction
"""
function build_U_matrix(
    dofs::SystemDofs,
    ops::AbstractVector{<:Operators},
    blocks=nothing
)
    isempty(ops) && error("Cannot build U matrix from empty operators list")

    T = typeof(first(ops).value)
    N = total_dim(dofs)

    # Helper function to check if (k,l) are in the same block
    function in_same_block(k, l, blocks)
        blocks === nothing && return true
        for block in blocks
            if k in block && l in block
                return true
            end
        end
        return false
    end

    # Use dictionary to accumulate U matrix elements
    # Key: (row, col), Value: matrix element
    U_dict = Dict{Tuple{Int, Int}, T}()

    for op in ops
        length(op.ops) == 4 || error("Two-body operator must have exactly 4 operators")
        all(o isa FermionOp for o in op.ops) || error("U matrix requires FermionOp")

        # Reorder to InterAll format: c†_i c_j c†_k c_l
        fermion_ops = Vector{FermionOp}(op.ops)
        sign, reordered = _reorder_to_interall(fermion_ops)

        i = qn2linear(dofs, reordered[1].qn)
        j = qn2linear(dofs, reordered[2].qn)
        k = qn2linear(dofs, reordered[3].qn)
        l = qn2linear(dofs, reordered[4].qn)

        v_val = sign * op.value

        # Apply 4-term formula: U_{ijkl} = V_{ijkl} + V_{klij} - V_{kjil} - V_{ilkj}
        # Each V_{ijkl} contributes to 4 different U elements:

        # 1. U_{ijkl} += V_{ijkl}
        row1 = (i - 1) * N + j
        col1 = (k - 1) * N + l
        if in_same_block(k, l, blocks)
            U_dict[(row1, col1)] = get(U_dict, (row1, col1), zero(T)) + v_val
        end

        # 2. U_{klij} += V_{ijkl} (this V_{ijkl} appears as V_{klij} in U_{klij}'s formula)
        row2 = (k - 1) * N + l
        col2 = (i - 1) * N + j
        if in_same_block(i, j, blocks)
            U_dict[(row2, col2)] = get(U_dict, (row2, col2), zero(T)) + v_val
        end

        # 3. U_{kjil} -= V_{ijkl} (this V_{ijkl} appears as V_{kjil} in U_{kjil}'s formula)
        row3 = (k - 1) * N + j
        col3 = (i - 1) * N + l
        if in_same_block(i, l, blocks)
            U_dict[(row3, col3)] = get(U_dict, (row3, col3), zero(T)) - v_val
        end

        # 4. U_{ilkj} -= V_{ijkl} (this V_{ijkl} appears as V_{ilkj} in U_{ilkj}'s formula)
        row4 = (i - 1) * N + l
        col4 = (k - 1) * N + j
        if in_same_block(k, j, blocks)
            U_dict[(row4, col4)] = get(U_dict, (row4, col4), zero(T)) - v_val
        end
    end

    # Convert dictionary to sparse matrix
    nnz_count = length(U_dict)
    rows = Vector{Int}(undef, nnz_count)
    cols = Vector{Int}(undef, nnz_count)
    vals = Vector{T}(undef, nnz_count)

    idx = 0
    for ((row, col), val) in U_dict
        if !iszero(val)
            idx += 1
            rows[idx] = row
            cols[idx] = col
            vals[idx] = val
        end
    end

    # Resize to actual non-zero count
    resize!(rows, idx)
    resize!(cols, idx)
    resize!(vals, idx)

    return sparse(rows, cols, vals, N^2, N^2)
end
