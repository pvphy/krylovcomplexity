module Hamiltonian

export apply_xxz!

@inline function sz(state::Int, i::Int)
    ((state >> (i-1)) & 1) == 1 ? 0.5 : -0.5
end

@inline function flip(state::Int, i::Int, j::Int)
    bi = (state >> (i-1)) & 1
    bj = (state >> (j-1)) & 1
    bi == bj && return -1
    return state ⊻ ((1<<(i-1)) | (1<<(j-1)))
end

"""
Compute out = H * v for the XXZ model.
"""
function apply_xxz!(
    out::Vector{Float64},
    v::Vector{Float64},
    basis::Vector{Int},
    index::Dict{Int,Int},
    L::Int,
    J::Float64,
    Δ::Float64
)
    fill!(out, 0.0)

    for (α,state) in enumerate(basis)
        amp = v[α]
        amp == 0 && continue

        for i in 1:L-1
            #diagonal
            out[α] += J * Δ * sz(state,i) * sz(state,i+1) * amp

            #off-diagonal
            new = flip(state,i,i+1)
            if new ≥ 0
                β = index[new]
                out[β] += 0.5 * J * amp
            end
        end
    end
end

end 
