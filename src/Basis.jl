module Basis

export generate_basis, state_index

"""
    generate_basis(L, Nup)

Generate all Sz-basis states with Nup spins up.
States are represented as integers.
"""
function generate_basis(L::Int, Nup::Int)
    basis = Int[]
    for state in 0:(1<<L)-1
        count_ones(state) == Nup && push!(basis, state)
    end
    return basis
end

"""
    state_index(basis)

Return dictionary mapping state → index.
"""
function state_index(basis::Vector{Int})
    idx = Dict{Int,Int}()
    for (i,s) in enumerate(basis)
        idx[s] = i
    end
    return idx
end

end # module
