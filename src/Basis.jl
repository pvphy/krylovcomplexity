module Basis

export generate_basis, state_index ,state_bits


function generate_basis(L::Int, Nup::Int)
    basis=Int[]
    for state in 0:(1<<L)-1
        count_ones(state) == Nup && push!(basis, state)
    end
    return basis
end

# state to index

function state_index(basis::Vector{Int})
    idx=Dict{Int,Int}()
    for (i,s) in enumerate(basis)
        idx[s]=i
    end
    return idx
end

# state to bits 
function state_bits(state::Int, L::Int)
    bits = zeros(Int, L)
    @inbounds for i in 1:L
        bits[L - i + 1]=(state >> (i - 1)) & 1
    end
    return bits
end

end # module
