module Hamiltonian

export apply_xxz!

@inline function sz(state::Int, i::Int)
    ((state >> (i-1)) & 1) == 1 ? 0.5 : -0.5
end

@inline function flip(state::Int, i::Int, j::Int)
    bi=(state >> (i-1)) & 1
    bj=(state >> (j-1)) & 1
    bi== bj && return -1
    return state ⊻ ((1<<(i-1)) | (1<<(j-1)))
end


function apply_xxz!(out::Vector{Float64},v::Vector{Float64},basis::Vector{Int},index::Dict{Int,Int},
    L::Int,J::Float64,delta::Float64,h)
    fill!(out, 0.0)

    for (ii,state) in enumerate(basis)
        amp=v[ii]
        amp==0 && continue

        for i in 1:L-1
            #diagonal
            out[ii] += J * delta*sz(state,i)*sz(state,i+1)*amp

            #off-diagonal
            new = flip(state,i,i+1)
            if new ≥ 0
                β = index[new]
                out[β] += 0.5*J*amp
            end
        end


        # onsite disorder
        for i in 1:L
            out[ii] += h[i]*sz(state,i)*amp
        end
    end
end

end 
