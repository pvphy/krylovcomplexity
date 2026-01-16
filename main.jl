####################################################################
#           computes krylov complexity
#               by   Prabhakar
#                  8/JAN/2026
####################################################################
using Random
# using MPI
using Printf
using Dates
using JLD2
using Statistics
using LinearAlgebra


include("src/Basis.jl")
include("src/Hamiltonian.jl")
include("src/Lanczos.jl")
include("src/krylov_time_evolution.jl")

using .Basis
using .Hamiltonian
using .Lanczos
using .KrylovTimeEvolution

function generate_disorder(rng, L::Int, W::Float64)
    return rand(rng,L) .* W .- W/2
end


function read_input_file(filename)
    values = Float64[]

    for line in eachline(filename)
        s=strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        val = split(s)[1]
        push!(values, parse(Float64, val))
    end

    return values
end


vals=read_input_file("input.txt")

L= Int(vals[1])
J=vals[2]          # Jx,Jy
delta =vals[3]          # Jz
m=Int(vals[4])
W=vals[5]            #disoredr stregth
seed=Int(vals[6])
Nup=L÷2              #no of up spins:Ndn=L-Nup



println("L                 = ",L)
println("Jx,Jy             = ",J)
println("delta (Jz)        = ",delta)
println("Lanczos vectors m = ",m)
println(" W                = ",W)
println("Seed              = ",seed)
println("Nup               = ",Nup)
println("Sz               = ",(Nup-(L-Nup))/2.0)

rng=MersenneTwister(seed)
h=generate_disorder(rng,L,W)


basis=generate_basis(L,Nup)
index=state_index(basis)
dim=length(basis)

println("dimension = ",dim)

applyH!(out, v) = apply_xxz!(out,v,basis,index,L,J,delta,h)

eigvals,kry_ham = lanczos(applyH!,dim;m=m,rng=rng)
println("Ground-state energy = ",minimum(eigvals))

                    
# state = basis[2]
# println(state_bits(state, L))

evolve_krylov(kry_ham;tmin=0.0,tmax=400.0,Nt=400,prefix = "run1")






