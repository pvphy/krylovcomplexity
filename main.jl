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

function disorder(rng, L::Int, W::Float64)
    return rand(rng,L) .* W .- W/2
end

function neel_state(L)
    state1 = 0
    for i in 1:L
        if isodd(i)
            state1 |= (1 << (L - i))
        end
    end
    return state1
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

vals=read_input_file("input.dat")

L= Int(vals[1])
J=vals[2]          # Jx,Jy
delta =vals[3]          # Jz
m=Int(vals[4])
W=vals[5]            #disoredr stregth
seed=Int(vals[6])
init_flag =Int(vals[7])
Nup=L÷2              #no of up spins:Ndn=L-Nup

BLAS.set_num_threads(4)

println("L                 =",L)
println("Jx,Jy             =",J)
println("Jz                =",delta)
println("Lanczos vectors m =",m)
println(" W                =",W)
println("Seed              =",seed)
println("Nup               =",Nup)
println("Sz                =",(Nup-(L-Nup))/2.0)
println("initial state     =",init_flag)


rng=MersenneTwister(seed)
h=disorder(rng,L,W)


basis=generate_basis(L,Nup)
index=state_index(basis)
dim=length(basis)

println("dimension = ",dim)

applyH!(out, v) = apply_xxz!(out,v,basis,index,L,J,delta,h)

if init_flag==194264
    eigvals,kry_ham=lanczos(applyH!,dim;m=m,rng=rng,init=:random)
    evolve_krylov(kry_ham;tmin=0.0,tmax=100.0,Nt=400,prefix="random",seed,L,J,delta)

elseif init_flag==121212
    psi00 = zeros(Float64,dim)
    ino=index[neel_state(L)]
    # println(ino)
    psi00[ino]=1.0

    eigvals,kry_ham=lanczos(applyH!,dim;m=m,rng=rng,init=:neel,v0=psi00)
    evolve_krylov(kry_ham;tmin=0.0,tmax=100.0,Nt=400,prefix="neel",seed,L,J,delta)
end

# println("Ground-state energy = ",minimum(eigvals))
# println(basis)                    
# state = basis[2]
# println(state_bits(682, L))








