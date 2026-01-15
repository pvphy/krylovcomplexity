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

L   = 20
J   = 1.0
delta   = 1.0
Nup = L ÷ 2         
m   = 80             


basis = generate_basis(L, Nup)
index = state_index(basis)
dim   = length(basis)

println("dimension = ", dim)

applyH!(out, v) = apply_xxz!(out, v, basis, index, L, J, delta)

eigvals,kry_ham = lanczos(applyH!, dim; m=m)
println("Ground-state energy = ", minimum(eigvals))


evolve_krylov(kry_ham;tmin = 0.0,tmax = 400.0,Nt   = 400,prefix = "run1")
