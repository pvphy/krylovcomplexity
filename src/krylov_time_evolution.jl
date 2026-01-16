module KrylovTimeEvolution
    using LinearAlgebra
    export evolve_krylov

    function evolve_krylov(kry_ham;tmin=0.0,tmax=50.0,Nt=50,prefix = "krylov")

        M=size(kry_ham, 1)


        phi0=zeros(ComplexF64, M)
        phi0[1]=1.0+0im


        eig=eigen(kry_ham)
        U=eig.vectors
        E=eig.values


        c0=U'*phi0


        times=range(tmin, tmax, length=Nt)


        nvals=collect(0:M-1)


        open("$(prefix)_probabilities.txt", "w") do io_p
        open("$(prefix)_complexity.txt", "w") do io_k


            for t in times

                phase=exp.(-1im .* E .* t)
                phi=U*(phase .* c0)
                pn=abs2.(phi)

                for n in 1:M
                    println(io_p, t, "  ", nvals[n], "  ", pn[n])
                end
                println(io_p) 
            
                Kt=sum(nvals .* pn)
                println(io_k, t, "  ", Kt)
            end

        end
        end

        return nothing
    end

end 
