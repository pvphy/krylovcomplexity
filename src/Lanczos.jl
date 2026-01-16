module Lanczos
    using LinearAlgebra
    export lanczos
    function lanczos(applyH!, dim::Int; m::Int=50,rng)
        v= randn(dim)
        v./= norm(v)

        w=zeros(dim)
        v_old=zeros(dim)

        alpha=zeros(m)
        beta=zeros(m-1)

        for j in 1:m
            applyH!(w,v)
            alpha[j] = dot(v,w)

            if j > 1
                w .-= beta[j-1] .*v_old
            end

            w .-= alpha[j] .*v

            if j < m
                beta[j]=norm(w)
                beta[j] < 1e-12 && break
                v_old .= v
                v .= w ./ beta[j]
            end
        end

        T = SymTridiagonal(alpha,beta)
        return eigen(T).values,T
    end

end 