module Orthopolys

using Polynomial

export
    jacobi,
    gegenbauer,
    prhermite,
    phhermite

function jacobi_coefmatrix{T<:Number}(n::Int, alpha::T, beta::T)
    m = zeros(T, n, n)
    
    if n == 0
        return m
    end
       
    m[1, 1] = one(T)
    
    if n == 1
        return m
    end
    
    m[1, 2] = (alpha - beta) / 2
    m[2, 2] = 1 + (alpha + beta) / 2
   
    t = alpha + beta + 4
   
    for l in 3:n
        a = 2 * (l - 1) * (t - 2) * (l + alpha + beta - 1)
        b = (t - 2) * (t - 1) * t
        c = (t - 1) * (alpha^2 - beta^2)
        d = -2t * (l + alpha - 2) * (l + beta - 2)
        
        m[1, l] = (c * m[1, l-1] + d * m[1, l-2]) / a
        for i in 2:n
            m[i, l] = (c * m[i, l-1] + d * m[i, l-2] + b * m[i-1, l-1]) / a
        end
        
        t += 2
    end

    m
end

function jacobi_coefmatrix{T<:Number}(n::Int, alpha::T, beta::T, normalization::Symbol)
    m = jacobi_coefmatrix(n, alpha, beta)
    if normalization == :std
        return m
    elseif normalization == :oneatone
        for i = 1:n
            m[:, i] /= sum(m[:, i])
        end
        return m
    else
        error("jacobi does not support $normalization normalization")
    end
end

jacobi{T<:Number}(n::Int, alpha::T, beta::T) =
    Poly(jacobi_coefmatrix(n, alpha, beta)[:, n])

jacobi{T<:Number}(n::Int, alpha::T, beta::T, normalization::Symbol) =
    Poly(jacobi_coefmatrix(n, alpha, beta, normalization)[:, n])

function gegenbauer_coefmatrix{T<:Number}(n::Int, alpha::T)
    m = zeros(T, n, n)
    
    if n == 0
        return m
    end
    
    m[1, 1] = one(T)
    
    if n == 1
        return m
    end
    
    m[2, 2] = 2 * alpha
   
    for l in 3:n
        m[1, l] = - (l + 2 * alpha - 3) * m[1, l-2] / (l - 1)
        for i in 2:n
            m[i, l] = (2 * (l + alpha - 2) * m[i-1, l-1] - (l + 2alpha - 3) * m[i, l-2]) / (l - 1)
        end
    end

    m
end

function gegenbauer_coefmatrix{T<:Number}(n::Int, alpha::T, normalization::Symbol)
    m = gegenbauer_coefmatrix(alpha, n)
    if normalization == :std
        return m
    elseif normalization == :oneatone
        for i = 1:n
            m[:, i] /= sum(m[:, i])
        end
        return m
    else
        error("gegenbauer does not support $normalization normalization")
    end
end

gegenbauer{T<:Number}(n::Int, alpha::T) =
    Poly(gegenbauer_coefmatrix(n, alpha)[:, n])

gegenbauer{T<:Number}(n::Int, alpha::T, normalization::Symbol) =
    Poly(gegenbauer_coefmatrix(n, alpha, normalization)[:, n])

function prhermite_coefmatrix(n::Int, T)
    m = zeros(T, n, n)
    
    if n == 0
        return m
    end
    
    m[1, 1] = one(T)
    
    if n == 1
        return m
    end
    
    m[2, 2] = one(T)
      
    for l in 3:n
        m[1, l] = - (l - 2) * m[1, l-2]
        for i in 2:n
            m[i, l] = m[i-1, l-1] - (l - 2) * m[i, l-2]
        end
    end

    m
end

prhermite_coefmatrix(n::Int) = prhermite_coefmatrix(n, Float64)

prhermite(n::Int, T) = Poly(prhermite_coefmatrix(n, T)[:, n])

function phhermite_coefmatrix(n::Int, T)
    m = zeros(T, n, n)
    
    if n == 0
        return m
    end
    
    m[1, 1] = one(T)
    
    if n == 1
        return m
    end
    
    m[2, 2] = 2 * one(T)
      
    for l in 3:n
        m[1, l] = - 2 * (l - 2) * m[1, l-2]
        for i in 2:n
            m[l, i] = 2 * (m[i-1, l-1] - (l - 2) * m[i, l-2])
        end
    end

    m
end

phhermite_coefmatrix(n::Int) = phhermite_coefmatrix(n, Float64)

phhermite{T<:Number}(n::Int, T) = Poly(phhermite_coefmatrix(n, T)[:, n])

end