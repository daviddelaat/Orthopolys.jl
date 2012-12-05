module Orthopolys

export three_term_recurrence,
       horner,
       jacobi_coefmatrix,
       jacobi_eval,
       gegenbauer_coefmatrix,
       gegenbauer_eval,
       prob_hermite_coefmatrix,
       prob_hermite_eval,
       phys_hermite_coefmatrix,
       phys_hermite_eval      
       
# Finds the n x n coefficient matrix of the polynomial defined by the recurrence relation:
# p_0(z) = p
# p_1(z) = q + r z
# p_k(z) = a(k) p_{k-1}(z) + b(k) p_{k-2}(z) + c(k) z p_{k-1)(z)
function three_term_recurrence{T<:Number}(p::T, q::T, r::T, a::Function, b::Function, c::Function, n::Int)
    m = zeros(T, n, n)
    
    if n == 0
        return m
    end
    
    m[1, 1] = p
    
    if n == 1
        return m
    end
    
    m[2, 1] = q
    m[2, 2] = r
   
    for l in 3:n
        m[l, 1] = b(l-1) * m[l-2, 1]
        for i in 2:n
            m[l, i] = a(l-1) * m[l-1, i] + b(l-1) * m[l-2, i] + c(l-1) * m[l-1, i-1]
        end
    end

    m
end

function horner{T<:Number}(coefficients::Vector{T}, x::T)
    r = zero(T)
    for coefficient in coefficients[end:-1:1]
        r = r * x + coefficient
    end
    r
end

function jacobi_coefmatrix{T<:Number}(n::Int, alpha::T, beta::T)
    m = zeros(T, n, n)
    
    if n == 0
        return m
    end
       
    m[1, 1] = one(T)
    
    if n == 1
        return m
    end
    
    m[2, 1] = (alpha - beta) / 2
    m[2, 2] = 1 + (alpha + beta) / 2
   
    t = alpha + beta + 4
   
    for l in 3:n
        a = 2 * (l - 1) * (t - 2) * (l + alpha + beta - 1)
        b = (t - 2) * (t - 1) * t
        c = (t - 1) * (alpha^2 - beta^2)
        d = -2t * (l + alpha - 2) * (l + beta - 2)
        
        m[l, 1] = (c * m[l-1, 1] + d * m[l-2, 1]) / a
        for i in 2:n
            m[l, i] = (c * m[l-1, i] + d * m[l-2, i] + b * m[l-1, i-1]) / a
        end
        
        t += 2
    end

    m
end

function jacobi_coefmatrix{T<:Number}(n::Int, alpha::T, beta::T, normalization::Symbol)
    m = jacobi_coefmatrix(alpha, beta, n)
    if normalization == :std
        return m
    elseif normalization == :oneatone
        for i = 1:n
            m[i, :] /= sum(m[i, :])
        end
        return m
    end
end

jacobi_eval{T<:Number}(n::Int, alpha::T, beta::T, x::T) =
    horner(vec(jacobi_coefmatrix(n, alpha, beta)[n, :]), x)

jacobi_eval{T<:Number}(n::Int, alpha::T, beta::T, x::T, normalization::Symbol) =
    horner(vec(jacobi_coefmatrix(n, alpha, beta, normalization)[n, :]), x)

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
        m[l, 1] = - (l + 2 * alpha - 3) * m[l-2, 1] / (l - 1)
        for i in 2:n
            m[l, i] = (2 * (l + alpha - 2) * m[l-1, i-1] - (l + 2alpha - 3) * m[l-2, i]) / (l - 1)
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
            m[i, :] /= sum(m[i, :])
        end
        return m
    end
end

gegenbauer_eval{T<:Number}(n::Int, alpha::T, x::T) =
    horner(vec(gegenbauer_coefmatrix(n, alpha)[n, :]), x)

gegenbauer_eval{T<:Number}(n::Int, alpha::T, x::T, normalization::Symbol) =
    horner(vec(gegenbauer_coefmatrix(n, alpha, normalization)[n, :]), x)

function prob_hermite_coefmatrix(n::Int, T)
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
        m[l, 1] = - (l - 2) * m[l-2, 1]
        for i in 2:n
            m[l, i] = m[l-1, i-1] - (l - 2) * m[l-2, i]
        end
    end

    m 
end

prob_hermite_eval{T<:Number}(n::Int, x::T) =
    horner(vec(prob_hermite_coefmatrix(n)[n, :]), x)

function phys_hermite_coefmatrix(n::Int, T)
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
        m[l, 1] = - 2 * (l - 2) * m[l-2, 1]
        for i in 2:n
            m[l, i] = 2 * (m[l-1, i-1] - (l - 2) * m[l-2, i])
        end
    end

    m 
end

phys_hermite_eval{T<:Number}(n::Int, x::T) =
    horner(vec(phys_hermite_coefmatrix(n)[n, :]), x)

end