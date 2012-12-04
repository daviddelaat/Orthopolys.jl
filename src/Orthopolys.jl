module Orthopolys

export three_term_recurrence,
       jacobi_coefmatrix,
       gegenbauer_coefmatrix,
       prob_hermite_coefmatrix,
       phys_hermite_coefmatrix
       
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
     
function jacobi_coefmatrix{T<:Number}(alpha::T, beta::T, n::Int)
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

function jacobi_coefmatrix{T<:Number}(alpha::T, beta::T, n::Int, normalization::Symbol)
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

function gegenbauer_coefmatrix{T<:Number}(alpha::T, n::Int)
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

function gegenbauer_coefmatrix{T<:Number}(alpha::T, n::Int, normalization::Symbol)
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

function prob_hermite_coefmatrix(T, n::Int)
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
        m[l, 1] = - (l - 1) * m[l-2, 1]
        for i in 2:n
            m[l, i] = m[l-1, i-1] - (l - 1) * m[l-2, i])
        end
    end

    m 
end

function phys_hermite_coefmatrix(T, n::Int)
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
        m[l, 1] = - 2 * (l - 1) * m[l-2, 1]
        for i in 2:n
            m[l, i] = 2 * (m[l-1, i-1] - (l - 1) * m[l-2, i])
        end
    end

    m 
end

end
