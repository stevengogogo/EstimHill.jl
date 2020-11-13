module Model

using Parameters

"""Hill function without basal activation"""
@with_kw struct hill
    k
    n
end

function (par::hill)(s)
    @unpack k,n = par
    return s^n / (k^n + s^n )
end

"""Hill function with basal activation"""
@with_kw struct hill_basal
    k
    n
    x_basal
    x_max
end
function (par::hill_basal)(s)
    @unpack k, n, x_base, x_max = par

    return x_base + (x_max - x_base)*(s^n/(k^n+s^n))
end

"""Model: Positive Feedback"""
@with_kw mutable struct positive_feedback
    x= 0.0
    s=0.1
    kK=1.
    kA=1.3
    m=3.
    x_50=0.5
    x_tot=1.0
    kp=1.0

    #@assert  sum(x) <= x_tot
end

function (p::positive_feedback)(x,p_,t)
    @unpack s, kK, kA, m, x_50, x_tot, kp = p
    dx = (kK * s + kA * (x^m/(x_50^m + x^m))) * (x_tot-x) - kp * x

    return dx
end


end


#=

@parameters t s kK kA m x_50 x_tot kp
@variables x(t)
@derivatives D'~t
=#
