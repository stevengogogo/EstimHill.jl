@warn "Opt.jl is not compatible with functors. See https://discourse.julialang.org/t/functor-usage-for-optim-jl/45778"

"
Optimization module for finding minimum and maximum of given function with single input (f(x))
"
module Opt

using Optim
using Parameters

export find_min, find_max, find_MinMax



@with_kw struct optRes
    f_opt
    x_arg
    summary
end

@with_kw struct minmax
    f_min
    f_max
    x_amin
    x_amax
    summary
end

"
Find minimimum of function within constraint domain

# Argument
- f::`function`: function with single input. That is f(x) = y
- domain::`vector{Number, 2}`: searching domain [lower, upper].

# Return
- optRes
"
function find_min(f, domain; method=BFGS())
    x_lower, x_upper = domain[1], domain[2]
    res = Optim.optimize(f, x_lower, x_upper, method)
    f_min = res.minimum
    x_amin = res.minimizer

    # Warning if the result is not converged
    if !res.converged
        @warn "Minimization failed. See summary\n$(res)"
    end

    return optRes(f_opt=f_min, x_arg=x_amin, summary=res)
end


"
Find maximum of function within constraint domain

# Argument
- f::`function`: function with single input. That is f(x) = y
- domain::`vector{Number, 2}`: searching domain [lower, upper].
"
function find_max(f, domain)
    inv_f(x; f=f) = f(x)*-1.0
    opres_min = find_min(inv_f, domain)
    opres_max = optRes(opres_min,
                       f_opt = -1.0 *opres_min.f_opt)
    return opres_max
end


function find_MinMax(f, domain)
    res_min = find_min(f, domain)
    res_max = find_max(f, domain)

    res = minmax(
            f_min=res_min.f_opt,
            f_max=res_max.f_opt,
            x_amin=res_min.x_arg,
            x_amax=res_max.x_arg,
            summary = @with_kw (min=res_min.summary, max=res_max.summary)
    )

    return res
end

end
