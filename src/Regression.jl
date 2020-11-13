module Regression

using Parameters
import LsqFit
import Dierckx
const lq = LsqFit
const dk = Dierckx

include("Utils.jl")
include("Opt.jl")
using .Opt
using .Utils

export curve_fit, is_monotonous


const DEF_RESOLUTION = 0.01

#=
function curve_fit(model, xdata, ydata, p0, lb, ub)
    fit = lq.curve_fit(model, xdata, ydata, p0, lower=lb, upper=ub)

    if !fit.converged
        @warn "The fitting is not converged with residual(=$(param.resid[end]))"
    end

    param = fit.param

    return param

end

function curve_fit(model, xdata, ydata, p0)
    fit = lq.curve_fit(model, xdata, ydata, p0)

    if !fit.converged
        @warn "The fitting is not converged with residual(=$(param.resid[end]))"
    end

    param = fit.param

    return param

end
=#


"""
Fit model with given x,y data.

# Argument
- model::`function`: function with input and parameters. model(x,p)
- xdata::`Array`:input data
- ydata:: `Array`: output data
- p0::'tuple{1,n}': n = number of inputs
- lb::`Array{1,n}`: lower bound of n inputs
- ub::`Array{1,n}`: upper bound of n inputs
- method::`Symbol`: Method of LsqFit. it can be `:forwarddiff` or `:finiteforward`
"""
function curve_fit(model_functor, xdata, ydata, p0; method=:forwarddiff)
    model_func =Utils.convert_functor_model(model_functor)

    fit =  lq.curve_fit(model_func, xdata, ydata, p0;autodiff=method)


    if !fit.converged
        @warn "The fitting is not converged"
    end

    return fit
end

function curve_fit(model_functor, xdata, ydata, p0, lb, ub; method=:forwarddiff)
    model_func = Utils.convert_functor_model(model_functor)
    fit = lq.curve_fit(model_func, xdata, ydata, p0, lower=lb, upper=ub; autodiff=method)

    if !fit.converged
        @warn "The fitting is not converged"
    end

    return fit
end


"""
Check the monotony of a series.
"""
function is_monotonous(series)
    is_mono_increase = all(series == accumulate(max,series) )
    is_mono_decrease = all(series == accumulate(min,series) )
    is_mono = is_mono_increase | is_mono_decrease
    return is_mono
end


"""
Check the mononous of a function with fixed domain.

# Argument
- f::`function`: the function to be tested for monotony.
- domain::`array of 2 numbers`: domain of function.
- resolution::`float`: resolution of screening
- reset_step:: `float`: reset resolution if resolution is larger than the domain range
"""
function is_monotonous(f, domain; resolution=0.01, reset_step=100.0)


    data = Utils.screen_domain(f, domain; resolution=0.01, reset_step=100.0)

    is_mono = is_monotonous(data.y)

    return is_mono
end


"""
Measure the Effective concentration with given ratio

# Arguement
- arrayX{Vector}: should be monotonous
- arrayY{Vector}: should be monotonous
- EC_ratio{Number}: The effective concentration
"""
function measure_EC( arrayX , arrayY , EC_ratio)

    if !is_monotonous(arrayX) | !is_monotonous(arrayY)
        @warn "The given curve is not monotonous"
    end

    if (EC_ratio < 0.0 )| (EC_ratio >= 1.0)
        @warn "Invalid EC_ratio($EC_ratio). EC_ratio ∈ [0,1)"
    end

    min_v, max_v = minimum(arrayY), maximum(arrayY)

    threshold_v = (max_v-min_v) * EC_ratio

    spl = dk.Spline1D(arrayX, arrayY .- threshold_v .- min_v)

    threshold_inputs = dk.roots(spl)

    return threshold_inputs
end

"""Measure Effective concentration with function input and screening resolution."""
function measure_EC(f, domain, EC_ratio, resolution)
    data = Utils.screen_domain(f, domain; resolution=resolution)

    threshold_inputs = measure_EC(data.x, data.y, EC_ratio)

    return threshold_inputs
end

end
