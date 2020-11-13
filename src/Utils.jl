module Utils

using Parameters, DifferentialEquations

export data, convert_functor_model, data, screen_domain, rangeStepsize


@with_kw mutable struct data
    x
    y
end


@with_kw struct MOD
    mod
end
function (set::MOD)(x, p)
    return @. set.mod(p...)(x)
end

function convert_functor_model(model_functor)
    return MOD(model_functor)
end

"Range with fixed step size"
rangeStepsize(start , stop ,resolution) = LinRange(start, stop, Int(abs(stop-start) ÷ resolution))


function screen_domainSS(f, domain; resolution=0.01, reset_step=100.0)
    x_lower, x_upper = domain[1], domain[2]

    if resolution> abs(x_upper-x_lower)
        res_old = resolution
        resolution = abs(x_upper-x_lower) / reset_step
        @warn "The resolution ($res_old) larger than the domain range($domain). Reset resolution to $(resolution)"
    end

    x_data = rangeStepsize(x_lower, x_upper, resolution)
    y_data = f.(x_data)

    return data(x=x_data, y=y_data)
end


"""Screening Domain with given input parameter"""
function screen_domainSS(mutable_functor, input_par::Symbol, domain, u0; resolution=0.01, method=DynamicSS(Rosenbrock23()))
    xdata = Utils.rangeStepsize(domain[1], domain[2], resolution)
    ydata = []
    for x in xdata
        sol = get_ss(mutable_functor, u0, input_par, method)(x)
        append!(ydata, sol)
    end
    return data(x=xdata, y=ydata)
end

@with_kw struct get_ss
    functor_model
    u0
    input_par
    method= DynamicSS(AutoTsit5(Rosenbrock23()), abstol=1e-8) # Recommended solver for positive_feedback model
end

function (set::get_ss)(x)
    @unpack functor_model, u0, input_par, method = set
    setfield!(functor_model, input_par, x)
    prob = SteadyStateProblem(functor_model, u0)
    sol = solve(prob, method)
    return sol.u
end




end
