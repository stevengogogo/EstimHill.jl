"
Measure Local sensitivity, Response coefficient and define fraction of activation
"

module Sens

using Calculus
using PyCall
using SymPy
using Parameters

sp = pyimport("sympy")

include("Regression.jl")

export localSens, responseCoef, actFrac

@with_kw struct sensitivity_result
    f_x
    f_dx
    x
    f_response_coefficient
    activated_fraction
end


"
Local Sensitivity Analysis

Compute the localsensitivity (df/ds) along the input s (number or vector).

# Arguments
- `f::function` : with single input f(s)
- `s`:: Single number (vector)
"
function localSens(f, s)
    df = derivative(f)
    return df(s)
end


function localSens(f, s)
    df = derivative(f)
    return df.(s)
end

function localSens(f, s::Sym)
    df = sp.diff(f,s)
    return sp.simplify(df)
end


"
Response Coefficient

Compute the response coefficient: \$R^{x}_{s} = \\frac{dx}{ds} \\frac{s}{x}\$.

where \$R^{x}_{s}\$ is the response coefficient of input \$s\$ and output \$x\$.

# Arguments
- `f::function` : with single input x = f(s)
- `s`:: Single number (can also be vector or sympy symbol).

# Warning

The result may be fluctuated due to the sinularity point.

# Reference:
1. Ferrell, J. E. & Ha, S. H. Ultrasensitivity part I: Michaelian responses and zero-order ultrasensitivity. Trends in Biochemical Sciences 39, 496–503 (2014).
2. Legewie, S., Blüthgen, N. & Herzel, H. Quantitative analysis of ultrasensitive responses. FEBS J. 272, 4071–4079 (2005).
"
function responseCoef(f, s)
    df = derivative(f)
    return df(s) * s / f(s)
end


function responseCoef(f, s)
    df = Calculus.derivative(f)
    return df.(s) .* s ./ f.(s)
end

function responseCoef(f,s::Sym)
    df = sp.diff(f,s)
    return sp.simplify(df * s / f)
end




"
Compute `activation fraction` of the stimulation process.

# Reference
1. Legewie, S., Blüthgen, N. & Herzel, H. Quantitative analysis of ultrasensitive responses. FEBS J. 272, 4071–4079 (2005).
"
actFrac(x, x_basal, x_max) = (x-x_basal) / (x_max - x_basal)


function relative_amplification_approach(f,x)
    df = Calculus.derivative(f) # Differential function
    f_x = f.(x) # f(x)
    dfdx = df.(x) # df/dx
    f_resp_coef = dfdx .* x ./ f_x # response coefficient

    y_max = maximum(f_x)
    y_min = minimum(f_x)

    act_fraction = actFrac.(f_x, y_min, y_max)

    return sensitivity_result(
        f_x=f_x,
        f_dx= dfdx,
        x=x,
        f_response_coefficient= f_resp_coef,
        activated_fraction=act_fraction
    )

end


end
