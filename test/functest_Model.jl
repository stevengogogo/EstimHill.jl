include("../src/Model.jl")
using Plots

using ModelingToolkit
using Calculus
using ForwardDiff
using OrdinaryDiffEq
using DifferentialEquations
##


## Use Functor

n = 10000
p1 = LinRange(0.001, 10, n)
out = zeros(n)

@time for (c, i) in enumerate(p1)
    R= Model.positive_feedback(s=i)
    prob = SteadyStateProblem(R, 0.0)

    #prob = ODEProblem(R, R.x, (0.0,100.0))
    local sol = solve(prob, DynamicSS(Rosenbrock23()))
    out[c] = sol.u
end

plot(p1, out)
plot!(xaxis=:log)


##

@parameters t s kK kA m x_50 x_tot kp
@variables x(t)
@derivatives D'~t

s = Model.positive_feedback(x,s, kK, kA, m, x_50, x_tot, kp)

dx = s(x, nothing, nothing)
