"Derive Effective Concentration with Regression.jl"

include("../src/Regression.jl")
include("../src/Model.jl")
include("../src/Utils.jl")
include("../src/Opt.jl")
import PyPlot
using Printf
using Test
using DifferentialEquations
using Parameters

const plt = PyPlot

## Model and Domain
sample_model = Model.hill(k=0.5, n=4.)
domain = [0.0,1.0]
data = Utils.screen_domain(sample_model, domain)

## Get Thresholds
thr_50 = Regression.measure_EC(data.x, data.y, 0.5)[1]
thr_90 = Regression.measure_EC(data.x, data.y, 0.9)[1]
thr_10 = Regression.measure_EC(data.x, data.y, 0.1)[1]
@test Regression.measure_EC(sample_model, domain,0.5,0.01)[1] ≈ thr_50 # Simulation method and data method


### Get extremes
argmin = Regression.measure_EC(data.x, data.y, 0.)[1]
argmax = Regression.measure_EC(data.x, data.y, 1-1e-10)[1]





### Visualization



fig, ax = plt.subplots(figsize=(10,7))
color = ["black", "green", "red", "blue", "gray"]
ax.plot(data.x, data.y, label="data", color=color[1])
ax.axvline(thr_50, label= string("EC50 = $(@sprintf("%.2f", thr_50))"), color=color[2])
ax.axvline(thr_90, label="EC90 = $(@sprintf("%.2f", thr_90))", color=color[3])
ax.axvline(thr_10, label="EC10 = $(@sprintf("%.2f", thr_10))", color=color[4])

ax.axvline(argmin, label="argMin = $(@sprintf("%.2f", min))", color=color[5])
ax.axvline(argmax, label="argMax = $(@sprintf("%.2f", max))", color=color[5])

ax.legend()
ax.set_title("Estimation of Effective Concentration (EC)")
ax.set_xlabel("Input")
ax.set_ylabel("Output (Monotonous)")
plt.display_figs()


## Positive Feedback

input = LinRange(0.1,1.0,1000)
output = zeros(length(input))



for (i, s) in enumerate(input)
    local mod = Model.positive_feedback(s=s)

    local prob = SteadyStateProblem(mod, mod.x)

    local sol = solve(prob, DynamicSS(Rosenbrock23()))

    output[i] = sol.u[1]
end

thr_50 = Regression.measure_EC(input, output, 0.5)[1]
thr_90 = Regression.measure_EC(input, output, 0.9)[1]
thr_10 = Regression.measure_EC(input, output, 0.1)[1]


### Get extremes
argmin = Regression.measure_EC(input, output, 0.)[1]
argmax = Regression.measure_EC(input, output, 1-1e-10)[1]

fig, ax = plt.subplots(figsize=(10,7))
color = ["black", "green", "red", "blue", "gray"]
ax.plot(input, output, label="data", color=color[1])
ax.axvline(thr_50, label= string("EC50 = $(@sprintf("%.2f", thr_50))"), color=color[2])
ax.axvline(thr_90, label="EC90 = $(@sprintf("%.2f", thr_90))", color=color[3])
ax.axvline(thr_10, label="EC10 = $(@sprintf("%.2f", thr_10))", color=color[4])

ax.axvline(argmin, label="argMin = $(@sprintf("%.2f", min))", color=color[5])
ax.axvline(argmax, label="argMax = $(@sprintf("%.2f", max))", color=color[5])

ax.set_xscale("log")
ax.legend()
ax.set_title("Estimation of Effective Concentration (EC)")
ax.set_xlabel("Input")
ax.set_ylabel("Output (Monotonous)")
plt.display_figs()

##

#=
Utils.screen_domain(Model.positive_feedback, 0.1, :s, [0.0,1.0])


df = Model.positive_feedback()

a = :s
df.a
df.:s

eval(quote df.$a=4 end)
=#
