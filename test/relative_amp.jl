"Apply Relative Ampification Approach"

include("../src/Model.jl")
include("../src/Utils.jl")
include("../src/Regression.jl")
include("../src/Sens.jl")
include("../src/EstimHill.jl")
using DifferentialEquations
using Parameters
using PyPlot
const plt = PyPlot
## Get x y dose response
mod = Model.positive_feedback()
data = Utils.screen_domainSS(mod, :s, [0.1,10.0], mod.x_50/2.0; resolution=0.001)
#data = Utils.screen_domainSS(mod, [0.,10.0]; resolution=0.001)
##
x = data.y
s = data.x
x_basal =minimum(x)
x_max = maximum(x)


fraction_series = Sens.actFrac.(x, x_basal, x_max)


mod_ = Utils.get_ss(functor_model=mod,
                    u0=0.5,
                    input_par=:s,
                    method=DynamicSS(AutoTsit5(Rosenbrock23()), abstol=1e-8))

@time dRdS_CAL = Sens.responseCoef(mod_, s )



fig, ax = plt.subplots(figsize=(7,7))
ax.plot(fraction_series, dRdS_CAL,label="Calculus.jl")
ax.set_xlabel("Activated Fraction")
ax.set_ylabel("Response Coefficient")



### Plotting
fig2,ax2 = plt.subplots(figsize=(7,7))
ax2.plot(s_data, mod_.(s_data))


plt.display_figs()


## Apply relative _amplifiation approach

modF = Model.positive_feedback()
modF_ = Utils.get_ss(functor_model=modF,
                    u0=0.5,
                    input_par=:s,
                    method=DynamicSS(AutoTsit5(Rosenbrock23()), abstol=1e-8))

modc = Model.hill(0.5,1)

x = Utils.rangeStepsize(0.01,10.0, 0.001) # domain of input

@time sensRes = Sens.relative_amplification_approach(modF_, x)
@time sensRes_c = Sens.relative_amplification_approach(modc, x)

actf, respcoef= sensRes.activated_fraction, sensRes.f_response_coefficient
actf_c, respcoef_c = sensRes_c.activated_fraction, sensRes_c.f_response_coefficient


### Plotting
fig3, ax3 = plt.subplots()

ax3.plot( actf, respcoef, label="Positive Feedback Model")
ax3.plot(  actf_c, respcoef_c, label="Hill(k=$(modc.k),n=$(modc.n))")

plt.legend()
plt.display_figs()


##
