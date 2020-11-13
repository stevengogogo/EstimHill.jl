include("../src/Regression.jl")
include("../src/Model.jl")

using BenchmarkTools
using Parameters
using Test
using Plots
using LsqFit

## test functions
test_f(x) = x


hill = Model.hill


## Test LinRange with `step size`
for i in Regression.rangeStepsize(0.2, -13.34, 0.1)
    @show i
end

## Test the time eclapse of `rangeStepsize`
# Result: same allocation
@time test_f.(Regression.rangeStepsize(0, -1300, 1))
@time test_f.(LinRange(0,-1300,1300))



## Test monotonous testing
x_inc = [1,2,3,4,5]
x_dec = [5,3,2,1,0.4]
x_non = [0,1,2,0,2]

@test Regression.is_monotonous(test_f, [4,0]; resolution=20) == true
@test Regression.is_monotonous(x_inc) == true
@test Regression.is_monotonous(x_dec) == true
@test Regression.is_monotonous(x_non) == false



## Compare functor type

# test model
@with_kw struct test_m
    p1
    p2
end
function (set::test_m)(x)
    @unpack p1, p2 = set
    return p1*exp(-x*p2)
end



@. exp_model(x, p) = p[1]*exp(-x*p[2])


# model input
f_ob = Regression.convert_functor_model(test_m) # convert functor to function
@test exp_model(3, [1,2]) == f_ob(3, [1,2])


# Create data
xdata = range(0, stop=10, length=20)
ydata = exp_model(xdata, [1.0 2.0]) + 0.01*randn(length(xdata))
p0_bounds = [1.2, 1.2] # we have to start inside the bounds

# Three methods of curve fit
@show fit_functor1 = Regression.curve_fit(f_ob, xdata, ydata,  p0_bounds)
@show fit_tutorial =LsqFit.curve_fit(exp_model, xdata, ydata,  p0_bounds)
@show fit_functor = Regression.curve_fit(test_m, xdata, ydata,  p0_bounds)

@test fit_functor1.param ≈ fit_tutorial.param ≈ fit_functor.param

plot(xdata,ydata, label="Real")
plot!(xdata, test_m(fit_functor1.param...).(xdata), label="Fitted with Regression.jl")
plot!(xdata, exp_model(xdata, fit_tutorial.param), label="Fitted with LsqFit.jl")

## Fit Hill Model
#` Don't use [0] to fit hill function
k, n = 0.5, 6.0
model_real = hill(k, n )

@. H(x, p ) = x^p[2] /(p[1]^p[2] + x^p[2])


xdata = range(0, stop=1, length=20)
ydata = model.(xdata) #+ 0.0001 * rand(length(xdata))

fit = Regression.curve_fit(Model.hill, xdata, ydata,[0.1,0.1])
fit_err = Regression.curve_fit(Model.hill, xdata, ydata,[0.,0.])
fit1 = LsqFit.curve_fit(H, xdata, ydata, [1.,1.])
fit2 = LsqFit.curve_fit(H, xdata, ydata, [1.,1.]; autodiff=:forwarddiff)
fit3 = LsqFit.curve_fit(H, xdata, ydata, [1.,1.]; autodiff=:finiteforward)

#` Plotting
plot(xdata, ydata, label="Noise Data")
plot!(xdata, Model.hill(fit.param...).(xdata) ,label="Fitted Curve")
plot!(xdata, Model.hill(fit_err.param...).(xdata) ,label="Error Curve with 0 init")
plot!(xdata, H(xdata, fit1.param) ,label="Fitted Curve:Default")
plot!(xdata, H(xdata, fit2.param) ,label="Fitted Curve:forwarddiff")
plot!(xdata, H(xdata, fit3.param) ,label="Fitted Curve:finiteforward")

@test k ≈ fit.param[1]
@test n ≈ fit.param[2]
