include("../src/Opt.jl")
include("../src/Model.jl")
h = Model.hill(k=0.5,n=3)

@show Opt.find_min(h,[0.,1.])

@show Opt.find_max(hill,[0.,30.])

@show Opt.find_MinMax(hill, [0.,30.])

using Optim

xx(w) = w

function t(x;functor)

return functor(x)

Optim.optimize(h, 0.1, 0.2)


##
using Optim
using Test

struct y
    p
end

function (s::y)(x)
    return x*s.p
end

model = y(1.0)
model2(x;p=1.0) = x*p

@show  model(2.0)
@show  model2(2.0)

@test typeof(y) <: Any # Because Optim.jl accepts :Any object

Optim.optimize(model, 0.1, 0.2);
Optim.optimize(model2, 0.1, 0.2);
