#' ---
#' title: Sensitivity analysis on Hill Model
#' author: Shao-Ting Chiu (r07945001@ntu.edu.tw)
#' date: 2020/08/27
#' ---

#' # Introduction
#' We are going to explore the

using PyCall
using Plots
using Parameters

include("../src/Sens.jl")
include("../src/Model.jl")
sp = pyimport("sympy")

#' Model: Hill Formula



s = sp.symbols("s", nonnegative=true) # Output / input
k, n = sp.symbols("k n", nonnegative=true) # dissociation coefficient constant and Hill coefficient

model = Model.hill(k=k, n=n)
X = model(s)

# Local Sensitivity
Sens.localSens(X, s)

rs = Sens.responseCoef(X, s)



frac = Sens.actFrac(X, X.subs(s,0), X.subs(s,1) )


rs
sp.solve(rs-n,frac)
