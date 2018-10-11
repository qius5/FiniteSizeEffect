module FiniteSizeEffect
Mzfc_vec(T, params) = Array{Float64}([Mzfc(t, params) for t in T])
using Example
using DelimitedFiles
using LinearAlgebra
include("naive_mf_cosine_fixI_Siwei.jl")
include("propagator_sub.jl")
end


