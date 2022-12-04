using FMMProject
using SafeTestsets

@safetestset "Utils" begin include("utils.jl") end

@safetestset "FMM" begin include("fmm.jl") end
