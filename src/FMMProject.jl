module FMMProject

using LinearAlgebra
using StaticArrays

export FMMLaplace2D,
       FMMLaplace3D,
       Point2D,
       Point3D

include("utils.jl")
include("tree.jl")
include("fmm.jl")
include("fmm_operators.jl")
include("fmm_assembly.jl")
include("fmm_algorithm.jl")
include("fmm_singlelevel.jl")

end # module
