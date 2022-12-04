
struct FMM{N,K<:AbstractKernel{N}}
    tree::TreeNode{N}
    levels::Vector{Vector{TreeNode{N}}}
    points::Vector{Point{N}}
    result::Vector{ComplexF64}
    maxpointsperleaf::Integer
    interactionrank::Integer
end

const FMMLaplace{N} = FMM{N,<:Laplace{N}} where {N}
const FMMLaplace2D = FMM{2,<:Laplace2D}
const FMMLaplace3D = FMM{3,<:Laplace3D}

kernel(::FMM{N,K}) where {N,K} = K
levels(f::FMM) = f.levels
nlevels(f::FMM) = length(levels(f))
points(f::FMM) = f.points 
maxpointsperleaf(f::FMM) = f.maxpointsperleaf
interaction_rank(f::FMM) = f.interaction_rank
eachlevel(f::FMM) = ((l-1,levels(f)[l]) for l in 1:nlevels(f))