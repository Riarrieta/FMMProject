
struct FMM{N,P,K<:AbstractKernel{N}}
    tree::TreeNode{N}
    levels::Vector{Vector{TreeNode{N}}}
    points::Vector{Point{N}}
    result::Vector{ComplexF64}
    maxpointsperleaf::Integer
end

const FMMLaplace{N,P} = FMM{N,P,<:Laplace{N}} where {N,P}
const FMMLaplace2D{P} = FMM{2,P,<:Laplace2D} where P
const FMMLaplace3D{P} = FMM{3,P,<:Laplace3D} where P

kernel(::FMM{N,P,K}) where {N,P,K} = K
levels(f::FMM) = f.levels
nlevels(f::FMM) = length(levels(f))
points(f::FMM) = f.points 
maxpointsperleaf(f::FMM) = f.maxpointsperleaf
eachlevel(f::FMM) = ((l-1,levels(f)[l]) for l in 1:nlevels(f))