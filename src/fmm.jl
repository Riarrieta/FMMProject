
struct FMM{N,K<:AbstractKernel{N}}
    tree::TreeNode{N}
    levels::Vector{Vector{TreeNode{N}}}
    points::Vector{Point{N}}
    result::Vector{ComplexF64}
    interactionrank::Integer    # P
    maxlevel::Integer           # L
    maxpointsperleaf::Integer   # M
end

const FMMLaplace{N} = FMM{N,<:Laplace{N}} where {N}
const FMMLaplace2D = FMM{2,<:Laplace2D}
const FMMLaplace3D = FMM{3,<:Laplace3D}

kernel(::FMM{N,K}) where {N,K} = K
levels(f::FMM) = f.levels
level(f::FMM,i::Int64) = levels(f)[i]
nlevels(f::FMM) = length(levels(f))
maxlevel(f::FMM) = f.maxlevel
points(f::FMM) = f.points 
maxpointsperleaf(f::FMM) = f.maxpointsperleaf
interaction_rank(f::FMM) = f.interaction_rank
eachlevel(f::FMM) = (l,level(f,l) for l in 1:nlevels(f))
leaves(f::FMM) = level(f,maxlevel(f))

is_new_node_a_leaf(fmm::FMM, npoints) = maxlevel(fmm)==nlevels(fmm)+1  # for the moment
