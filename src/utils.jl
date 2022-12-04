
abstract type AbstractKernel{N} end
abstract type Laplace{N} <: AbstractKernel{N} end

abstract type Laplace2D <: Laplace{2} end
(::Type{Laplace2D})(x,y) = log(norm(x-y))

abstract type Laplace3D <: Laplace{3} end
(::Type{Laplace3D})(x,y) = 1/(4π*norm(x-y))

dimension(::Type{AbstractKernel{N}}) where N = N

const Point{N} = SVector{N,Float64}
const Point2D  = Point{2}
const Point3D  = Point{3}

from_point2d_to_complex(p::Point2D) = ComplexF64(p...)
from_point2d_to_complex(l::Vector{Point2D}) = reinterpret(ComplexF64,l)

minus1exp(p) = (-1)^isodd(p)

# Box of dimension N
struct Box{N}
    lc::Point{N}  # lower corner 
    uc::Point{N}  # upper corner 
    ce::Point{N}  # center
    a::Float64    # half side length
end
function Box(ce::Point{N},a::Float64) where N
    @assert a ≥ 0
    avec = ntuple(i -> a, Val(N))
    uc = ce .+ avec
    lc = ce .- avec
    return Box(lc,uc,ce,a)
end

corners(b::Box) = b.lc,b.uc
center(b::Box) = b.ce
halfside(b::Box) = b.a
dist(b1::Box{N},b2::Box{N}) where N = norm(center(b1)-center(b2),Inf)

isinbox(p::Point{N},b::Box{N}) where N = all(b.lc .≤ p .≤ b.uc)

function compute_bounding_box(points::Vector{Point{N}}) where N
    lc = Point{N}(ntuple(i -> Inf,N))
    uc = Point{N}(ntuple(i -> -Inf,N))
    for p in points
        lc = min.(lc,p)
        uc = max.(uc,p)
    end
    ce = (lc + uc) / 2   # center
    a = maximum(uc - lc) / 2 # half side
    a = a * 1.1     # increase half side a bit
    return Box(ce,a)
end

# Split a vector of points in 2^N vectors of points
# with respect to a central point 'ce'
function split_points(points::Vector{Point{N}},
                      point_indices::Vector{Int64},
                      ce::Point{N}) where N
    vecpoints = [Point{N}[] for _ in 1:2^N]
    matpoints = reshape(vecpoints,ntuple(_->2,N))
    vecindices = [Int64[] for _ in 1:2^N]
    matindices = reshape(vecindices,ntuple(_->2,N))
    for (global_index,p) in zip(point_indices,points)
        index = ntuple(N) do i
            1 + (p[i] > ce[i])
        end
        push!(matpoints[index...],p)
        push!(matindices[index...],global_index)
    end
    return vecpoints,vecindices
end
split_points(points::Vector{Point{N}},
             point_indices::Vector{Int64},
             b::Box{N}) where N = split_points(points,point_indices,center(b))

# Split a Box in 2^N equal Box'es
# should be consistent with the function 'split_points'
function split_box(b::Box{N}) where N
    ce = center(b)
    a = halfside(b)
    children_a  = a/2
    children_ce = [ce for _ in 1:2^N]
    mat_children_ce = reshape(children_ce,ntuple(_->2,N))
    for index in CartesianIndices(mat_children_ce)
        offset = ntuple(N) do i
            children_a*(-1)^isodd(index[i])
        end
        mat_children_ce[index] += Point{N}(offset)
    end
    children_boxes = [Box(c,children_a) for c in children_ce]
    return children_boxes
end


