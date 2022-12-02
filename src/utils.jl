
abstract type AbstractKernel end
abstract type Laplace <: AbstractKernel end

const Point{N} = SVector{N,Float64}
const Point2D  = Point{2}
const Point3D  = Point{3}

from_point2d_to_complex(l::Vector{Point2D}) = reinterpret(ComplexF64,l)

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


