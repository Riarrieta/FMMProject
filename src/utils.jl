
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



