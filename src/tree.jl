
struct TreeNode{N}
    parent::Union{Nothing,TreeNode{N}}
    children::Vector{TreeNode{N}}   # list of children
    nlist::Vector{TreeNode{N}}      # neighbor list (leaves only)
    ilist::Vector{TreeNode{N}}      # interaction list
    box::Box{N}
    qhat::Vector{ComplexF64}        # outcoming expansion
    vhat::Vector{ComplexF64}        # incoming expansion
    # Translation operators
    Tofo::Vector{LowerTriangular{ComplexF64,Matrix{ComplexF64}}}  # outgoing-from-outgoing
    Tifo::Vector{Matrix{ComplexF64}}  # incoming-from-outgoing
    Tifi::UpperTriangular{ComplexF64,Matrix{ComplexF64}}  # incoming-from-incoming
    # Leaves-only translation operators
    Tofs::Matrix{ComplexF64}
    Ttfi::Matrix{ComplexF64}
    # Leaves-only data
    points::Vector{Point{N}}           # points xᵢ in R^N
    points_indices::Vector{Int64}    # global indices of points xᵢ 
end

const TreeNode2D = TreeNode{2}
const TreeNode3D = TreeNode{3}

Base.parent(t::TreeNode) = t.parent
children(t::TreeNode) = t.children
neighbor_list(t::TreeNode) = t.nlist
interaction_list(t::TreeNode) = t.ilist
points(t::TreeNode) = t.points
points_indices(t::TreeNode) = t.points_indices
box(t::TreeNode) = t.box
npoints(t::TreeNode) = length(points(t))

# returns local_index,global_index,point
eachpoint(t::TreeNode) = ((l,points_indices(t)[l],points(t)[l]) for l in 1:npoints(t))  

corners(t::TreeNode) = corners(t.box)
center(t::TreeNode) = center(t.box)
halfside(t::TreeNode) = halfside(t.box)
dist(t1::TreeNode,t2::TreeNode) = dist(t1.box,t2.box)

isroot(t::TreeNode) = isnothing(parent(t))
isleaf(t::TreeNode) = isempty(children(t))

is_a_well_separated_from_b(a::TreeNode,b::TreeNode) = dist(a,b)≥3*halfside(b)
function are_a_and_b_in_interaction_list(a::TreeNode,b::TreeNode)
    return is_a_well_separated_from_b(a,b) &&
           !is_a_well_separated_from_b(parent(a),parent(b))
end










