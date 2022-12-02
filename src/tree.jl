

struct FMMStruct{N,K<:AbstractKernel}
    tree
    levels

end

struct TreeNode{N}
    parent::Union{Nothing,TreeNode{N}}
    children::Vector{TreeNode{N}}   # list of children
    ilist::Vector{TreeNode{N}}      # interaction list
    points::Vector{Point{N}}
    box::Box{N}
end

parent(t::TreeNode) = t.parent
children(t::TreeNode) = t.children
interaction_list(t::TreeNode) = t.ilist
points(t::TreeNode) = t.points
box(t::TreeNode) = t.box

corners(t::TreeNode) = corners(t.box)
center(t::TreeNode) = center(t.box)
halfside(t::TreeNode) = halfside(t.box)
dist(t1::TreeNode,t2::TreeNode) = dist(t1.box,t2.box)

isroot(t::TreeNode) = isnothing(parent(t))
isleaf(t::TreeNode) = isempty(children(t))

is_a_well_separated_from_b(a::TreeNode,b::TreeNode) = dist(a,b)≥3*halfside(b)










