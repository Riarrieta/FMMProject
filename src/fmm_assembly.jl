
function initialize_tree_structures(fmm::FMM, npoints)
    K = kernel(fmm)
    P = interaction_rank(fmm)
    isleaf = is_new_node_a_leaf(fmm, npoints)
    return initialize_tree_structures(K, npoints, P; isleaf)
end

function assemble_root_tree_node(points,K,P,M)
    npoints = length(points)
    parent = nothing
    box = compute_bounding_box(points)
    children,nlist,ilist,qhat,vhat,
             Tofo,Tifo,Tifi,Tofs,Ttfi = initialize_tree_structures(K,npoints,P;isleaf=false)
    points_indices = collect(1:npoints)
    copy_points = deepcopy(points)  # make copy of points
    tree = TreeNode(parent,children,nlist,ilist,box,qhat,vhat,
                    Tofo,Tifo,Tifi,Tofs,Ttfi,copy_points,points_indices)
    return tree                    
end

function assemble_child_tree_node(fmm::FMM{N},
                                  parent::TreeNode{N},
                                  points::Vector{Point{N}},
                                  points_indices::Vector{Int64},
                                  box::Box{N}) where N
    npoints = length(points)
    children,nlist,ilist,qhat,vhat,Tofo,Tifo,Tifi,Tofs,Ttfi = initialize_tree_structures(fmm, npoints)
    tree = TreeNode(parent,children,nlist,ilist,box,qhat,vhat,
                    Tofo,Tifo,Tifi,Tofs,Ttfi,points,points_indices)
    return tree                    
end

# for the moment, tree nodes split until max level is reached
shouldsplit(fmm::FMM{N},tree::TreeNode{N}) where N = true 

function split_tree!(fmm::FMM{N},tree::TreeNode{N}) where N
    shouldsplit(fmm,tree) || return
    xpoints = points(tree)
    xindices = points_indices(tree)
    nodebox = box(tree)
    children_list = children(tree)
    children_points,children_indices = split_points(xpoints,xindices,nodebox)
    children_boxes = split_box(nodebox)
    # empty unnecessary data in current tree
    empty!(xpoints)
    empty!(xindices)
    # assemble children
    for (cpoints,cindices,cbox) in zip(children_points,children_indices,children_boxes)
        child = assemble_child_tree_node(fmm,tree,cpoints,cindices,cbox)
        push!(children_list,child)
    end
end

function assemble_interaction_and_neigh_lists!(new_level::Vector{TreeNode{N}},
                                              is_last_level::Bool) where N
    nnodes = length(new_level)
    for i in 1:nnodes
        node_i = new_level[i]
        for j in i+1:nnodes
            node_j = new_level[j]
            if are_a_and_b_in_interaction_list(node_i,node_j)
                push!(interaction_list(node_i),node_j)
                push!(interaction_list(node_j),node_i)
            end
            if is_last_level && !is_a_well_separated_from_b(node_i,node_j)
                push!(neighbor_list(node_i),node_j)
                push!(neighbor_list(node_j),node_i)
            end
        end
    end
end

function assemble_fmm!(fmm::FMM{N}) where N
    for l in 1:maxlevel(fmm)-1
        level_trees = level(fmm,l)
        # add new level
        new_level = TreeNode{N}[]
        push!(levels(fmm),new_level)
        # split trees and add children to new level
        for tree in level_trees
            split_tree!(fmm,tree)
            children_list = children(tree)
            append!(new_level,children_list)
        end
        is_last_level = l==maxlevel(fmm)-1
        assemble_interaction_and_neigh_lists!(new_level,is_last_level)
    end
    @assert maxlevel(fmm) == nlevels(fmm)
end

function FMM(K,points::Vector{Point{N}},P,L,M=100) where N
    @assert P > 0
    @assert L > 1
    @assert M > 0
    tree = assemble_root_tree_node(points,K,P,M)
    levels = [[tree]]
    result = zeros(ComplexF64,length(points))
    fmm = FMM{2,K}(tree,levels,points,result,P,L,M)
    assemble_fmm!(fmm)
    return fmm
end

###### Laplace 2D ################
FMMLaplace2D(points::Vector{Point2D},P,L,M=100) = FMM(Laplace2D,points,P,L,M)
