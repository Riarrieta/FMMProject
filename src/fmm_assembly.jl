
function assemble_root_tree_node(points,K,P,M)
    npoints = length(points)
    parent = Nothing
    box = compute_bounding_box(points)
    children,ilist,qhat,vhat,Tofo,Tifo,Tifi,Tofs,Ttfi = initialize_tree_structures(K,npoints,P,M)
    points_indices = collect(1:npoints)
    tree = TreeNode(parent,children,ilist,box,qhat,vhat,
                    Tofo,Tifo,Tifi,Tofs,Ttfi,deepcopy(points),points_indices)
    return tree                    
end

shouldsplit(fmm::FMM{N},tree::TreeNode{N}) where N = npoints(tree)>maxpointsperleaf(fmm)

function split_tree(fmm::FMM{N},tree::TreeNode{N},level=0) where N
    shouldsplit(fmm,tree) || return
    xpoints = points(tree)
    xindices = points_indices(tree)
    split_points(points::Vector{Point{N}},ce::Point{N})
    
end

function FMM(K,points::Vector{Point{N}},P,M) where N
    tree = assemble_root_tree_node(points,K,P,M)
    levels = [[tree]]
    result = zeros(ComplexF64,length(points))
    fmm = FMM{2,P,K}(tree,levels,points,result,M)
    return fmm
end

###### Laplace 2D ################
FMMLaplace(points::Vector{Point2D},P,M) = FMM(Laplace2D,points,P,M)
