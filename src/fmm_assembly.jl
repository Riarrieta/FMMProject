
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


###### Laplace 2D ################
function FMMLaplace(points::Vector{Point2D},P,M)
    K = Laplace2D   # kernel
    tree = assemble_root_tree_node(points,K,P,M)
    levels = [[tree]]
    result = zeros(ComplexF64,length(points))
    fmm = FMM{2,P,K}(tree,levels,points,result,M)
    return fmm
end