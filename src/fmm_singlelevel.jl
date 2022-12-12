function convert_to_singlelevel!(fmm::FMMLaplace2D)
    L = maxlevel(fmm)
    Nleaves = 4^(L-1)
    for leaf in leaves(fmm)
        # use interaction list as far-field list
        empty!(interaction_list(leaf))
        # we remove the neighbor list and the leaf itself
        append!(interaction_list(leaf),setdiff(leaves(fmm),(neighbor_list(leaf)...,leaf)))
        @assert length(leaves(fmm))-9 ≤ length(interaction_list(leaf)) ≤ length(leaves(fmm))-4
        @assert Nleaves == 1 + length(interaction_list(leaf)) + length(neighbor_list(leaf))
        # assemble incoming-from-outgoing operators
        compute_Tifo_ops!(fmm,leaf)
    end
end

function singlelevel_forwardmap!(fmm::FMM,qvec::AbstractVector)
    @assert length(qvec) == npoints(fmm) == length(result(fmm))
    # construct outgoing expansions for leaves
    qbuffer = deepcopy(result(fmm))  # use 'result' as buffer for charges
    for leaf in leaves(fmm)
        compute_outgoing_expansion!(leaf,qbuffer,qvec)
    end
    # construct incoming expansion
    for leaf in leaves(fmm)
        fill!(leaf.vhat,0)
        @assert length(leaf.Tifo) == length(interaction_list(leaf))
        for (i,Tifo) in zip(interaction_list(leaf),leaf.Tifo)
            mul!(leaf.vhat,Tifo,i.qhat,true,true)    # τ.vhat = Tifo*i.qhat + τ.vhat
        end
    end
    # construct potentials
    compute_potential!(fmm,qvec,qbuffer)
    return result(fmm)
end 

# Without incoming versions
function compute_potential_without_incoming!(fmm::FMM,leaf::TreeNode,qvec,qbuffer)
    # Fix: function barriers
    Nτ = npoints(leaf)
    iszero(Nτ) && return
    K = kernel(fmm)
    r = result(fmm)
    P = interaction_rank(fmm)
    leaf_result_incoming = zeros(ComplexF64,Nτ)
    # compute interactions 
    for (ilocal_x,iglobal_x,x) in eachpoint(leaf)
        # far field
        leaf_result = zero(ComplexF64)
        d = zero(ComplexF64)
        @assert length(leaf.Tifo) == length(interaction_list(leaf))
        for (far,Tifo) in zip(interaction_list(leaf),leaf.Tifo)
            qhat = far.qhat
            cσ = center(far)
            d = x-cσ |> from_point2d_to_complex
            leaf_result += log(d)*qhat[1] + sum(1/(d)^(p-1)*qhat[p] for p in 2:P)
        end
        r[iglobal_x] += leaf_result 
        # self-interactions
        for (ilocal_y,iglobal_y,y) in eachpoint(leaf)
            qy = qvec[iglobal_y]  # charge
            r[iglobal_x] += K(x,y)*qy
        end
        # neighbor-interactions
        for neigh in neighbor_list(leaf)
            for (ilocal_y,iglobal_y,y) in eachpoint(neigh)
                qy = qvec[iglobal_y]  # charge
                r[iglobal_x] += K(x,y)*qy
            end
        end
    end
end
function compute_potential_without_incoming!(fmm,qvec,qbuffer)
    fill!(result(fmm),0)  # reset result
    for leaf in leaves(fmm)
        compute_potential_without_incoming!(fmm,leaf,qvec,qbuffer)
    end
end
function singlelevel_forwardmap_without_incoming!(fmm::FMM,qvec::AbstractVector)
    @assert length(qvec) == npoints(fmm) == length(result(fmm))
    # construct outgoing expansions for leaves
    qbuffer = deepcopy(result(fmm))  # use 'result' as buffer for charges
    for leaf in leaves(fmm)
        compute_outgoing_expansion!(leaf,qbuffer,qvec)
    end
    # construct potentials
    compute_potential_without_incoming!(fmm,qvec,qbuffer)
    return result(fmm)
end 
