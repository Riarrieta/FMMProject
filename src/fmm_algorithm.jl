
function Base.:*(fmm::FMM,qvec::AbstractVector)
    @assert length(qvec) == npoints(fmm) == length(result(fmm))
    qbuffer = deepcopy(result(fmm))
    upward_pass!(fmm,qvec,qbuffer)
    downward_pass!(fmm)
    compute_potential!(fmm,qvec,qbuffer)
    return result(fmm)
end

function compute_outgoing_expansion!(τ::TreeNode,qbuffer,qvec)
    if isleaf(τ)
        Nτ = npoints(τ)
        iszero(Nτ) && return
        qτ = @view qbuffer[1:Nτ]   # buffer
        qτ .= qvec[points_indices(τ)]  # copy charges in buffer
        mul!(τ.qhat,τ.Tofs,qτ)  #  τ.qhat = τ.Tofs*qτ
    else
        fill!(τ.qhat,0)
        for (child,Tofo) in zip(children(τ),τ.Tofo)
            qhat_child = child.qhat
            mul!(τ.qhat,Tofo,qhat_child,true,true)  # τ.qhat = Tofo*qhat_child + τ.qhat
        end
    end
end
function upward_pass!(fmm::FMM,qvec::AbstractVector,qbuffer)
    for (_,level) in Iterators.reverse(eachlevel(fmm))
        for tree in level
            compute_outgoing_expansion!(tree,qbuffer,qvec)
        end
    end
end

function compute_incoming_expansion!(τ::TreeNode)
    σ = parent(τ)
    mul!(τ.vhat,τ.Tifi,σ.vhat)  # τ.vhat = τ.Tifi*σ.vhat
    for (i,Tifo) in zip(interaction_list(τ),τ.Tifo)
        mul!(τ.vhat,Tifo,i.qhat,true,true)    # τ.vhat = Tifo*i.qhat + τ.vhat
    end
end
function downward_pass!(fmm::FMM)
    # set incoming expansions to zero in level 2
    for tree in level(fmm,2)
        fill!(tree.vhat,0)
    end
    for (l,level) in eachlevel(fmm)
        l ≤ 2 && continue  # only level 3 onwards
        for tree in level
            compute_incoming_expansion!(tree)
        end
    end
end

function compute_potential!(fmm::FMM,leaf::TreeNode,qvec,qbuffer)
    # Fix: function barriers
    Nτ = npoints(leaf)
    iszero(Nτ) && return
    leaf_result = @view qbuffer[1:Nτ]   # buffer
    mul!(leaf_result,leaf.Ttfi,leaf.vhat)  # leaf_result = leaf.Ttfi*leaf.vhat
    K = kernel(fmm)
    r = result(fmm)
    # compute interactions 
    for (ilocal_x,iglobal_x,x) in eachpoint(leaf)
        # target-from-incoming
        r[iglobal_x] += leaf_result[ilocal_x]  
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
function compute_potential!(fmm::FMM,qvec,qbuffer)
    fill!(result(fmm),0)  # reset result
    for leaf in leaves(fmm)
        compute_potential!(fmm,leaf,qvec,qbuffer)
    end
end
