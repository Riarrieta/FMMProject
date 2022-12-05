
function Base.:*(fmm::FMM,qvec::AbstractVector)
    @assert length(qvec) == npoints(fmm) == length(result(fmm))
    upward_pass!(fmm,qvec)
    downward_pass!(fmm)


    return result(fmm)
end

function compute_outgoing_expansion!(τ::TreeNode,qbuffer,qvec)
    fill!(τ.qhat,0)
    Nτ = npoints(τ)
    iszero(Nτ) && return
    qτ = @view qbuffer[1:Nτ]   # buffer
    qτ .= qvec[points_indices(τ)]  # copy charges in buffer
    if isleaf(τ)
        mul!(τ.qhat,τ.Tofs,qτ)  #  τ.qhat = τ.Tofs*qτ
    else
        for (child,Tofo) in zip(children(τ),τ.Tofo)
            qhat_child = child.qhat
            mul!(τ.qhat,Tofo,qhat_child,true,true)  # τ.qhat = Tofo*qhat_child + τ.qhat
        end
    end
end
function upward_pass!(fmm::FMM,qvec::AbstractVector)
    qbuffer = # use 'result' as buffer for charges
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

function compute_potential!(fmm::FMM,leaf::TreeNode,qvec)
    # Fix: function barriers
    leaf_result = leaf.Ttfi*leaf.vhat   # FIX: do not allocate
    K = kernel(fmm)
    xpoints = points(leaf)
    xindices = points_indices(leaf)
    r = result(fmm)
    # compute interactions 
    for (x,ix) in zip(xpoints,xindices)
        # target-from-incoming
        r[ix] += leaf_result[ix]  
        # self-interactions
        for (y,iy) in zip(xpoints,xindices)
            qy = qvec[iy]  # charge
            r[ix] += K(x,y)*qy
        end
        # neighbor-interactions
        for neigh in neighbor_list(leaf)
            ypoints = points(neigh)
            yindices = points_indices(neigh)
            for (y,iy) in zip(ypoints,yindices)
                qy = qvec[iy]  # charge
                r[ix] += K(x,y)*qy
            end
        end
    end
end
function compute_potential!(fmm::FMM,qvec)
    fill!(result(fmm),0)  # reset result
    for leaf in leaves(fmm)
        compute_potential!(fmm,leaf,qvec)
    end
end
