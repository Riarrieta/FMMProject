function convert_to_singlelevel!(fmm::FMMLaplace2D)
    for leaf in leaves(fmm)
        # use interaction list as far-field list
        empty!(interaction_list(leaf))
        # we remove the neighbor list and the leaf itself
        append!(interaction_list(leaf),setdiff(leaves(fmm),(neighbor_list(leaf)...,leaf)))
        @assert length(leaves(fmm))-9 ≤ length(interaction_list(leaf)) ≤ length(leaves(fmm))-4
    end
end

function singlelevel_forwardmap!(fmm::FMM,qvec::AbstractVector)
    @assert length(qvec) == npoints(fmm) == length(result(fmm))
    # construct outgoing expansions for leaves
    qbuffer = result(fmm)  # use 'result' as buffer for charges
    for leaf in leaves(fmm)
        compute_outgoing_expansion!(leaf,qbuffer,qvec)
    end
    # contruct incoming expansion
    for leaf in leaves(fmm)
        fill!(leaf.vhat,0)
        for (i,Tifo) in zip(interaction_list(leaf),leaf.Tifo)
            mul!(leaf.vhat,Tifo,i.qhat,true,true)    # τ.vhat = Tifo*i.qhat + τ.vhat
        end
    end
    # contruct potentials
    compute_potential!(fmm,qvec)
    return result(fmm)
end 