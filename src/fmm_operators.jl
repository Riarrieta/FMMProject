
###### Laplace 2D ################
function initialize_tree_structures(::Type{Laplace2D}, npoints, P; isleaf)
    children = TreeNode2D[]   # list of children
    nlist = TreeNode2D[]      # neighbor list
    ilist = TreeNode2D[]      # interaction list
    qhat = zeros(ComplexF64,P)        # outcoming expansion
    vhat = zeros(ComplexF64,P)         # incoming expansion
    # Translation operators
    Tofo = LowerTriangular{ComplexF64,Matrix{ComplexF64}}[]  # outgoing-from-outgoing
    Tifo = Matrix{ComplexF64}[]  # incoming-from-outgoing
    Tifi = UpperTriangular(zeros(ComplexF64,P,P))  # incoming-from-incoming
    if isleaf
        Tofs = zeros(ComplexF64,P,npoints)
        Ttfi = zeros(ComplexF64,npoints,P)
    else
        Tofs = zeros(ComplexF64,0,0)
        Ttfi = zeros(ComplexF64,0,0)
    end
    return children,nlist,ilist,qhat,vhat,Tofo,Tifo,Tifi,Tofs,Ttfi
end

function compute_Tofo_ops!(fmm::FMMLaplace2D,t::TreeNode{2})
    P = interaction_rank(fmm)
    childlist = children(t)
    empty!(t.Tofo)
    append!(t.Tofo, LowerTriangular(zeros(ComplexF64,P,P)) for _ in childlist)
    cτ = center(t) |> from_point2d_to_complex
    for (child,mat) in zip(childlist,t.Tofo)
        cσ = center(child) |> from_point2d_to_complex
        d = cσ - cτ
        for j in 1:P
            s = j-1
            for i in j:P
                r = i-1
                mat[i,j] = binomial(r,s)*d^(r-s)
            end
        end
    end
end

function compute_Tifo_ops!(Tifo::Matrix{ComplexF64},
                           cσ::ComplexF64,
                           cτ::ComplexF64,
                           P::Int64)
    d = cσ - cτ
    for j in 1:P
        p = j-1
        for i in 1:P
            r = i-1
            if p==0 && r==0
                Tifo[i,j] = log(-d)
            elseif r==0
                Tifo[i,j] = minus1exp(p)/d^p
            elseif p==0
                Tifo[i,j] = -1/(r*d^r)
            else
                Tifo[i,j] = minus1exp(p)*binomial(r+p-1,p-1)/d^(r+p)
            end
        end
    end
end
function compute_Tifo_ops!(fmm::FMMLaplace2D,τ::TreeNode{2})
    P = interaction_rank(fmm)
    ilist = interaction_list(τ)
    cτ = center(τ) |> from_point2d_to_complex
    empty!(τ.Tifo)
    append!(τ.Tifo, zeros(ComplexF64,P,P) for _ in ilist)
    for (σ,mat) in zip(ilist,τ.Tifo)
        cσ = center(σ) |> from_point2d_to_complex
        compute_Tifo_ops!(mat,cσ,cτ,P)
    end
end

function compute_Tifi_ops!(Tifi::UpperTriangular,
                           cσ::ComplexF64,
                           cτ::ComplexF64,
                           P::Int64)
    d = cσ - cτ
    for j in 1:P
        p = j-1
        for i in 1:j
            r = i-1
            Tifi[i,j] = binomial(p,r)*d^(p-r)
        end
    end
end
function compute_Tifi_ops!(fmm::FMMLaplace2D,σ::TreeNode{2})
    isroot(σ) && return
    P = interaction_rank(fmm)
    τ = parent(σ)
    cσ = center(σ) |> from_point2d_to_complex
    cτ = center(τ) |> from_point2d_to_complex
    compute_Tifi_ops!(σ.Tifi,cσ,cτ,P)
end

function compute_Tofs_ops!(Tofs::Matrix{ComplexF64},
                           xlist::AbstractVector{ComplexF64},
                           cσ::ComplexF64,
                           P::Int64)
    Nσ = length(xlist)
    for j in 1:Nσ
        xj = xlist[j]
        d = xj-cσ
        Tofs[1,j] = one(ComplexF64)
        for i in 2:P
            p = i-1
            Tofs[i,j] = -1/p*d^p
        end
    end
end
function compute_Tofs_ops!(fmm::FMMLaplace2D,σ::TreeNode{2})
    P = interaction_rank(fmm)
    cσ = center(σ) |> from_point2d_to_complex
    xlist = points(σ) |> from_point2d_to_complex
    compute_Tofs_ops!(σ.Tofs,xlist,cσ,P)
end

function compute_Ttfi_ops!(Ttfi::Matrix{ComplexF64},
                           xlist::AbstractVector{ComplexF64},
                           cτ::ComplexF64,
                           P::Int64)
    Nτ = length(xlist)
    for j in 1:P
        p = j-1
        for i in 1:Nτ
            xi = xlist[i] 
            d = xi-cτ
            Ttfi[i,j] = d^p
        end
    end
end
function compute_Ttfi_ops!(fmm::FMMLaplace2D,τ::TreeNode{2})
    P = interaction_rank(fmm)
    cτ = center(τ) |> from_point2d_to_complex
    xlist = points(τ) |> from_point2d_to_complex
    compute_Ttfi_ops!(τ.Ttfi,xlist,cτ,P)
end
