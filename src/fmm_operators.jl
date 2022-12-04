
###### Laplace 2D ################
function initialize_tree_structures(::Type{Laplace2D}, npoints, P, M)
    children = TreeNode2D[]   # list of children
    ilist = TreeNode2D[]      # interaction list
    qhat = zeros(ComplexF64,P)        # outcoming expansion
    vhat = zeros(ComplexF64,P)         # incoming expansion
    # Translation operators
    Tofo = LowerTriangular{ComplexF64,Matrix{ComplexF64}}[]  # outgoing-from-outgoing
    Tifo = Matrix{ComplexF64}[]  # incoming-from-outgoing
    Tifi = UpperTriangular(zeros(ComplexF64,P,P))  # incoming-from-incoming
    if npoints ≤ M   # leaf
        Tofs = zeros(ComplexF64,P,npoints)
        Ttfi = zeros(ComplexF64,npoints,P)
    else  # no leaf
        Tofs = zeros(ComplexF64,0,0)
        Ttfi = zeros(ComplexF64,0,0)
    end
    return children,ilist,qhat,vhat,Tofo,Tifo,Tifi,Tofs,Ttfi
end

function compute_Tofo_ops(::FMMLaplace2D,t::TreeNode{2})
    isleaf(t) && return
    P = interaction_rank(t)
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

function compute_Tifo_ops(::FMMLaplace2D,t::TreeNode{2})
    P = interaction_rank(t)
    ilist = interaction_list(τ)
    cτ = center(τ) |> from_point2d_to_complex
    empty!(τ.Tifo)
    append!(τ.Tifo, zeros(ComplexF64,P,P) for _ in ilist)
    for (σ,mat) in zip(ilist,τ.Tifo)
        cσ = center(σ) |> from_point2d_to_complex
        d = cσ - cτ
        for j in 1:P
            p = j-1
            for i in 1:P
                r = i-1
                if p==0 && r==0
                    mat[i,j] = log(-d)
                elseif r==0
                    mat[i,j] = minus1exp(p)/d^p
                elseif p==0
                    mat[i,j] = -1/(r*d^r)
                else
                    mat[i,j] = minus1exp(p)*binomial(r+p-1,p-1)/d^(r+p)
                end
            end
        end
    end
end

function compute_Tifi_ops(::FMMLaplace2D,t::TreeNode{2})
    isroot(σ) && return
    P = interaction_rank(t)
    τ = parent(σ)
    cσ = center(σ) |> from_point2d_to_complex
    cτ = center(τ) |> from_point2d_to_complex
    d = cσ - cτ
    mat = σ.Tifi
    for j in 1:P
        p = j-1
        for i in 1:j
            r = i-1
            mat[i,j] = binomial(p,r)*d^(p-r)
        end
    end
end

function compute_Tofs_ops(::FMMLaplace2D,t::TreeNode{2})
    P = interaction_rank(t)
    cσ = center(σ) |> from_point2d_to_complex
    xlist = points(σ)
    Nσ = npoints(σ)
    mat = σ.Tofs
    for j in 1:Nσ
        xj = xlist[j]
        d = xj-cσ
        mat[1,j] = one(ComplexF64)
        for i in 2:P
            p = i-1
            mat[i,j] = -1/p*d^p
        end
    end
end

function compute_Ttfi_ops(::FMMLaplace2D,t::TreeNode{2})
    P = interaction_rank(t)
    cτ = center(τ) |> from_point2d_to_complex
    xlist = points(τ) 
    Nτ = npoints(τ)
    mat = τ.Ttfi
    for j in 1:P
        p = j-1
        for i in 1:Nτ
            xi = xlist[i]
            d = xi-cτ
            mat[i,j] = d^p
        end
    end
end
