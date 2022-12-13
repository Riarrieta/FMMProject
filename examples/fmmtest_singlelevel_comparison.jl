using Random
using LinearAlgebra
using FMMProject
import FMMProject: Laplace2D,
                   convert_to_singlelevel!,
                   singlelevel_forwardmap!,
                   singlelevel_forwardmap_without_incoming!

##
Random.seed!(1)
npoints = 1000
points = rand(Point2D,npoints)
qcharges = randn(npoints)#ones(npoints)
qcharges /= norm(qcharges,1)  # normalize charges

ntest = 256
idxs  = rand(1:npoints,ntest)

## Forward map
K = Laplace2D
pot = [sum(K(points[i],y)*qy for (y,qy) in zip(points,qcharges)) for i in idxs]

## fmm
L = 6 # max level   Lapprox = log(npoints/10^2)/log(4)+1
P = 5 # interaction rank  max=34
fmm = FMMLaplace2D(points,P,L);
fmm_sl = FMMLaplace2D(points,P,L);
convert_to_singlelevel!(fmm_sl);

## Multi level
pot_approx = fmm*qcharges

pot_test = pot_approx[idxs]

err = norm(real.(pot_test)-pot,Inf)/norm(pot,Inf)
@info "" err  

## Single level
pot_approx = singlelevel_forwardmap!(fmm_sl,qcharges)

pot_test = pot_approx[idxs]

err2 = norm(real.(pot_test)-pot,Inf)/norm(pot,Inf)
@info "" err2

##
vl = []
tl = []
p1 = []
p2 = []
@assert length(FMMProject.leaves(fmm)) == length(FMMProject.leaves(fmm_sl))
for (l1,l2) in zip(FMMProject.leaves(fmm),FMMProject.leaves(fmm_sl))
    @assert l1.qhat == l2.qhat
    @assert l1.Tifi == l2.Tifi
    @assert l1.Ttfi == l2.Ttfi
    push!(vl,norm(l1.vhat-l2.vhat,Inf))
    push!(tl,norm(real.(l1.Ttfi*l1.vhat-l2.Ttfi*l2.vhat),Inf))
    append!(p1,real.(l1.Ttfi*l1.vhat))
    append!(p2,real.(l2.Ttfi*l2.vhat))
end

## try leaf
leave_index = 1
@assert FMMProject.leaves(fmm)[leave_index].box.ce == FMMProject.leaves(fmm_sl)[leave_index].box.ce

leaf = FMMProject.leaves(fmm)[leave_index]
parent = leaf.parent
farfield = setdiff(FMMProject.leaves(fmm),(leaf,leaf.ilist...,leaf.nlist...));
farparent = (f.parent for f in farfield) |> Set |> collect;
@assert length(farfield) == length(FMMProject.leaves(fmm))-1-length(leaf.ilist)-length(leaf.nlist)

## outgoing-from-outgoing
for farparent in farparent
    @assert length(farparent.children) == length(farparent.Tofo)
    farparent.qhat .= sum(Tofo*c.qhat for (Tofo,c) in zip(farparent.Tofo,farparent.children))
end

## incoming-from-outgoing far field  # works!
parent_vhat = sum(Tifo*c.qhat for (Tifo,c) in zip(parent.Tifo,parent.ilist))
vhat = leaf.Tifi*parent_vhat  # vhat = zeros(ComplexF64,P)

## incoming-from-outgoing interaction list   # works!
for (i,Tifo) in zip(leaf.ilist,leaf.Tifo)
    mul!(vhat,Tifo,i.qhat,true,true)    # vhat = Tifo*i.qhat + vhat
end

## targets-from-incoming (only far field contribution)
pot1 = leaf.Ttfi*leaf.vhat # multi level, farfield+interaction
pot2 = leaf.Ttfi*vhat   # try
pot3 = leaf.Ttfi*FMMProject.leaves(fmm_sl)[leave_index].vhat  # single level, farfield+interaction

@info "multi vs single" norm(real.(pot1-pot3),Inf)
@info "try vs single" norm(real.(pot2-pot3),Inf)
@info "try vs multi" norm(real.(pot2-pot1),Inf)

## compare single leaf


