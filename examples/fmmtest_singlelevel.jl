using Random
using LinearAlgebra
using FMMProject
import FMMProject: Laplace2D,convert_to_singlelevel!,singlelevel_forwardmap!

##
Random.seed!(1)
npoints = 1000
points = rand(Point2D,npoints)
qcharges = ones(npoints)
qcharges /= norm(qcharges,1)  # normalize charges

ntest = 256
idxs  = rand(1:npoints,ntest)

## Forward map
K = Laplace2D
pot = [sum(K(points[i],y)*qy for (y,qy) in zip(points,qcharges)) for i in idxs]

## fmm
L = 4 # max level   Lapprox = log(npoints/10^2)/log(4)+1
P = 5 # interaction rank  max=34
fmm = FMMLaplace2D(points,P,L);
convert_to_singlelevel!(fmm);

##
pot_approx = singlelevel_forwardmap!(fmm,qcharges)

pot_test = pot_approx[idxs]

err = norm(real.(pot_test)-pot,Inf)/norm(pot,Inf)
@info "" err  
