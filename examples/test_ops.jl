using Test
using Random
using LinearAlgebra
using StaticArrays
using FMMProject
import FMMProject: Laplace2D,from_point2d_to_complex,minus1exp,compute_Tifo_ops!,
                   compute_Tifo_ops!,compute_Tofo_ops!,compute_Tifi_ops!,
                   compute_Tofs_ops!,compute_Ttfi_ops!

##

ζ = 2    # distance between source/target boxes
N = M = 1000   # number of sources/targets
qcharges = ones(N)    # charges at sources (unit charges)
xpoints =  rand(Point2D,M) # target points in [0,1]²
ypoints =  [rand(Point2D)+SVector(ζ,0) for _ in 1:N]   # source points in [ζ,ζ+1]×[0,1]
K = Laplace2D

ϕmat = [K(x,y) for x in xpoints, y in ypoints]
pot = ϕmat*qcharges  # exact potential

## setup
cx = Point2D(0.5,0.5) |> from_point2d_to_complex # center targets
cσ = Point2D(0.5+ζ,0.5) |> from_point2d_to_complex # center sources
P = 15

## Tofs op, working!
Nσ = N
Tofs = zeros(ComplexF64,P,Nσ)
compute_Tofs_ops!(Tofs,from_point2d_to_complex(ypoints),cσ,P)
qhat = Tofs*qcharges

pot_approx1 = [log(x-cσ)*qhat[1] + sum(1/(x-cσ)^(p-1)*qhat[p] for p in 2:P) for x in from_point2d_to_complex(xpoints)]
err1 = norm(pot_approx1-pot,Inf)/norm(pot,Inf)
err1_real = norm(real.(pot_approx1)-pot,Inf)/norm(pot,Inf)
@info "" err1 err1_real

## Tifo op and Ttfi op, working!
cτ = cx #+ (0.1-0.0im) # center of expansion

Tifo = zeros(ComplexF64,P,P)
compute_Tifo_ops!(Tifo,cσ,cτ,P)
vhat = Tifo*qhat

xlist = xpoints |> from_point2d_to_complex
Nτ = M
Ttfi = zeros(ComplexF64,Nτ,P)
compute_Ttfi_ops!(Ttfi,xlist,cτ,P)

pot_approx2 = Ttfi*vhat
err2 = norm(pot_approx2-pot,Inf)/norm(pot,Inf)
err2_real = norm(real.(pot_approx2)-pot,Inf)/norm(pot,Inf)
@info "" err2 err2_real

##



