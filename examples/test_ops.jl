using Test
using Random
using LinearAlgebra
using StaticArrays
using FMMProject
import FMMProject: Laplace2D,from_point2d_to_complex,minus1exp,compute_Tifo_ops!,
                   compute_Tifo_ops!,compute_Tofo_ops!,compute_Tifi_ops!,
                   compute_Tofs_ops!,compute_Ttfi_ops!

##

ζ = 6    # distance between source/target boxes
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
P = 5

## Tofs op, working!
Nσ = N
Tofs = zeros(ComplexF64,P,Nσ)
compute_Tofs_ops!(Tofs,from_point2d_to_complex(ypoints),cσ,P)
qhat = Tofs*qcharges

pot_approx1 = [log(x-cσ)*qhat[1] + sum(1/(x-cσ)^(p-1)*qhat[p] for p in 2:P) for x in from_point2d_to_complex(xpoints)]
err1_real = norm(real.(pot_approx1)-pot,Inf)/norm(pot,Inf)
@info "Tofs" err1_real

## Tifo op and Ttfi op, working!
cτ = cx  # center of expansion

Tifo = zeros(ComplexF64,P,P)
compute_Tifo_ops!(Tifo,cσ,cτ,P)
vhat = Tifo*qhat

xlist = xpoints |> from_point2d_to_complex
Nτ = M
Ttfi = zeros(ComplexF64,Nτ,P)
compute_Ttfi_ops!(Ttfi,xlist,cτ,P)

pot_approx2 = Ttfi*vhat
err2_real = norm(real.(pot_approx2)-pot,Inf)/norm(pot,Inf)
@info "Tifo + Ttfi" err2_real

## Tofo op, new center, working!
new_cσ = cσ + (0.1 - 0.1im)

new_Tofs = zeros(ComplexF64,P,Nσ)
compute_Tofs_ops!(new_Tofs,from_point2d_to_complex(ypoints),new_cσ,P)
qhat_new_cσ = new_Tofs*qcharges

pot_approx3 = [log(x-new_cσ)*qhat_new_cσ[1] + sum(1/(x-new_cσ)^(p-1)*qhat_new_cσ[p] for p in 2:P) for x in from_point2d_to_complex(xpoints)]
err3_real = norm(real.(pot_approx3)-pot,Inf)/norm(pot,Inf)
#@info "" err3_real

Tofo = LowerTriangular(zeros(ComplexF64,P,P))
compute_Tofo_ops!(Tofo,cσ,new_cσ,P)
approx_qhat_new_cσ = Tofo*qhat

err4 = norm(approx_qhat_new_cσ-qhat_new_cσ,Inf)/norm(qhat_new_cσ,Inf)
@info "Tofo" err4   # must be exact

## Tifi op, new center
new_cτ = cx + (0.1-0.1im) # center of expansion

new_Tifo = zeros(ComplexF64,P,P)
compute_Tifo_ops!(new_Tifo,cσ,new_cτ,P)
new_vhat = new_Tifo*qhat  # with respect to new_cτ

new_Ttfi = zeros(ComplexF64,Nτ,P)
compute_Ttfi_ops!(new_Ttfi,xlist,new_cτ,P)   # with respect to new_cτ

pot_approx5 = new_Ttfi*new_vhat
err5_real = norm(real.(pot_approx5)-pot,Inf)/norm(pot,Inf)
#@info "" err5_real

Tifi = UpperTriangular(zeros(ComplexF64,P,P))
compute_Tifi_ops!(Tifi,cτ,new_cτ,P)
approx_vhat = Tifi*new_vhat;   # with respect to cτ

err6 = norm(Ttfi*approx_vhat-new_Ttfi*new_vhat,Inf)/norm(new_Ttfi*new_vhat,Inf)
@info "Tifi" err6  # must be exact


