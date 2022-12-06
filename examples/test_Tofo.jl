using LinearAlgebra
using Random
using FMMProject
import FMMProject: compute_Tofo_ops!

N = 100
P = 3
cσ = rand(ComplexF64)
cτ = rand(ComplexF64)
ypoints = rand(ComplexF64,N)
qcharges = randn(N)

##
qhatσ = [sum((-1/p)*(ypoints[j]-cσ)^p*qcharges[j] for j in 1:N) for p in 0:P-1] #
qhatσ[1] = sum(qcharges)

qhatτ = [sum((-1/p)*(ypoints[j]-cτ)^p*qcharges[j] for j in 1:N) for p in 0:P-1]
qhatτ[1] = sum(qcharges)

##
Tofo_manual = [s/p*binomial(p,s)*(cσ-cτ)^(p-s) for p in 0:P-1, s in 0:P-1]  # s/p*
Tofo_manual[:,1] = [(-1/p)*(cσ-cτ)^p for p in 0:P-1]
Tofo_manual[1,:] = [i==1 for i in 1:P]

#compute_Tofo_ops!(Tofo,cσ,cτ,P)

@info "correct" norm(qhatτ-Tofo_manual*qhatσ,Inf)
@info "incorrect" norm(qhatσ-Tofo_manual*qhatτ,Inf)

