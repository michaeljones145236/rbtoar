# Restarted Block TOAR
A pure Julia implementation of a restarted block version of the TOAR algorithm for the solution of large sparse quadratic eigenvalue problems. Currently a WIP, some features are missing and it is largely untested.

## Usage
```julia
    quadEigRBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix, req::Int=100, tol::Float64=1e-12, kℓₘₐₓ::Int, ℓ::Int; step::Int=10, σ::Union{Float64,ComplexF64}=0.0+0.0im, inv::Bool=true, keep::Function=every, dtol::Float64=1e-10, rrv::Bool=false, arpack::Bool=true, flvd::Bool=true, verb::Int=0, check_singular::Bool=false)
```
The function builds a second-order Krylov subspace, restarting when necessary, until the desired number of eigenpairs can be recovered with an acceptable backward error.

### Arguments
-`M`, `D` and `K` are the mass, damping and stiffness matrices respectively from the QEP `(λ²M+λD+K)x=0`. These should typically be given in `SparseMatrixCSC` format, although any `AbstractMatrix` type is allowed.
-`req` is the number of eigenpairs required. This should not be set larger than half of `kℓₘₐₓ` or else the algorithm may stagnate, that is, it may get stuck in a loop of expanding the subspace and restarting without ever reaching `req` eigenpairs. If the domain of interest is particularly awkward, it may be beneficial to set `req` even lower in relation to `kℓₘₐₓ`.
-`tol` is the tolerance for the backward error residual in computed eigenpairs. It is not reccommended to set this smaller than $10^{-13}$. Backward error residuals are computed as 
