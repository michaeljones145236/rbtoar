# Restarted Block TOAR
A pure Julia implementation of a restarted block version of the TOAR algorithm for the solution of large sparse quadratic eigenvalue problems. Currently a WIP, some features are missing and it is largely untested.

## Usage
```julia
    quadEigRBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix, req::Int=100, tol::Float64=1e-12, kℓₘₐₓ::Int, ℓ::Int; step::Int=10, σ::Union{Float64,ComplexF64}=0.0+0.0im, inv::Bool=true, keep::Function=every, dtol::Float64=1e-10, rrv::Bool=false, arpack::Bool=true, flvd::Bool=true, verb::Int=0, check_singular::Bool=false)
```
