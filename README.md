# Restarted Block TOAR
A pure Julia implementation of a restarted block version of the TOAR algorithm for the solution of large sparse quadratic eigenvalue problems. Currently a WIP, some features are missing and it is largely untested.

## Usage
```julia
    quadEigRBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix, req::Int=100, tol::Float64=1e-12, kℓₘₐₓ::Int, ℓ::Int; step::Int=10, σ::Union{Float64,ComplexF64}=0.0+0.0im, inv::Bool=true, keep::Function=every, dtol::Float64=1e-10, rrv::Bool=false, arpack::Bool=true, flvd::Bool=true, verb::Int=0, check_singular::Bool=false)
```
The function builds a second-order Krylov subspace, restarting when necessary, until the desired number of eigenpairs can be recovered with an acceptable backward error.

### Arguments
  - `M`, `D` and `K` are the mass, damping and stiffness matrices respectively from the QEP $(\lambda^2M+\lambda D+K)x=0$. These should typically be given in `SparseMatrixCSC` format, although any `AbstractMatrix` type is allowed. Strictly speaking, these are the only three required arguments, although typically you would want to specify at least `req` and `tol`.
  - `req` is the number of eigenpairs required. This should not be set larger than half of `kℓₘₐₓ` or else the algorithm may stagnate, that is, it may get stuck in a loop of expanding the subspace and restarting without ever reaching `req` eigenpairs. If the domain of interest is particularly awkward, it may be beneficial to set `req` even lower in relation to `kℓₘₐₓ`.
  - `tol` is the tolerance for the backward error residual in computed eigenpairs. It is not reccommended to set this smaller than $10^{-13}$. Backward error residuals of a computed eigenpair $(\tilde{\lambda},\tilde{x})$ are defined as

$$\rho=\frac{\Vert(\lambda^2M+\lambda D+K)\tilde{x}\Vert_1}{|\tilde{\lambda}|^2\Vert M\Vert_2+|\tilde{\lambda}|\Vert D\Vert_1+\Vert K\Vert_1}.$$

  - `kℓₘₐₓ` is the maximum subspace size the algorithm is permitted to build. This should be at least two times `req`, maybe even larger for a difficult domain of interest. However, if memory requirements are a problem, you may wish to lower `kℓₘₐₓ`. Setting this argument to a very high value has the effect of reducing the algorithm to non-restarted block TOAR.
  - `ℓ` is the block size/width. This should normally be set to at most 5. For larger block sizes, it might be a good idea to set `step` lower. The default value, 1, reduces the algorithm to the standard non-block restarted TOAR algorithm.
  - `step` is an optional argument specifying how many blocks should be added to the second-order Krylov subspace between each check for convergence of sufficiently many eigenpairs. Setting this too high could waste time building a larger subspace than neccessary; setting it too low could waste time solving the reduced-order QEP too many times when little progress is made enriching the subspace.
  - `σ` is the shift point in the shift-and-invert spectral transformation.

