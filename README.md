# Restarted Block TOAR
A pure Julia implementation of a restarted block version of the TOAR algorithm for the solution of large sparse quadratic eigenvalue problems. Currently a WIP, some features are missing and it is largely untested.

## Installation
Download `BTOAR.jl` to your Julia working directory and type
```julia-repl
julia> include("BTOAR.jl")
```
You may have to install the following dependencies:
  - `Arpack.jl`
  - `LowRankApprox.jl`
  - `QuadEig.jl`

## Usage
```julia
    quadEigRBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix, req::Int=100, tol::Float64=1e-12, kℓₘₐₓ::Int, ℓ::Int; step::Int=10, σ::Union{Float64,ComplexF64}=0.0+0.0im, inv::Bool=true, keep::Function=every, dtol::Float64=1e-10, rrv::Bool=false, arpack::Bool=true, flvd::Bool=true, verb::Int=0, check_singular::Bool=false)
```
The function builds a second-order Krylov subspace, restarting when necessary, until the desired number of eigenpairs can be recovered with an acceptable backward error.

### Arguments
  - `M`, `D` and `K` are the mass, damping and stiffness matrices respectively from the QEP $(\lambda^2M+\lambda D+K)x=0$. These should typically be given in `SparseMatrixCSC` format, although any `AbstractMatrix` type is allowed. Strictly speaking, these are the only three required arguments, although typically you would want to specify at least `req` and `tol`.
  - `req` is the number of eigenpairs required. This should not be set larger than half of `kℓₘₐₓ` or else the algorithm may stagnate, that is, it may get stuck in a loop of expanding the subspace and restarting without ever reaching `req` eigenpairs. If the domain of interest is particularly awkward, it may be beneficial to set `req` even lower in relation to `kℓₘₐₓ`.
  - `tol` is the tolerance for the backward error residual in computed eigenpairs. It is not reccommended to set this smaller than $10^{-13}$. Backward error residuals of a computed eigenpair $(\tilde{\lambda},\tilde{x})$ are defined as

$$\rho=\frac{\Vert(\tilde{\lambda}^2M+\tilde{\lambda}D+K)\tilde{x}\Vert_2}{|\tilde{\lambda}|^2\Vert M\Vert_1+|\tilde{\lambda}|\Vert D\Vert_1+\Vert K\Vert_1}.$$

  - `kℓₘₐₓ` is the maximum subspace size the algorithm is permitted to build. This should be at least two times `req`, maybe even larger for a difficult domain of interest. However, if memory requirements are a problem, you may wish to lower `kℓₘₐₓ`. Setting this argument to a very high value has the effect of reducing the algorithm to non-restarted block TOAR.
  - `ℓ` is the block size/width. This should normally be set to at most 5. For larger block sizes, it might be a good idea to set `step` lower. The default value, 1, reduces the algorithm to the standard non-block restarted TOAR algorithm.
  - `step` is an optional argument specifying how many blocks should be added to the second-order Krylov subspace between each check for convergence of sufficiently many eigenpairs. Setting this too high could waste time building a larger subspace than neccessary; setting it too low could waste time solving the reduced-order QEP too many times when little progress is made enriching the subspace.
  - `σ` is the shift point in $\mathbb{C}$ for the shift-and-invert spectral transformation. This should normally be set to a value well within the domain of interest, or convergence of wanted eigenvalues may be slow.
  - `inv` is whether to invert the QEP after shifting. This should almost always be set to `true` (the default value if unspecified) so that eigenvalues closest to `σ` converge first. If set to `false`, eigenvalues furthest from the shift point will converge first.
  - `keep` is the function that specifies the domain of interest. This function should accept a single positional argument of type `ComplexF64` and return `true` if it is in the domain of interest and false otherwise. By default, the domain of interest is taken to be the full complex plane (so `keep` always returns `true` for all inputs).
  - `dtol` is an internal numerical tolerance for breakdown/deflation detection. Normally this should not be changed, except in the case of badly behaved QEPs if you know what you're doing.
  - `rrv` is whether to use refined Ritz vectors. Refined Ritz vectors are not currently implemented and may be scrapped entirely, as it is not clear if they are ever worthwhile.
  - `arpack` is meaningless as long as refined Ritz vectors are not implemented.
  - `flvd` is whether to use Fan, Lin & Van Dooren scaling on the QEP. Scaling is applied based on matrix 1-norms.
  - `verb` is the level of verbosity. It can take three values:
    - `0`: no verbosity
    - `1`: some verbosity
    - `2`: full verbosity

    Full verbosity can have a large performance impact, as it computes certain expensive quantities each iteration so as to provide more information for troubleshooting.
  - `check_singular` is whether to check if the QEP is numerically singular or not. The default value is `false`, because it is a potentially expensive test that may fail to detect nonsingularity.

### Returns
  - `λ`: array of approximate eigenvalues.
  - `V`: matrix of approximate eigenvectors (as columns).
  - `ρ`: array of associated backward error residuals (as defined above).

### Examples
Basic usage example:
```julia
M = rand(ComplexF64, 1000, 1000) #small dense matrices just for simple demonstration,
D = rand(ComplexF64, 1000, 1000) #not the standard use-case
K = rand(ComplexF64, 1000, 1000)

#we ask for 30 eigenpairs with residuals below 10^-10
#this will find those near the origin
λ, V, ρ = quadEigRBTOAR(M, D, K, 30, 1e-10)

#plot results, colouring eigenvalues according to residual
using Plots
scatter(λ, marker_z=log.(ρ), c=:rainbow2, xlabel="Re(λ)", ylabel="Im(λ)", legend=false)
```
To find 14 eigenvalues close to $69-420i$:
```julia
λ, V, ρ = quadEigRBTOAR(M, D, K, 14, σ=69-420im)
```
To find 40 largest-magnitude eigenvalues to a low accuracy:
```julia
λ, V, ρ = quadEigRBTOAR(M, D, K, 40, 1e-6, inv=false) CHANGE INV NAME
```
To print information during execution:
```julia

```
To use a larger block size:
```julia

```
**Don't** do this:
```julia
#what could go wrong? :3

```

Still need to cite sources, change name of inv, sort out checksingular(), add verbosity, verify examples, test code on more problems and refine defaults if necessary, check stupid edge-cases
