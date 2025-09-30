# Restarted Block TOAR
A pure Julia implementation of a restarted (see [2, Section 5]) block version of the TOAR [1] algorithm for the solution of large sparse quadratic eigenvalue problems. Currently a WIP, some features are missing and it is largely untested.

## Installation
Download `BTOAR.jl` to your Julia working directory and type
```julia-repl
julia> include("BTOAR.jl")
```
You may have to install the following dependencies:
  - `LowRankApprox.jl`
  - `QuadEig.jl`

## Usage
```julia
    quadEigRBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix; req::Int=100, tol::Float64=1e-10, kℓ_max::Int, ℓ::Int, step::Int=10, σ::Union{Float64,ComplexF64}=0.0+0.0im, which::Symbol=:SM, keep::Function=every, dtol::Float64=1e-10, rrv::Int=0, flvd::Bool=true, verb::Int=0, check_singular::Bool=false, give_up::Int=10)
```
The function builds a second-order Krylov subspace, restarting when necessary, until the desired number of eigenpairs can be recovered with an acceptable backward error.

### Arguments
  - `M`, `D` and `K` are the mass, damping and stiffness matrices respectively from the QEP $(\lambda^2M+\lambda D+K)x=0$. These should typically be given in `SparseMatrixCSC` format, although any `AbstractMatrix` type is allowed. Strictly speaking, these are the only three required arguments, although typically you would want to specify at least `req` and `tol`.
  - `req` is the number of eigenpairs required (default: 100). This should be set smaller than `kℓ_max` in practically all circumstances (perhaps *much* smaller) or else the algorithm may stagnate, that is, it may get stuck in a loop of expanding the subspace and restarting without ever reaching `req` eigenpairs.
  - `tol` is the tolerance for the backward error residual in computed eigenpairs. By default, it is $10^{-10}$. It is not reccommended to set this smaller than $10^{-13}$. Some QEPs permit finding more accurate eigenpairs than others; if you are not finding any "good" eigenpairs with the default `tol`, it may be that you have a difficult QEP and need to set `tol` lower. Backward error residuals of a computed eigenpair $(\tilde{\lambda},\tilde{x})$ are defined as

$$\rho=\frac{\Vert(\tilde{\lambda}^2M+\tilde{\lambda}D+K)\tilde{x}\Vert_2}{|\tilde{\lambda}|^2\Vert M\Vert_1+|\tilde{\lambda}|\Vert D\Vert_1+\Vert K\Vert_1}.$$

  - `kℓ_max` is the maximum subspace size the algorithm is permitted to build. By default it is 300. This should be larger than `req`, perhaps much larger for a difficult QEP. However, if memory requirements are a problem, you may have to lower `kℓ_max`. Setting this argument to a very high value has the effect of reducing the algorithm to non-restarted block TOAR.
  - `ℓ` is the block size/width. This should normally be set to at most 5. The default value, 1, reduces the algorithm to the standard non-block restarted TOAR algorithm. For larger block sizes, it might be a good idea to set `step` lower.
  - `step` is an optional argument specifying the minimum number of blocks to be added to the second-order Krylov subspace between each check for convergence of sufficiently many eigenpairs. Setting this too high could waste time building a larger subspace than neccessary; setting it too low could waste time solving the reduced-order QEP too many times when little progress is made enriching the subspace. `step` defaults to 10 if unspecified, which is probably unsuitable for larger block sizes.
  - `σ` is the shift point in $\mathbb{C}$ for the shift-and-invert spectral transformation. This should normally be set to a value well within the domain of interest, or convergence of wanted eigenvalues may be slow. By default there is no shift, that is, `σ` is set to 0.
  - `which` specifies which eigenvalues to target. The default is `:SM`, which targets eigenvalues closest to `σ`. The other accepted value is `:LM`, which targets eigenvalues furthest from `σ`.
  - `keep` is the function that specifies the domain of interest. This function should accept a single positional argument of type `ComplexF64` and return `true` if it is in the domain of interest and false otherwise. By default, the domain of interest is taken to be the full complex plane (so `keep` always returns `true` for all inputs).
  - `dtol` is an internal numerical tolerance for breakdown/deflation detection (default $10^{-10}$). Normally this should not be changed, except in the case of badly behaved QEPs if you know what you're doing.
  - `rrv` is the number of inverse power iterations to apply in the Ritz vector refinement. The default is 0, equivalent to doing nothing. **Currently not implemented, has no effect.**
  - `flvd` is whether to use Fan, Lin & Van Dooren [3] scaling on the QEP. Scaling is applied based on matrix 1-norms. By default, scaling is always applied, as it can sometimes help a lot and is quite cheap.
  - `verb` is the level of verbosity. It can take three values:
    - `0`: no verbosity
    - `1`: some verbosity
    - `2`: full verbosity

    Full verbosity can have a large performance impact, as it computes certain expensive quantities each iteration so as to provide more information for troubleshooting. By default, `verb` is set to 0, so the function will "not speak unless asked to".
  - `check_singular` is whether to check if the QEP is numerically singular or not. The default value is `false`, because it is a potentially expensive test that may fail to detect nonsingularity.
  - `give_up` specifies how many restarts to allow before terminating the algorithm without finishing, i.e. "giving up". When this happens, the algorithm will still return whatever it was able to find and raise a warning. The default is to allow (a reasonably generous) 10 restarts. If the algorithm does not succeed within 10 restarts, it is unlikely it ever will without changes to other parameters (like `req`, `tol`, `kℓ_max`) especially for low values of `tol`.

### Returns
  - `λ`: array of approximate eigenvalues.
  - `V`: matrix of approximate eigenvectors (as columns).
  - `ρ`: array of associated backward error residuals (as defined above).

### Examples
Basic usage example:
```julia
#add diagonal to make singularity unlikely
M = sprand(ComplexF64, 10000, 10000, 1e-3) + Diagonal(randn(ComplexF64, 10000))
D = sprand(ComplexF64, 10000, 10000, 1e-3) + Diagonal(randn(ComplexF64, 10000))
K = sprand(ComplexF64, 10000, 10000, 1e-3) + Diagonal(randn(ComplexF64, 10000))

#we ask for 30 eigenpairs with residuals below 10^-10
#this will find those near the origin
λ, V, ρ = quadEigRBTOAR(M, D, K, req=30, tol=1e-10)

#plot results, colouring eigenvalues according to residual
using Plots
display(scatter(λ, marker_z=log.(ρ), c=:rainbow2, xlabel="Re(λ)", ylabel="Im(λ)", legend=false))
```
To find 14 eigenvalues close to $69-420i$ (if there are any):
```julia
λ, V, ρ = quadEigRBTOAR(M, D, K, req=14, σ=69.0-420.0im)
```
To find 40 largest-magnitude eigenvalues to a low accuracy:
```julia
λ, V, ρ = quadEigRBTOAR(M, D, K, req=40, tol=1e-6, which=:LM)
```
To print information during execution:
```julia
λ, V, ρ = quadEigRBTOAR(M, D, K, tol=1e-10, verb=1)
```
To use a larger block size:
```julia
λ, V, ρ = quadEigRBTOAR(M, D, K, tol=1e-10, ℓ=4, step=3)
```
To find eigenvalues in the top-right quadrant of $\mathbb{C}$:
```julia
doi(λ) = (real.(λ) > 0) && (imag.(λ) > 0)
λ, V, ρ = quadEigRBTOAR(M, D, K, tol=1e-10, σ=0.07+0.07im, keep=doi)
```
**Don't** do this:
```julia
#what could go wrong? :3
doi(λ) = abs(λ) .> 1e3
λ, V, ρ = quadEigRBTOAR(M, D, K, tol=1e-16, req=120, kℓ_max=100, ℓ=20, dtol=1e-3, keep=doi)
```

## References
[1] D. Lu, Y. Su, and Z. Bai. "Stability Analysis of the Two-Level Orthogonal Arnoldi Procedure". In: SIAM J. Matrix Anal. Appl. 37 (2016), pp. 195–214. URL: https://epubs.siam.org/doi/10.1137/151005142

[2] C. Campos and J. E. Roman. "Restarted Q-Arnoldi-Type Methods Exploiting Symmetry in Quadratic Eigenvalue Problems". In: BIT Numerical Mathematics 56 (2016), pp. 1213–1236. URL: https://link.springer.com/article/10.1007/s10543-016-0601-5

[3] H. Fan and W. Lin and P. Van Dooren. "Normwise Scaling of Second Order Polynomial Matrices". In: SIAM J. Matrix Anal. Appl. 26 (2004), pp. 252-256. URL: https://epubs.siam.org/doi/abs/10.1137/S0895479803434914
