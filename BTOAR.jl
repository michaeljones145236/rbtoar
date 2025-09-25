using LinearAlgebra,LowRankApprox,QuadEig,Printf

function fmean(v::Vector,f::Function,f⁻¹::Function)
    return f⁻¹(sum(f.(v))/size(v)[1])
end

function b(r::Union{Int,UnitRange{Int}},ℓ::Int) #function to help out with long block indices with constant block size
    if typeof(r) <: Int
        return (r-1)ℓ+1:r*ℓ
    else
        return (minimum(r)-1)ℓ+1:maximum(r)*ℓ
    end
end

#this function doesn't actually compute a rrqr factorisation because R will not necessarily be upper trapezoidal, but we don't really care
function rrqr(A::Matrix,tol::Float64)
    Q,R,p = pqr(A,rtol=tol) #use the so-called `partial QR factorisation'
    return Q,R*ColumnPermutation(p)' #return Q and (unfortunately non-upper trapezoidal) R
end

every(λ) = true #this function is just here as the default for the keep function in quadEigRBTOAR

deflate_distant(λ) = [(sum(abs(λ[i]) .< abs.(λ)) < 2*size(λ,1)/3) for i in 1:size(λ,1)] #we actually want to deflate those closest to 0 in the shifted-inverted QEP

"""
    BTOAR(M⁻¹::Function, D::Function, K::Function, R::Matrix{Complex{Float64}}, k::Int, deftol::Float64=NaN, verb::Int=0)
    
Compute an orthogonal basis for the `k`th second-order Krylov subspace Gₖ(M⁻¹D,M⁻¹K;R).

# Arguments
 -`M⁻¹::Function`: function that provides left-multiplication of `n`×`ℓ` matrices by the inverse of M.\n
 -`D::Function`: function that provides left-multiplication of `n`×`ℓ` matrices by the matrix D.\n
 -`K::Function`: function that provides left-multiplication of `n`×`ℓ` matrices by the matrix K.\n
 -`R::Matrix{Complex{Float64}}`: starting `n`×`ℓ` block vector R.\n
 -`k::Int`: degree of the second-order Krylov subspace Gₖ(A,B;R) required.\n
 -`deftol::Float64`: internal numerical tolerance for deflation detection, defaults to `1e-10`.\n
 -`verb::Int`: level of verbosity. 0: no verbosity, 1: some verbosity, 2: full verbosity. `verb=2` has a large performance impact.\n

# Returns
 -`Qₖ::Matrix{Complex{Float64}}`: orthonormal basis for second-order Krylov subspace Gₖ(A,B;R).\n
 -`Uₖ::Matrix{Complex{Float64}}`: orthonormal matrix such that `Vₖ`=(`I₂`⊗`Qₖ`)`Uₖ`, where `Vₖ` is an orthonormal basis for the first-order Krylov subspace Kₖ([A B;I 0];[`R`;0]).\n
 -`Hₖ::Matrix{Complex{Float64}}`: block-upper Hessenberg matrix such that [A B;I 0]`Vₖ₋₁`=`VₖHₖ`, where `Vₖ` is as defined above.\n
 -`ρ₁::Float64`: 1-norm orthonormality residual ‖`QₖᴴQₖ`-I‖₁.\n
 -`ρ₂::Float64`: 1-norm orthonormality residual ‖`UₖᴴUₖ`-I‖₁.\n
 -`ρ₃::Float64`: 1-norm TOAR relation residual ‖[A B;I 0]*(`I₂`⊗`Qₖ₋₁`)`Uₖ₋₁` - (`I₂`⊗`Qₖ`)`UₖHₖ`‖₁ / ‖(`I₂`⊗`Qₖ`)`UₖHₖ`‖₁.\n
"""
function BTOAR(M⁻¹::Function,D::Function,K::Function,R::Matrix,k::Int;deftol::Float64=NaN,verb::Int=0)
    n,ℓ = size(R) #n and ℓ are not function arguments, but are taken implicitly from R

    if deftol ≡ NaN #if deftol isn't set, we give it a reasonable default
        deftol = ℓ*eps(Float64) #machine epsilon times the smallest dimension of the matrix is a standard numerical rank tolerance
    end
    if deftol < eps(Float64) #warning about setting deflation tolerance too low
        @warn "deftol should not be set lower than ϵ≈2.22×10⁻¹⁶ (deftol=$deftol)\nSetting deftol to ϵ"
        deftol = eps(Float64)
    end
    if deftol ≥ 1.0 #not sure what the highest reasonable deftol would be
        @warn "deftol way too large (deftol=$deftol)\nSetting deftol to ϵℓ=$(ℓ*eps(Float64))"
        deftol = ℓ*eps(Float64)
    end
    if k*ℓ > n #warning about setting k too large
        @warn "kℓ greater than n, expect deflation (kℓ = $(k*ℓ), n = $n)"
    end
    if verb > 2
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 2"
        verb = 2
    end
    if verb < 0
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 0"
        verb = 0
    end
    
    m = zeros(Int,k) #we preallocate the array of ranks of Rⱼ-QⱼSⱼ (note that m[1] is m₀, not m₁)
    m[1] = rank(R[1:ℓ,:],rtol=deftol) #to save time, we start by taking the rank of a small full-width submatrix of R
    if m[1] < ℓ #if the small full-width submatrix was rank deficient (unlikely most of the time) 
        m[1] = rank(R,rtol=deftol) #test the rank of the full block vector R
        if m[1] < ℓ #if R is rank-deficient
            R = Matrix(rrqr(R)[1]) #we reduce R to an orthonormal matrix with the same span, using the RRQR factorisation
            ℓ = m[1] #the new ℓ after R has been reduced
            if verb > 0
                print("Starting vector R rank deficient, reducing ℓ to $ℓ\n")
            end
        end
    end
    Q,R = qr(R) #standard QR factorisation because R now must be full-rank
    Qⱼ = Matrix{ComplexF64}(Q) #qr() outputs Q in a special form without explicitly forming the matrix, so we have to tell it to
    Uⱼ = [I(ℓ);zeros(Complex{Float64},ℓ,ℓ)] #because V₁ = [Q₁;0]
    Hₖ = zeros(Complex{Float64},k*ℓ,(k-1)*ℓ) #we can preallocate H because its final size is known now
    
    for j in 1:k-1 #main for loop
        Rⱼ = M⁻¹(D(Qⱼ*-Uⱼ[1:sum(m),(j-1)ℓ+1:j*ℓ]) + K(Qⱼ*-Uⱼ[sum(m)+1:2sum(m),(j-1)ℓ+1:j*ℓ])) #take next Rⱼ block vector
        Sⱼ = zeros(Complex{Float64},sum(m[1:j]),ℓ) #preallocate Sⱼ
        for i in 1:sum(m[1:j]) #doing it this way seems to greatly reduce error for Ansys QEPs
            Sⱼ[i:i,:] = Qⱼ[:,i:i]'*Rⱼ #coefficients of components of Rⱼ parallel to columns of Qⱼ
            Rⱼ -= Qⱼ[:,i:i]*Sⱼ[i:i,:] #subtract off parts of Rⱼ that are parallel to columns of Qⱼ
        end
        for i in 1:sum(m[1:j]) #partial reorthogonalisation
            Sᵣₑₛ = Qⱼ[:,i:i]'*Rⱼ #not the full residual (block) vector
            Rⱼ -= Qⱼ[:,i:i]*Sᵣₑₛ #correct Rⱼ
            Sⱼ[i:i,:] += Sᵣₑₛ #correct Sⱼ
        end
        Qʰ,Rʰ = rrqr(Rⱼ,deftol) #MIGHT WANT TO DO CHEAP RANK TEST HERE FOR EFFICIENCY DEPENDING ON PERFORMANCE OF RRQR
        m[j+1] = size(Qʰ,2) #record the rank of Rⱼ-Qⱼ*Sⱼ
        if m[j+1] < ℓ #if we have deflation
            if verb == 1 #if we have medium verbosity
                print("🟨j=$j,mⱼ=$(m[j+1])\n") #tell the user
            elseif verb == 2 #if we have high verbosity
                print("🟨Deflation at j=$j (mⱼ=$(m[j+1])) ") #tell the user more
                #here, we want to estimate the maximum value of deftol that would have caused the deflation to pass
                tdt = deftol #trial deftol starts as deftol
                for i in 1:20 #don't need a very fine approximation
                    if size(rrqr(Rⱼ,tdt)[1],2) == m[j+1] #if too high
                        tdt *= (eps(Float64) / deftol)^(2.0^-i)
                    else #if too low
                        tdt *= (deftol / eps(Float64))^(2.0^-i)
                    end
                end
                print("(mdt=$tdt)\n")
            end
        elseif verb == 2 #if verbosity is high and there is no deflation
            print("🟩No deflation at j=$j ") #we inform the user
            #here, we want to estimate the minimum value of deftol that would have caught a "deflation"
            tdt = deftol #trial deftol starts as deftol
            for i in 1:20 #don't need a very fine approximation
                if size(rrqr(Rⱼ,tdt)[1],2) < m[j+1] #if too high
                    tdt *= (deftol / 1.0)^(2.0^-i)
                else #if too high
                    tdt *= (1.0 / deftol)^(2.0^-i)
                end
            end
            print("(mdt=$tdt)\n")
        elseif verb == 1 #finally, if verbosity is medium and there is no deflation
            print("🟩") #just a green square with no newline character
        end
        Qⱼ = [Qⱼ Qʰ] #append the newly computed columns to Qⱼ
        Uʰ = Uⱼ[1:sum(m[1:j]),(j-1)ℓ+1:j*ℓ] #we copy this because we want to modify it in the orthogonalisation without modifying Uⱼ
        for i = 1:j #second-level orthogonalisation
            Hₖ[(i-1)ℓ+1:i*ℓ,(j-1)ℓ+1:j*ℓ] = Uⱼ[1:sum(m[1:j]),(i-1)ℓ+1:i*ℓ]'*Sⱼ + Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),(i-1)ℓ+1:i*ℓ]'*Uʰ #fill in new block column of Hₖ
            Sⱼ -= Uⱼ[1:sum(m[1:j]),(i-1)ℓ+1:i*ℓ]*Hₖ[(i-1)ℓ+1:i*ℓ,(j-1)ℓ+1:j*ℓ] #orthogonalise Sⱼ against Uⱼ,₁
            Uʰ -= Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),(i-1)ℓ+1:i*ℓ]*Hₖ[(i-1)ℓ+1:i*ℓ,(j-1)ℓ+1:j*ℓ] #and Uʰ against Uⱼ,₂
        end
        for i = 1:j #second-level reorthogonalisation
            Hᵣₑₛ = Uⱼ[1:sum(m[1:j]),(i-1)ℓ+1:i*ℓ]'*Sⱼ + Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),(i-1)ℓ+1:i*ℓ]'*Uʰ #not the full residual
            Sⱼ -= Uⱼ[1:sum(m[1:j]),(i-1)ℓ+1:i*ℓ]*Hᵣₑₛ #reorthogonalise Sⱼ against Uⱼ,₁
            Uʰ -= Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),(i-1)ℓ+1:i*ℓ]*Hᵣₑₛ #and Uʰ against Uⱼ,₂
            Hₖ[(i-1)ℓ+1:i*ℓ,(j-1)ℓ+1:j*ℓ] += Hᵣₑₛ #correct block column of Hₖ
        end
        Qᵗ,Rᵗ = qr([Sⱼ;Rʰ;Uʰ]) #standard QR factorisation, not RRQR
        Qᵗ = Matrix(Qᵗ) #we must force explicit formation of Qᵗ
        if rank(Rᵗ,rtol=deftol) < ℓ #this means (at least partial) breakdown of the concurrent block Arnoldi procedure
            if verb == 2
                rk = rank(Rᵗ,rtol=deftol)
                print("🟥Breakdown at j=$j (rank=$rk)\n")
            elseif verb == 1
                print("🟥j=$j\n")
            end
            @warn "breakdown at j=$j, returning only Qⱼ" #warn the user, to prevent silent failure from causing problems
            return Qⱼ,nothing,nothing,nothing,nothing,nothing #return Qⱼ because it might still be useful (probably not)
        end
        Hₖ[j*ℓ+1:(j+1)ℓ,(j-1)ℓ+1:j*ℓ] = Rᵗ #fill in bottom-right block entry of Hₖ
        Uⱼ = [[[Uⱼ[1:sum(m[1:j]),:];zeros(m[j+1],j*ℓ);Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),:]] Qᵗ];zeros(m[j+1],(j+1)ℓ)] #form new columns and rows of Uⱼ
        if verb == 2
            print("TOAR relation residual: ")
            display(opnorm([-(M⁻¹(D(Qⱼ[:,1:sum(m[1:j])]*Uⱼ[1:sum(m[1:j]),1:j*ℓ]) + K(Qⱼ[:,1:sum(m[1:j])]*Uⱼ[sum(m)+1:sum(m)+sum(m[1:j]),1:j*ℓ])));Qⱼ[:,1:sum(m[1:j])]*Uⱼ[1:sum(m[1:j]),1:j*ℓ]] - [Qⱼ*Uⱼ[1:sum(m),:];Qⱼ*Uⱼ[sum(m)+1:2sum(m),:]]*Hₖ[1:(j+1)ℓ,1:j*ℓ],1) / opnorm([Qⱼ*Uⱼ[1:sum(m),:];Qⱼ*Uⱼ[sum(m)+1:2sum(m),:]]*Hₖ[1:(j+1)ℓ,1:j*ℓ],1)) #this was soul-crushing to debug
            print("Qⱼ residual: ")
            display(opnorm(Qⱼ'*Qⱼ-I,1))
            print("Uⱼ residual: ")
            display(opnorm(Uⱼ'*Uⱼ-I,1))
            print("\n")
        end
    end
    j = k-1 #rather than going through the expression for ρ₃ and changing every j to k-1, I just set this and copy-paste it
    ρ₁ = opnorm(Qⱼ'*Qⱼ-I,1) # }
    ρ₂ = opnorm(Uⱼ'*Uⱼ-I,1) # } - these two residuals are easier to type out than ρ₃
    ρ₃ = opnorm([-(M⁻¹(D(Qⱼ[:,1:sum(m[1:j])]*Uⱼ[1:sum(m[1:j]),1:j*ℓ]) + K(Qⱼ[:,1:sum(m[1:j])]*Uⱼ[sum(m)+1:sum(m)+sum(m[1:j]),1:j*ℓ])));Qⱼ[:,1:sum(m[1:j])]*Uⱼ[1:sum(m[1:j]),1:j*ℓ]] - [Qⱼ*Uⱼ[1:sum(m),:];Qⱼ*Uⱼ[sum(m)+1:2sum(m),:]]*Hₖ[1:(j+1)ℓ,1:j*ℓ],1) / opnorm([Qⱼ*Uⱼ[1:sum(m),:];Qⱼ*Uⱼ[sum(m)+1:2sum(m),:]]*Hₖ[1:(j+1)ℓ,1:j*ℓ],1) #this was soul-crushing to debug
    if verb == 1
        print("🟦\nTOAR relation residual: ")
        display(ρ₃)
        print("Qⱼ residual: ")
        display(opnorm(Qⱼ'*Qⱼ-I,1))
        print("Uⱼ residual: ")
        display(opnorm(Uⱼ'*Uⱼ-I,1))
        print("\n")
    elseif verb == 2
        print("🟦Terminated successfully.\n\n")
    end
    return Qⱼ,Uⱼ,Hₖ,ρ₁,ρ₂,ρ₃
end

"""
    quadEigBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix, k::Int, ℓ::Int; σ::Union{Float64,Complex{Float64}}=0.0+0.0im, inv::Bool=true, dtol::Float64=1e-10, rrv::Int=0, flvd::Bool=true, verb::Int=0, check_singular::Bool=false)
    
Compute some eigenpairs of the QEP `(λ²M + λD + K)x=0` using the block TOAR algorithm.

# Arguments:
 -`M::AbstractMatrix`: mass matrix of the QEP.\n
 -`D::AbstractMatrix`: damping matrix of the QEP.\n
 -`K::AbstractMatrix`: stiffness matrix of the QEP.\n
 -`k::Int`: number of block TOAR iterations.\n
 -`ℓ::Int`: block size/width.\n
 -`σ::Union{Float64,Complex{Float64}}`: spectral transformation shift to use (default `0`).\n
 -`inv::Bool`: whether to use spectral inversion (default yes).\n
 -`dtol::Float64`: optional internal numerical deflation/breakdown tolerance for BTOAR.\n
 -`rrv::Int`: the number of inverse power iterations to use in Ritz vector refinement (default `0`).\n
 -`flvd::Bool`: whether to use Fan, Lin and Van Dooren scaling on the QEP (default yes).\n
 -`verb::Int`: level of verbosity. 0: no verbosity, 1: some verbosity, 2: full verbosity. `verb=2` has a large performance impact.\n
 -`check_singular::Bool`: whether to test for numerical singularity and warn if it is found. Default `false`, setting to `true` can have a significant performance impact.\n

# Returns:
 -`λ::Vector{ComplexF64}`: Array of Ritz values (approximate eigenvalues).\n
 -`V::Matrix{ComplexF64}`: Array of Ritz vectors (approximate eigenvectors).\n
 -`ρ::Vector{ComplexF64}`: Array of backward error residuals for the eigenpairs.\n
"""
function quadEigBTOAR(M::AbstractMatrix,D::AbstractMatrix,K::AbstractMatrix,k::Int,ℓ::Int;σ::Union{Float64,Complex{Float64}}=0.0+0.0im,inv::Bool=true,dtol::Float64=1e-10,rrv::Int=0,flvd::Bool=true,verb::Int=0,check_singular::Bool=false)
    n = size(M,1) #take n implicitly
    if false in (n .== [size(M,2);size(D,1);size(D,2);size(K,1);size(K,2)]) #M, D and K must all be n×n
        error("M, D and K must all be n×n")
    end
    
    if verb > 2
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 2"
        verb = 2
    end
    if verb < 0
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 0"
        verb = 0
    end
    
    if check_singular
        if (cond(M,1) > 1e10) && (cond(K,1) > 1e10) #only if both of these two are singular can the whole QEP be
            test_λ = [1e-9,1e-6,1e-3,1.0,1e3,1e6,1e9,-1e-9,-1e-6,-1e-3,-1.0,-1e3,-1e6,-1e9] #only test real values, this is not exhaustive
            could_be_singular = true
            for i in test_λ
                if cond(i^2*M+i*D+K) < 1e10
                    could_be_singular = false
                    break
                end
            end
        end
        if could_be_singular
            @warn "QEP could be numerically singular."
        elseif verb == 2
            print("QEP is not numerically singular 👍\n\n")
        end
    end

    Mₙₒᵣₘ = opnorm(M,1) #precompute matrix norms for efficiency
    Dₙₒᵣₘ = opnorm(D,1)
    Kₙₒᵣₘ = opnorm(K,1)
    badly_scaled_QEP = maximum(abs.(log.(10,[Mₙₒᵣₘ/Dₙₒᵣₘ,Mₙₒᵣₘ/Kₙₒᵣₘ,Dₙₒᵣₘ/Kₙₒᵣₘ]))) > 8 #the QEP is considered badly scaled if the ratio of any of the 1-norms of its coefficient matrices is more than 10⁸
    #Fan, Lin & Van Dooren scaling (2004)
    γ = flvd ? sqrt(Kₙₒᵣₘ/Mₙₒᵣₘ) : 1.0 #eigenvalue scaling
    δ = flvd ? 2/(Kₙₒᵣₘ+γ*Dₙₒᵣₘ) : 1.0 #uniform scaling
    if flvd
        still_badly_scaled = maximum(abs.(log.(10,[γ*Mₙₒᵣₘ/Dₙₒᵣₘ,γ^2*Mₙₒᵣₘ/Kₙₒᵣₘ,γ*Dₙₒᵣₘ/Kₙₒᵣₘ]))) > 8 #cancelled factors of γ and excluded δ since it makes no difference
        if still_badly_scaled
            @warn "QEP still badly scaled after Fan, Lin & Van Dooren scaling: γ²δ‖M‖₁=$(γ^2*δ*Mₙₒᵣₘ), γδ‖D‖₁=$(γ*δ*Dₙₒᵣₘ), δ‖K‖₁=$(δ*Kₙₒᵣₘ)."
        end
    else
        if badly_scaled_QEP
            @warn "QEP is badly scaled: ‖M‖₁=$Mₙₒᵣₘ, ‖D‖₁=$Dₙₒᵣₘ, ‖K‖₁=$Kₙₒᵣₘ. Consider setting flvd=true."
        end
    end
    if verb > 0 #some or all verbosity
        print("== SCALING INFORMATION ==\n")
        if flvd
            print("    Fan, Lin & Van Dooren scaling applied with γ=$γ, δ=$δ.\n    Pre-scaling matrix norms:\n        ‖M‖₁=$Mₙₒᵣₘ\n        ‖D‖₁=$Dₙₒᵣₘ\n        ‖K‖₁=$Kₙₒᵣₘ\n    Scaled matrix norms:\n        γ²δ‖M‖₁=$(γ^2*δ*Mₙₒᵣₘ)\n        γδ‖D‖₁=$(γ*δ*Dₙₒᵣₘ)\n        δ‖K‖₁=$(δ*Kₙₒᵣₘ)\n\n")
        else
            print("    No scaling applied.\n    Matrix norms:\n        ‖M‖₁=$Mₙₒᵣₘ\n        ‖D‖₁=$Dₙₒᵣₘ\n        ‖K‖₁=$Kₙₒᵣₘ\n\n")
        end
    end
    
    Kₛ = δ*K + σ*δ*D + σ^2*δ*M #scaled and shifted matrices (I multiplied out the cancelling factors of γ)
    Dₛ = δ*γ*D + 2σ*δ*γ*M
    Mₛ = δ*γ^2*M
    
    Mᴸᵁ = lu(inv ? Kₛ : Mₛ) #factorise M for fast linear solves
    M⁻¹(x) = Mᴸᵁ\x #use factorised version of M
    M¹(x) = (inv ? Kₛ : Mₛ)*x
    D¹(x) = Dₛ*x #we call the functions D¹ and K¹ because the names D and K are taken
    K¹(x) = (inv ? Mₛ : Kₛ)*x
    
    if verb > 0
        print("== BTOAR ALGORITHM ==\n\n")
    end
    Qₖ,_,_,_,_,_ = BTOAR(M⁻¹,D¹,K¹,rand(ComplexF64,n,ℓ),k,deftol=dtol,verb=verb) #run the BTOAR algorithm
    m = size(Qₖ,2) #there might have been deflations
    
    if false #rrv (temporary bodge)
        MQₖ = M¹(Qₖ) #save intermediate step for refined ritz vectors
        DQₖ = D¹(Qₖ)
        KQₖ = K¹(Qₖ)
        Mₘ = Qₖ'*MQₖ #reduced order QEP
        Dₘ = Qₖ'*DQₖ
        Kₘ = Qₖ'*KQₖ
    else
        Mₘ = Qₖ'*M¹(Qₖ) #reduced order QEP
        Dₘ = Qₖ'*D¹(Qₖ)
        Kₘ = Qₖ'*K¹(Qₖ)
    end
    LP = linearize(Kₘ,Dₘ,Mₘ) #linearize the QEP
    Z = zeros(ComplexF64,m,2m) #preallocate space for reduced-order eigenvectors
    λ,V = eigen(LP.A,LP.B) #generalised eigenproblem solver
    if false #rrv (temporary bodge)
        MQₖᵀMQₖ = MQₖ'MQₖ #pre-form matrices for faster formation of PQ
        MQₖᵀDQₖ = MQₖ'DQₖ #these multiplications are (kℓ×n)(n×kℓ)
        MQₖᵀKQₖ = MQₖ'KQₖ
        DQₖᵀMQₖ = DQₖ'MQₖ
        DQₖᵀDQₖ = DQₖ'DQₖ
        DQₖᵀKQₖ = DQₖ'KQₖ #yeah, there are a lot of them but it's most efficient way if m >> 9
        KQₖᵀMQₖ = KQₖ'MQₖ
        KQₖᵀDQₖ = KQₖ'DQₖ
        KQₖᵀKQₖ = KQₖ'KQₖ
        for i in 1:2m
            #additions here are only of kℓ×kℓ matrices (cheap)
            PQᵀPQ = λ[i]'^2*(λ[i]^2*MQₖᵀMQₖ + λ[i]*MQₖᵀDQₖ + MQₖᵀKQₖ) + λ[i]'*(λ[i]^2*DQₖᵀMQₖ + λ[i]*DQₖᵀDQₖ + DQₖᵀKQₖ) + λ[i]^2*KQₖᵀMQₖ + λ[i]*KQₖᵀDQₖ + KQₖᵀKQₖ #matrix to find least dominant right singular vector of
            Z[:,i] = arpack ? eigs(PQᵀPQ,nev=1,which=:LM,ritzvec=true,v0=V[1:m,i],sigma=0.0,tol=1e-50)[2] : svd(λ[i]^2*MQₖ + λ[i]*DQₖ + KQₖ).V[:,m]
        end
    else
        for i in 1:2m #extract quadratic eigenvectors
            Z[:,i] = (abs(λ[i]) > 1) ? V[1:m,i] : V[m+1:2m,i] #this should be more stable I think
            Z[:,i] /= norm(Z[:,i]) #normalise
        end
    end
    
    λ .^= (inv ? -1 : 1) #if QEP was inverted, uninvert eigenvalues
    λ .+= σ/γ #unshift eigenvalues
    λ .*= γ #unscale eigenvalues
    
    X = Qₖ*Z #transform eigenvectors back to full-sized space (they stay normalised)
    
    ρ = zeros(2m) #preallocate space for residuals
    for i in 1:2m #calculate backward error residuals
        ρ[i] = norm(λ[i]^2*M*X[:,i]+λ[i]*D*X[:,i]+K*X[:,i]) / (abs2(λ[i])*Mₙₒᵣₘ+abs(λ[i])*Dₙₒᵣₘ+Kₙₒᵣₘ)
    end
    if verb == 2 #if maximum verbosity, we say a bit about the quality of results
        print("== RESIDUAL INFORMATION ==\n    num. residuals below 10⁻¹⁵: $(sum(ρ .≤ 1e-15))\n    num. residuals below 10⁻¹²: $(sum(ρ .≤ 1e-12))\n    num. residuals below 10⁻⁹: $(sum(ρ .≤ 1e-9))\n    num. residuals below 10⁻⁶: $(sum(ρ .≤ 1e-6))\n\n")
    end
    return λ,X,ρ #return eigenvalues, eigenvectors, residuals
end

"""
    restartBTOAR(Q□::Matrix, U□::Matrix, H□::Matrix, keep::Function, verb::Int=0)

Purge the second-order Krylov basis `Q□` of some of its Ritz vector dimensions, while keeping selected ones.

# Arguments
 -`Q□::Matrix`: second-order Krylov basis from BTOAR to be purged.\n
 -`U□::Matrix`: auxiliary matrix from BTOAR.\n
 -`H□::Matrix`: auxiliary matrix from BTOAR.\n
 -`keep::Function`: function `keep(λ)` that returns `true` if an eigenvalue should be locked and `false` if it should be purged.\n
 -`verb::Int`: verbosity level. 0: no verbosity, 1: some verbosity, 2: full verbosity.\n

# Returns
 -`Qᵣ::Matrix`: reduced-dimension second-order Krylov basis.\n
 -`Uₚ₊₁::Matrix`: reduced auxiliary matrix.\n
 -`Hₚ₊₁::Matrix`: reduced auxiliary matrix.\n
"""
function restartBTOAR(Q□::Matrix,U□::Matrix,H□::Matrix,keep::Function,verb::Int=0)
    if verb > 2
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 2"
        verb = 2
    end
    if verb < 0
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 0"
        verb = 0
    end

    schurfact = schur(H□[1:size(H□,2),:]) #Schur factorisation of the top square
    #if !(false in keep(schurfact.values)) #if we don't need to get rid of anything, no need to do any work
    #    if verb == 2
    #        print("Deflated nothing.\n")
    #    end
    #    if verb > 0 #brief verbose output of number of discarded/retained dimensions
    #        print("Dimensions retained: $(size(Qᵣ,2))\n")
    #        print("Dimensions purged: 0\n\n")
    #    end
    #    return Q□,U□,H□
    #end
    ordschur!(schurfact,keep(schurfact.values)) #reorder Schur factorisation to put desired eigenvalues in top left
    
    #determine some constants implicitly
    p = sum(keep(schurfact.values))
    ℓ = size(H□,1)-size(H□,2)
    
    Hₚ₊₁ = [schurfact.T[1:p,1:p];H□[size(H□,2)+1:size(H□,1),:]*schurfact.Z[:,1:p]] #form new H
    
    U□⁽¹⁾ = U□[1:size(Q□,2),:] #deconstruct U□ for readability
    U□⁽²⁾ = U□[size(Q□,2)+1:2size(Q□,2),:]
    U,Σ,V = psvd([[U□⁽¹⁾[1:size(U□⁽¹⁾,1)-ℓ,1:size(U□,2)-ℓ]*schurfact.Z[:,1:p];zeros(ℓ,p)] U□⁽¹⁾[:,size(U□,2)-ℓ+1:size(U□,2)] [U□⁽²⁾[1:size(U□⁽²⁾,1)-ℓ,1:size(U□,2)-ℓ]*schurfact.Z[:,1:p];zeros(ℓ,p)] U□⁽²⁾[:,size(U□,2)-ℓ+1:size(U□,2)]],rank=p+2ℓ) #that took a while to type
    
    Uₚ₊₁ = kron(I(2),Diagonal(Σ))*[V'[:,1:p+ℓ];V'[:,p+ℓ+1:2(p+ℓ)]]
    
    Qᵣ = Q□*U #dominant flop cost
    
    if verb > 0 #brief verbose output of number of discarded/retained dimensions
        print("Dimensions retained: $(size(Qᵣ,2))\n")
        print("Dimensions purged: $(size(Q□,2)-size(Qᵣ,2))\n\n")
    end
    
    return Qᵣ,Uₚ₊₁,Hₚ₊₁
end

"""
    continueBTOAR(M⁻¹::Function, D::Function, K::Function, Qᵣ::Matrix{Complex{Float64}}, Uₚ₊₁::Matrix{Complex{Float64}}, Hₚ₊₁::Matrix{Complex{Float64}}, k::Int, ℓ::Int; deftol::Float64=1e-10)

Continue the BTOAR algorithm after a restart.

# Arguments
 -`M⁻¹::Function`: function that provides left-multiplication of `n`×`ℓ` matrices by the inverse of M.\n
 -`D::Function`: function that provides left-multiplication of `n`×`ℓ` matrices by the matrix D.\n
 -`K::Function`: function that provides left-multiplication of `n`×`ℓ` matrices by the matrix K.\n
 -`Qᵣ::Matrix{Complex{Float64}}`: locked second-order Krylov basis.\n
 -`Uₚ₊₁::Matrix{Complex{Float64}}`: locked auxiliary matrix.\n
 -`Hₚ₊₁::Matrix{Complex{Float64}}`: locked auxiliary matrix.\n
 -`k::Int`: number of new blocks to add to the subspace.\n
 -`ℓ::Int`: block size/width.\n
 -`deftol::Float64`: internal numerical tolerance for deflation detection, defaults to `1e-10`.\n
 -`verb::Int`: verbosity option. 0: no verbosity, 1: some verbosity, 2: full verbosity. Full verbosity has a large performance impact.\n

# Returns
 -`Q::Matrix`: second-order Krylov subspace basis.\n
 -`U::Matrix`: auxiliary matrix.\n
 -`H::Matrix`: auxiliary matrix.\n
"""
function continueBTOAR(M⁻¹::Function,D::Function,K::Function,Qᵣ::Matrix,Uₚ₊₁::Matrix,Hₚ₊₁::Matrix,k::Int,ℓ::Int;deftol::Float64=1e-10,verb::Int=0)
    n = size(Qᵣ,1) #take n implicitly
    
    if deftol < eps(Float64) #warning about setting deflation tolerance too low
        @warn "deftol should not be set lower than ϵ≈2.22×10⁻¹⁶ (deftol=$deftol)\nSetting deftol to ϵ"
        deftol = eps(Float64)
    end
    if deftol ≥ 1.0 #not sure what the highest reasonable deftol would be
        @warn "deftol way too large (deftol=$deftol)\nSetting deftol to ϵℓ=$(ℓ*eps(Float64))"
        deftol = ℓ*eps(Float64)
    end
    if k*ℓ+size(Qᵣ,2) > n #warning about setting k too large
        @warn "kℓ+r greater than n, expect deflation (kℓ+r = $(k*ℓ+size(Qᵣ,2)), n = $n)"
    end
    if verb > 2
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 2"
        verb = 2
    end
    if verb < 0
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 0"
        verb = 0
    end
    
    #no need to set up quantities, they already exist, we just need to make copies (inefficient, should really write an in-place function)
    Q = Matrix{ComplexF64}(Qᵣ)
    U = Matrix{ComplexF64}(Uₚ₊₁)
    H = Matrix{ComplexF64}(Hₚ₊₁)
    
    for j in 1:k #main for loop
        Rⱼ = M⁻¹(D(Q*-U[1:Int(size(U,1)/2),size(U,2)-ℓ+1:size(U,2)]) + K(Q*-U[Int(size(U,1)/2)+1:size(U,1),size(U,2)-ℓ+1:size(U,2)])) #take next Rⱼ block vector
        Sⱼ = zeros(Complex{Float64},size(Q,2),ℓ) #preallocate Sⱼ
        for i in 1:size(Q,2) #doing it this way seems to greatly reduce error for Ansys QEPs
            Sⱼ[i:i,:] = Q[:,i:i]'*Rⱼ #coefficients of components of Rⱼ parallel to columns of Qⱼ
            Rⱼ -= Q[:,i:i]*Sⱼ[i:i,:] #subtract off parts of Rⱼ that are parallel to columns of Qⱼ
        end
        for i in 1:size(Q,2) #partial reorthogonalisation
            Sᵣₑₛ = Q[:,i:i]'*Rⱼ #not the full residual (block) vector
            Rⱼ -= Q[:,i:i]*Sᵣₑₛ #correct Rⱼ
            Sⱼ[i:i,:] += Sᵣₑₛ #correct Sⱼ
        end
        Qʰ,Rʰ = rrqr(Rⱼ,deftol) #MIGHT WANT TO DO CHEAP RANK TEST HERE FOR EFFICIENCY DEPENDING ON PERFORMANCE OF RRQR
        m = size(Qʰ,2) #record the rank of Rⱼ-Qⱼ*Sⱼ
        if m < ℓ #if we have deflation
            if verb == 1 #if we have medium verbosity
                print("🟨j=$j,mⱼ=$m\n") #tell the user
            elseif verb == 2 #if we have high verbosity
                print("🟨Deflation at j=$j (mⱼ=$m) ") #tell the user more
                #here, we want to estimate the maximum value of deftol that would have caused the deflation to pass
                tdt = deftol #trial deftol starts as deftol
                for i in 1:20 #don't need a very fine approximation
                    if size(rrqr(Rⱼ,tdt)[1],2) == m #if too high
                        tdt *= (eps(Float64) / deftol)^(2.0^-i)
                    else #if too low
                        tdt *= (deftol / eps(Float64))^(2.0^-i)
                    end
                end
                print("(mdt=$tdt)\n")
            end
        elseif verb == 2 #if verbosity is high and there is no deflation
            print("🟩No deflation at j=$j ") #we inform the user
            #here, we want to estimate the minimum value of deftol that would have caught a "deflation"
            tdt = deftol #trial deftol starts as deftol
            for i in 1:20 #don't need a very fine approximation
                if size(rrqr(Rⱼ,tdt)[1],2) < m #if too high
                    tdt *= (deftol / 1.0)^(2.0^-i)
                else #if too high
                    tdt *= (1.0 / deftol)^(2.0^-i)
                end
            end
            print("(mdt=$tdt)\n")
        elseif verb == 1 #finally, if verbosity is medium and there is no deflation
            print("🟩") #just a green square with no newline character
        end
        
        Q = [Q Qʰ] #append the newly computed columns to Qⱼ
        Uʰ = U[1:Int(size(U,1)/2),size(U,2)-ℓ+1:size(U,2)] #we copy this because we want to modify it in the orthogonalisation without modifying Uⱼ
        H = [H zeros(size(H,1),ℓ);zeros(ℓ,size(H,2)) zeros(ℓ,ℓ)] #expand H (zeros for now)
        for i = 1:size(U,2) #second-level orthogonalisation
            H[i:i,size(H,2)-ℓ+1:size(H,2)] = U[1:Int(size(U,1)/2),i:i]'*Sⱼ + U[Int(size(U,1)/2)+1:size(U,1),i:i]'*Uʰ #fill in new block column of Hₖ
            Sⱼ -= U[1:Int(size(U,1)/2),i:i]*H[i:i,size(H,2)-ℓ+1:size(H,2)] #orthogonalise Sⱼ against Uⱼ,₁
            Uʰ -= U[Int(size(U,1)/2)+1:size(U,1),i:i]*H[i:i,size(H,2)-ℓ+1:size(H,2)] #and Uʰ against Uⱼ,₂
        end
        for i = 1:size(U,2) #second-level reorthogonalisation
            Hᵣₑₛ = U[1:Int(size(U,1)/2),i:i]'*Sⱼ + U[Int(size(U,1)/2)+1:size(U,1),i:i]'*Uʰ #not the full residual
            Sⱼ -= U[1:Int(size(U,1)/2),i:i]*Hᵣₑₛ #reorthogonalise Sⱼ against Uⱼ,₁
            Uʰ -= U[Int(size(U,1)/2)+1:size(U,1),i:i]*Hᵣₑₛ #and Uʰ against Uⱼ,₂
            H[i:i,size(H,2)-ℓ+1:size(H,2)] += Hᵣₑₛ #correct block column of Hₖ
        end
        Qᵗ,Rᵗ = qr([Sⱼ;Rʰ;Uʰ]) #standard QR factorisation, not RRQR
        Qᵗ = Matrix(Qᵗ) #we must force explicit formation of Qᵗ
        if rank(Rᵗ,rtol=deftol) < ℓ #this means (at least partial) breakdown of the concurrent block Arnoldi procedure
            if verb == 2
                rk = rank(Rᵗ,rtol=deftol)
                print("🟥Breakdown at j=$j (rank=$rk)\n")
            elseif verb == 1
                print("🟥j=$j\n")
            end
            @warn "breakdown at j=$j, returning only Qⱼ" #warn the user, to prevent silent failure from causing problems
            return Q,nothing,nothing #return Qⱼ because it might still be useful (probably not)
        end
        H[size(H,1)-ℓ+1:size(H,1),size(H,2)-ℓ+1:size(H,2)] = Rᵗ #fill in bottom-right block entry of Hₖ
        U = [[[U[1:Int(size(U,1)/2),:];zeros(size(Qᵗ,1)-size(U,1),size(U,2));U[Int(size(U,1)/2)+1:size(U,1),:]] Qᵗ];zeros(size(Qᵗ,1)-size(U,1),size(U,2)+size(Qᵗ,2))] #form new columns and rows of Uⱼ
    end
    if verb > 0
        print("🟦\n\n")
    end
    return Q,U,H
end

"""
    quadEigRBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix, req::Int=100, tol::Float64=1e-12, kℓₘₐₓ::Int, ℓ::Int; step::Int=10, σ::Union{Float64,ComplexF64}=0.0+0.0im, smallest::Bool=true, keep::Function=every, dtol::Float64=1e-10, rrv::Int=0, flvd::Bool=true, verb::Int=0, check_singular::Bool=false, give_up::Int=10)

Compute some eigenpairs of the QEP `(λ²M + λD + K)x=0` using the restarted block TOAR algorithm.

# Arguments
 -`M::AbstractMatrix`: mass matrix from QEP.\n
 -`D::AbstractMatrix`: damping matrix from QEP.\n
 -`K::AbstractMatrix`: stiffness matrix from QEP.\n
 -`req::Int`: required number of eigenpairs. Make sure this is at most `kℓₘₐₓ/2`. Note that the number of returned eigenpairs will often be slightly larger than `req`.\n
 -`tol::Float64`: maximum permissible backward error residual `ρ` for an eigenpair to be returned.\n
 -`kℓₘₐₓ::Int`: maximum subspace size before restart. Defaults to `300`, reduce this if memory consumption is an issue but set significantly larger than `req`.\n
 -`ℓ::Int`: block size/width. Defaults to `1`. It is not advised to set this higher than `5`.\n
 -`step::Int`: minimum number of blocks to add to the subspace between checks for convergence. Defaults to `10`, you may wish to set this lower for a higher `ℓ`.\n
 -`σ::Union{Float64,ComplexF64}`: shift point for shift-and-invert transformation. Defaults to `0.0`. Should be set within the domain of interest.\n
 -`smallest::Bool`: whether to invert the QEP. Inverting will find eigenvalues closest to `σ`, not inverting will find those furthest away. Defaults to `true`.\n
 -`keep::Function`: function that accepts a `ComplexF64` eigenvalue and returns whether it is within the domain of interest. Defaults to always true.\n
 -`dtol::Float64`: internal numerical tolerance for deflation/breakdown detection. Don't change this unless you know what you're doing.\n
 -`rrv::Int`: the number of inverse power iterations to use in Ritz vector refinement (default `0`). Not currently implemented.\n
 -`flvd::Bool`: whether to apply Fan, Lin & Van Dooren scaling to the QEP. Default `true`.\n
 -`verb::Int`: verbosity level. 0: no verbosity, 1: some verbosity, 2: full verbosity. Full verbosity has a large performance impact. (Not implemented yet.)\n
 -`check_singular::Bool`: whether to check if the QEP is close to being singular, default `true`. This test can be expensive and could give false positives for some QEPs.\n
 -`give_up::Int`: how many restarts to allow before terminating in failure.\n

# Returns
 -`λ::Vector`: array of Ritz values.\n
 -`X::Matrix`: array of Ritz vectors.\n
 -`ρ::Vector`: array of backward error residuals for returned eigenpairs `λ`,`X`.\n
"""
function quadEigRBTOAR(M::AbstractMatrix,D::AbstractMatrix,K::AbstractMatrix;req::Int=100,tol::Float64=1e-12,kℓₘₐₓ::Int=300,ℓ::Int=1,step::Int=10,σ::Union{Float64,ComplexF64}=0.0+0.0im,which::Symbol=:SM,keep::Function=every,dtol::Float64=1e-10,rrv::Int=0,arpack::Bool=true,flvd::Bool=true,verb::Int=0,check_singular::Bool=false,give_up::Int=10)
    n = size(M,1) #take n implicitly
    if false in (n .== [size(M,2);size(D,1);size(D,2);size(K,1);size(K,2)]) #M, D and K must all be n×n
        error("M, D and K must all be n×n")
    end
    
    if verb > 2
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 2"
        verb = 2
    end
    if verb < 0
        @warn "valid values of verb are 0, 1 and 2\nSetting verb to 0"
        verb = 0
    end
    if 2req > kℓₘₐₓ
        @warn "req should not be larger than kℓₘₐₓ/2; the algorithm may stagnate"
    end
    if ℓ*step > 50
        @warn "ℓ*step = $(ℓ*step), consider setting step lower to avoid building more subspace than necessary"
    end
    if step < 1
        @error "step must be positive"
    end
    if ℓ > 5
        @warn "it is not reccommended to set ℓ greater than 5 (ℓ = $ℓ)"
    end
    if ℓ < 1
        @error "ℓ must be positive"
    end
    if (dtol > 1e-6) || (dtol < 1e-15)
        @warn "bad value for dtol (dtol = $dtol)"
    end
    if which ∉ [:SM;:LM]
        @error "valid values of which are :SM and :LM"
    end
    if give_up < 3
        @warn "probably best to set give_up higher than $give_up"
    end
    if give_up < 1
        @error "give_up must be positive (give_up=$give_up)"
    end
    if give_up > 20
        @warn "probably best to set give_up lower than $give_up"
    end
    if n < 1000
        @warn "this eigensolver is not designed for small (n < 10³) QEPs. Consider using QuadEig.jl"
    end

    if rrv != 0
        @warn "refined Ritz vectors are not currently implemented"
    end

    inv = which == :SM #Fran wanted me to rename this argument
    
    if check_singular
        #normally a condition number of 1e10 isn't enough to be called numerically singular in Float64, but for TOAR we have to careful
        if (cond(M,1) > 1e10) && (cond(K,1) > 1e10) #only if both of these two are singular can the whole QEP be
            could_be_singular = true
            if verb == 2
                print("\n⚠️M and K both close to singular.\n\n")
            end
        end
    end

    Mₙₒᵣₘ = opnorm(M,1) #precompute matrix norms for efficiency
    Dₙₒᵣₘ = opnorm(D,1)
    Kₙₒᵣₘ = opnorm(K,1)
    badly_scaled_QEP = maximum(abs.(log.(10,[Mₙₒᵣₘ/Dₙₒᵣₘ,Mₙₒᵣₘ/Kₙₒᵣₘ,Dₙₒᵣₘ/Kₙₒᵣₘ]))) > 8 #the QEP is considered badly scaled if the ratio of any of the 1-norms of its coefficient matrices is more than 10⁸
    #Fan, Lin & Van Dooren scaling (2004)
    γ = flvd ? sqrt(Kₙₒᵣₘ/Mₙₒᵣₘ) : 1.0 #eigenvalue scaling
    δ = flvd ? 2/(Kₙₒᵣₘ+γ*Dₙₒᵣₘ) : 1.0 #uniform scaling
    if flvd
        still_badly_scaled = maximum(abs.(log.(10,[γ*Mₙₒᵣₘ/Dₙₒᵣₘ,γ^2*Mₙₒᵣₘ/Kₙₒᵣₘ,γ*Dₙₒᵣₘ/Kₙₒᵣₘ]))) > 8 #cancelled factors of γ and excluded δ since it makes no difference
        if still_badly_scaled
            @warn "QEP still badly scaled after Fan, Lin & Van Dooren scaling: γ²δ‖M‖₁=$(@sprintf("%.2g",γ^2*δ*Mₙₒᵣₘ)), γδ‖D‖₁=$(@sprintf("%.2g",γ*δ*Dₙₒᵣₘ)), δ‖K‖₁=$(@sprintf("%.2g",δ*Kₙₒᵣₘ))."
        end
    else
        if badly_scaled_QEP
            @warn "QEP is badly scaled: ‖M‖₁=$(@sprintf("%.2g",Mₙₒᵣₘ)), ‖D‖₁=$(@sprintf("%.2g",Dₙₒᵣₘ)), ‖K‖₁=$(@sprintf("%.2g",Kₙₒᵣₘ)). Consider setting flvd=true."
        end
    end

    if verb > 0 #some or all verbosity
        print("== SCALING INFORMATION ==\n")
        if flvd
            print("    Fan, Lin & Van Dooren scaling applied with γ=$(@sprintf("%.2g",γ)), δ=$(@sprintf("%.2g",δ)).\n    Pre-scaling matrix norms:\n        ‖M‖₁=$(@sprintf("%.2g",Mₙₒᵣₘ))\n        ‖D‖₁=$(@sprintf("%.2g",Dₙₒᵣₘ))\n        ‖K‖₁=$(@sprintf("%.2g",Kₙₒᵣₘ))\n    Scaled matrix norms:\n        γ²δ‖M‖₁=$(@sprintf("%.2g",γ^2*δ*Mₙₒᵣₘ))\n        γδ‖D‖₁=$(@sprintf("%.2g",γ*δ*Dₙₒᵣₘ))\n        δ‖K‖₁=$(@sprintf("%.2g",δ*Kₙₒᵣₘ))\n\n")
        else
            print("    No scaling applied.\n    Matrix norms:\n        ‖M‖₁=$(@sprintf("%.2g",Mₙₒᵣₘ))\n        ‖D‖₁=$(@sprintf("%.2g",Dₙₒᵣₘ))\n        ‖K‖₁=$(@sprintf("%.2g",Kₙₒᵣₘ))\n\n")
        end
    end
    
    Kₛ = δ*K + σ*δ*D + σ^2*δ*M #scaled and shifted matrices (I multiplied out the cancelling factors of γ)
    Dₛ = δ*γ*D + 2σ*δ*γ*M
    Mₛ = δ*γ^2*M

    if check_singular && could_be_singular #short-circuiting and so if check_singular is false, could_be_singular doesn't need to exist
        if cond(Kₛ,1) > 1e10
            @warn "QEP may be close to singular"
        end
    end
    
    Mᴸᵁ = lu(inv ? Kₛ : Mₛ) #factorise M for fast linear solves
    M⁻¹(x) = Mᴸᵁ\x #use factorised version of M
    M¹(x) = (inv ? Kₛ : Mₛ)*x
    D¹(x) = Dₛ*x #we call the functions D¹ and K¹ because the names D and K are taken
    K¹(x) = (inv ? Mₛ : Kₛ)*x

    
    transformed_keep(λ) = keep.((λ .- σ) .^ -1 ./ γ) #to keep the right eigenvalues regardless of spectral transformation
    transformed_keep_deflate_distant(λ) = transformed_keep(λ) .&& deflate_distant(λ) #to keep the right eigenvalues and deflate the third most distant
    
    if verb > 0
        print("== START OF BTOAR ALGORITHM ==\n\n")
    end

    ############################## EVERYTHING BEFORE THIS POINT IS BASICALLY THE SAME AS IN quadEigBTOAR() #################################

    #initialise intermediate matrices retained for efficiency
    MQ = zeros(ComplexF64,n,0)
    DQ = zeros(ComplexF64,n,0)
    KQ = zeros(ComplexF64,n,0)

    #initialise reduced order QEP as 0×0 matrices
    QᵀMQ = zeros(ComplexF64,0,0)
    QᵀDQ = zeros(ComplexF64,0,0)
    QᵀKQ = zeros(ComplexF64,0,0)

    Q,U,H,_,_,_ = BTOAR(M⁻¹,D¹,K¹,rand(ComplexF64,n,ℓ),maximum([step,Int(floor(req/ℓ))]),deftol=dtol,verb=verb) #initialise by building extra large step
    m = size(Q,2) #there might have been deflations
    good = 0 #number of acceptable eigenpairs computed
    restarts = 0 #restart counter for give_up

    #main loop
    while true
        #form the reduced-order QEP efficiently, taking advantage of prior computations and retained matrices
        #I'm pretty sure this is the most efficient way to do this
        mₗₐₛₜ = size(QᵀMQ,1) #m before the next step is taken
        MQ = [MQ M¹(Q[:,mₗₐₛₜ+1:m])] #next step-column of MQ
        QᵀMQ = [QᵀMQ Q[:,1:mₗₐₛₜ]'MQ[:,mₗₐₛₜ+1:m];Q[:,mₗₐₛₜ+1:m]'MQ[:,1:mₗₐₛₜ] Q[:,mₗₐₛₜ+1:m]'MQ[:,mₗₐₛₜ+1:m]] #expand
        DQ = [DQ D¹(Q[:,mₗₐₛₜ+1:m])]
        QᵀDQ = [QᵀDQ Q[:,1:mₗₐₛₜ]'DQ[:,mₗₐₛₜ+1:m];Q[:,mₗₐₛₜ+1:m]'DQ[:,1:mₗₐₛₜ] Q[:,mₗₐₛₜ+1:m]'DQ[:,mₗₐₛₜ+1:m]]
        KQ = [KQ K¹(Q[:,mₗₐₛₜ+1:m])]
        QᵀKQ = [QᵀKQ Q[:,1:mₗₐₛₜ]'KQ[:,mₗₐₛₜ+1:m];Q[:,mₗₐₛₜ+1:m]'KQ[:,1:mₗₐₛₜ] Q[:,mₗₐₛₜ+1:m]'KQ[:,mₗₐₛₜ+1:m]]

        #solve reduced-order problem
        LP = linearize(QᵀKQ,QᵀDQ,QᵀMQ) #linearize the QEP
        Z = zeros(ComplexF64,m,2m) #preallocate space for reduced-order eigenvectors
        λ,V = eigen(LP.A,LP.B)
        for i in 1:2m #extract quadratic eigenvectors
            Z[:,i] = (abs(λ[i]) > 1) ? V[1:m,i] : V[m+1:2m,i] #this should be more stable I think
            Z[:,i] /= norm(Z[:,i]) #normalise
        end

        #uninvert, unscale, project back to full space
        λ .^= (inv ? -1 : 1) #if QEP was inverted, uninvert eigenvalues
        λ .+= σ/γ #unshift eigenvalues
        λ .*= γ #unscale eigenvalues
        X = Q*Z #transform eigenvectors back to full-sized space (they stay normalised)

        #compute residuals
        ρ = zeros(2m)
        for i in 1:2m #calculate backward error residuals
            ρ[i] = norm(λ[i]^2*M*X[:,i]+λ[i]*D*X[:,i]+K*X[:,i]) / (abs2(λ[i])*Mₙₒᵣₘ+abs(λ[i])*Dₙₒᵣₘ+Kₙₒᵣₘ)
        end
        good = sum((ρ .< tol) .&& keep.(λ)) #number of acceptable residuals in domain of interest
        if verb == 1
            print("Subspace size: $m / $kℓₘₐₓ\nGood eigenpairs: $good\n\n")
        elseif verb == 2
            print("Subspace size: $m / $kℓₘₐₓ\nTotal good eigenpairs: $(sum((ρ .< tol)))\nGood eigenpairs in DoI: $good\n\n")
        end

        if good ≥ req #if we have found enough acceptable eigenpairs
            if verb == 2
                print("$good good eigenpairs found, returning.")
            end
            good_ones = (ρ .< tol) .&& keep.(λ) #basically nothing to recompute this: O(kℓ)
            λ = [λ[i] for i in 1:size(λ,1) if good_ones[i]]
            X = hcat([X[:,i] for i in 1:size(X,2) if good_ones[i]]...)
            ρ = [ρ[i] for i in 1:size(ρ,1) if good_ones[i]]
            return λ,X,ρ #we're done here
        elseif m+step*ℓ > kℓₘₐₓ #if another step could expand the subspace too far
            if restarts == give_up
                @warn "restart limit exceeded, not enough eigenpairs found"
                if verb > 0
                    print("⚠️Restart limit ($give_up) exceeded, giving up. (good eigenpairs: $good)\n\n")
                end
                good_ones = (ρ .< tol) .&& keep.(λ) #basically nothing to recompute this: O(kℓ)
                λ = [λ[i] for i in 1:size(λ,1) if good_ones[i]]
                X = hcat([X[:,i] for i in 1:size(X,2) if good_ones[i]]...)
                ρ = [ρ[i] for i in 1:size(ρ,1) if good_ones[i]]
                return λ,X,ρ #return at least what we have
            end
            restarts += 1
            if verb > 0
                print("== RESTART =="*"\n"^(3-verb)) #fancy way of getting the number of newlines right
            end
            if sum(keep.(λ)) + step*ℓ > kℓₘₐₓ #if deflating by keep() would not be enough
                Q,U,H = restartBTOAR(Q,U,H,transformed_keep_deflate_distant,verb)
                if verb == 2
                    print("Deflating according to keep() and ⅓ most distant.\n\n")
                end
            else #if deflating by keep() is enough
                Q,U,H = restartBTOAR(Q,U,H,transformed_keep,verb)
                if verb == 2
                    print("Deflating only according to keep().\n\n")
                end
            end
            MQ = M¹(Q)
            DQ = D¹(Q)
            KQ = K¹(Q)
            QᵀMQ = Q'*M¹(Q) #no way to avoid this cost really
            QᵀDQ = Q'*D¹(Q)
            QᵀKQ = Q'*K¹(Q)
        else #if we are free to grow the subspace by step
            if verb > 0
                print("== CONTINUING BTOAR ALGORITHM ==\n")
            end
            Q,U,H = continueBTOAR(M⁻¹,D¹,K¹,Q,U,H,maximum([step,minimum([Int(floor((req-good)/ℓ)),Int(floor((kℓₘₐₓ-m)/ℓ))])]),ℓ;deftol=dtol,verb=verb) #grow the subspace by as much as we can without overflowing or overdoing it
        end
        m = size(Q,2)
    end
end
