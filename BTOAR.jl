using LinearAlgebra,LowRankApprox,QuadEig,Printf

function fmean(v::Vector,f::Function,f⁻¹::Function)
    return f⁻¹(sum(f.(v))/size(v)[1])
end

function b(r::Union{Int,UnitRange{Int}},l::Int) #function to help out with long block indices with constant block size
    if typeof(r) <: Int
        return (r-1)l+1:r*l
    else
        return (minimum(r)-1)l+1:maximum(r)*l
    end
end

#this function doesn't actually compute a rrqr factorisation because R will not necessarily be upper trapezoidal, but we don't really care
function rrqr(A::Matrix,tol::Float64)
    Q,R,p = pqr(A,rtol=tol) #use the so-called `partial QR factorisation'
    return Q,R*ColumnPermutation(p)' #return Q and (unfortunately non-upper trapezoidal) R
end

hermitify(A) = (A+A')/2 #for the Cholesky factorisation, which thinks a matrix that's not EXACTLY Hermitian isn't acceptable

every(λ) = true #this function is just here as the default for the keep function in quadEigRBTOAR

"""
    BTOAR(M⁻¹::Function, D::Function, K::Function, R::Matrix{Complex{Float64}}, k::Int, deftol::Float64=NaN, verb::Int=0)
    
Compute an orthogonal basis for the `k`th second-order Krylov subspace Gₖ(M⁻¹D,M⁻¹K;R).

# Arguments
 -`M⁻¹::Function`: function that provides left-multiplication of `n`×`l` matrices by the inverse of M.\n
 -`D::Function`: function that provides left-multiplication of `n`×`l` matrices by the matrix D.\n
 -`K::Function`: function that provides left-multiplication of `n`×`l` matrices by the matrix K.\n
 -`R::Matrix{Complex{Float64}}`: starting `n`×`l` block vector R.\n
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
    n,l = size(R) #n and l are not function arguments, but are taken implicitly from R

    if deftol ≡ NaN #if deftol isn't set, we give it a reasonable default
        deftol = l*eps(Float64) #machine epsilon times the smallest dimension of the matrix is a standard numerical rank tolerance
    end
    if deftol < eps(Float64) #warning about setting deflation tolerance too low
        @warn "deftol should not be set lower than ϵ≈2.22×10⁻¹⁶ (deftol=$deftol)\nSetting deftol to ϵ"
        deftol = eps(Float64)
    end
    if deftol ≥ 1.0 #not sure what the highest reasonable deftol would be
        @warn "deftol way too large (deftol=$deftol)\nSetting deftol to ϵl=$(l*eps(Float64))"
        deftol = l*eps(Float64)
    end
    if k*l > n #warning about setting k too large
        @warn "kl greater than n, expect deflation (kl = $(k*l), n = $n)"
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
    m[1] = rank(R[1:l,:],rtol=deftol) #to save time, we start by taking the rank of a small full-width submatrix of R
    if m[1] < l #if the small full-width submatrix was rank deficient (unlikely most of the time) 
        m[1] = rank(R,rtol=deftol) #test the rank of the full block vector R
        if m[1] < l #if R is rank-deficient
            R = Matrix(rrqr(R)[1]) #we reduce R to an orthonormal matrix with the same span, using the RRQR factorisation
            l = m[1] #the new l after R has been reduced
            if verb > 0
                print("Starting vector R rank deficient, reducing l to $l\n")
            end
        end
    end
    Q,R = qr(R) #standard QR factorisation because R now must be full-rank
    Qⱼ = Matrix{ComplexF64}(Q) #qr() outputs Q in a special form without explicitly forming the matrix, so we have to tell it to
    Uⱼ = [I(l);zeros(Complex{Float64},l,l)] #because V₁ = [Q₁;0]
    Hₖ = zeros(Complex{Float64},k*l,(k-1)*l) #we can preallocate H because its final size is known now
    
    for j in 1:k-1 #main for loop
        Rⱼ = M⁻¹(D(Qⱼ*-Uⱼ[1:sum(m),(j-1)l+1:j*l]) + K(Qⱼ*-Uⱼ[sum(m)+1:2sum(m),(j-1)l+1:j*l])) #take next Rⱼ block vector
        Sⱼ = zeros(Complex{Float64},sum(m[1:j]),l) #preallocate Sⱼ
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
        if m[j+1] < l #if we have deflation
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
        Uʰ = Uⱼ[1:sum(m[1:j]),(j-1)l+1:j*l] #we copy this because we want to modify it in the orthogonalisation without modifying Uⱼ
        for i = 1:j #second-level orthogonalisation
            Hₖ[(i-1)l+1:i*l,(j-1)l+1:j*l] = Uⱼ[1:sum(m[1:j]),(i-1)l+1:i*l]'*Sⱼ + Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),(i-1)l+1:i*l]'*Uʰ #fill in new block column of Hₖ
            Sⱼ -= Uⱼ[1:sum(m[1:j]),(i-1)l+1:i*l]*Hₖ[(i-1)l+1:i*l,(j-1)l+1:j*l] #orthogonalise Sⱼ against Uⱼ,₁
            Uʰ -= Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),(i-1)l+1:i*l]*Hₖ[(i-1)l+1:i*l,(j-1)l+1:j*l] #and Uʰ against Uⱼ,₂
        end
        for i = 1:j #second-level reorthogonalisation
            Hᵣₑₛ = Uⱼ[1:sum(m[1:j]),(i-1)l+1:i*l]'*Sⱼ + Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),(i-1)l+1:i*l]'*Uʰ #not the full residual
            Sⱼ -= Uⱼ[1:sum(m[1:j]),(i-1)l+1:i*l]*Hᵣₑₛ #reorthogonalise Sⱼ against Uⱼ,₁
            Uʰ -= Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),(i-1)l+1:i*l]*Hᵣₑₛ #and Uʰ against Uⱼ,₂
            Hₖ[(i-1)l+1:i*l,(j-1)l+1:j*l] += Hᵣₑₛ #correct block column of Hₖ
        end
        Qᵗ,Rᵗ = qr([Sⱼ;Rʰ;Uʰ]) #standard QR factorisation, not RRQR
        Qᵗ = Matrix(Qᵗ) #we must force explicit formation of Qᵗ
        if rank(Rᵗ,rtol=deftol) < l #this means (at least partial) breakdown of the concurrent block Arnoldi procedure
            if verb == 2
                rk = rank(Rᵗ,rtol=deftol)
                print("🟥Breakdown at j=$j (rank=$rk)\n")
            elseif verb == 1
                print("🟥j=$j\n")
            end
            @warn "breakdown at j=$j, returning only Qⱼ" #warn the user, to prevent silent failure from causing problems
            return Qⱼ,nothing,nothing,nothing,nothing,nothing #return Qⱼ because it might still be useful (probably not)
        end
        Hₖ[j*l+1:(j+1)l,(j-1)l+1:j*l] = Rᵗ #fill in bottom-right block entry of Hₖ
        Uⱼ = [[[Uⱼ[1:sum(m[1:j]),:];zeros(m[j+1],j*l);Uⱼ[sum(m[1:j])+1:2sum(m[1:j]),:]] Qᵗ];zeros(m[j+1],(j+1)l)] #form new columns and rows of Uⱼ
        if verb == 2
            print("TOAR relation residual: ")
            display(opnorm([-(M⁻¹(D(Qⱼ[:,1:sum(m[1:j])]*Uⱼ[1:sum(m[1:j]),1:j*l]) + K(Qⱼ[:,1:sum(m[1:j])]*Uⱼ[sum(m)+1:sum(m)+sum(m[1:j]),1:j*l])));Qⱼ[:,1:sum(m[1:j])]*Uⱼ[1:sum(m[1:j]),1:j*l]] - [Qⱼ*Uⱼ[1:sum(m),:];Qⱼ*Uⱼ[sum(m)+1:2sum(m),:]]*Hₖ[1:(j+1)l,1:j*l],1) / opnorm([Qⱼ*Uⱼ[1:sum(m),:];Qⱼ*Uⱼ[sum(m)+1:2sum(m),:]]*Hₖ[1:(j+1)l,1:j*l],1)) #this was soul-crushing to debug
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
    ρ₃ = opnorm([-(M⁻¹(D(Qⱼ[:,1:sum(m[1:j])]*Uⱼ[1:sum(m[1:j]),1:j*l]) + K(Qⱼ[:,1:sum(m[1:j])]*Uⱼ[sum(m)+1:sum(m)+sum(m[1:j]),1:j*l])));Qⱼ[:,1:sum(m[1:j])]*Uⱼ[1:sum(m[1:j]),1:j*l]] - [Qⱼ*Uⱼ[1:sum(m),:];Qⱼ*Uⱼ[sum(m)+1:2sum(m),:]]*Hₖ[1:(j+1)l,1:j*l],1) / opnorm([Qⱼ*Uⱼ[1:sum(m),:];Qⱼ*Uⱼ[sum(m)+1:2sum(m),:]]*Hₖ[1:(j+1)l,1:j*l],1) #this was soul-crushing to debug
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
    quadEigBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix, k::Int, l::Int; σ::Union{Float64,Complex{Float64}}=0.0+0.0im, inv::Bool=true, dtol::Float64=1e-10, rrv::Int=0, flvd::Bool=true, verb::Int=0, check_singular::Bool=false)
    
Compute some eigenpairs of the QEP `(λ²M + λD + K)x=0` using the block TOAR algorithm.

# Arguments:
 -`M::AbstractMatrix`: mass matrix of the QEP.\n
 -`D::AbstractMatrix`: damping matrix of the QEP.\n
 -`K::AbstractMatrix`: stiffness matrix of the QEP.\n
 -`k::Int`: number of block TOAR iterations.\n
 -`l::Int`: block size/width.\n
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
function quadEigBTOAR(M::AbstractMatrix,D::AbstractMatrix,K::AbstractMatrix,k::Int,l::Int;σ::Union{Float64,Complex{Float64}}=0.0+0.0im,inv::Bool=true,dtol::Float64=1e-10,rrv::Int=0,flvd::Bool=true,verb::Int=0,check_singular::Bool=false)
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
    Qₖ,_,_,_,_,_ = BTOAR(M⁻¹,D¹,K¹,rand(ComplexF64,n,l),k,deftol=dtol,verb=verb) #run the BTOAR algorithm
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
        MQₖᵀDQₖ = MQₖ'DQₖ #these multiplications are (kl×n)(n×kl)
        MQₖᵀKQₖ = MQₖ'KQₖ
        DQₖᵀMQₖ = DQₖ'MQₖ
        DQₖᵀDQₖ = DQₖ'DQₖ
        DQₖᵀKQₖ = DQₖ'KQₖ #yeah, there are a lot of them but it's most efficient way if m >> 9
        KQₖᵀMQₖ = KQₖ'MQₖ
        KQₖᵀDQₖ = KQₖ'DQₖ
        KQₖᵀKQₖ = KQₖ'KQₖ
        for i in 1:2m
            #additions here are only of kl×kl matrices (cheap)
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
    l = size(H□,1)-size(H□,2)
    
    Hₚ₊₁ = [schurfact.T[1:p,1:p];H□[size(H□,2)+1:size(H□,1),:]*schurfact.Z[:,1:p]] #form new H
    
    U□⁽¹⁾ = U□[1:size(Q□,2),:] #deconstruct U□ for readability
    U□⁽²⁾ = U□[size(Q□,2)+1:2size(Q□,2),:]
    U,Σ,V = svd([[U□⁽¹⁾[1:size(U□⁽¹⁾,1)-l,1:size(U□,2)-l]*schurfact.Z[:,1:p];zeros(l,p)] U□⁽¹⁾[:,size(U□,2)-l+1:size(U□,2)] [U□⁽²⁾[1:size(U□⁽²⁾,1)-l,1:size(U□,2)-l]*schurfact.Z[:,1:p];zeros(l,p)] U□⁽²⁾[:,size(U□,2)-l+1:size(U□,2)]]) #that took a while to type
    U = U[:,1:minimum([p+2l,size(U,2)])] #this minimum() is in case the number of rows of [W X Y Z] is too small (it is (k+1)l, which could be smaller than p+2l if l is big enough)
    Σ = Σ[1:minimum([p+2l,size(Σ,1)])]
    V = V[:,1:p+2l]
    
    Uₚ₊₁ = kron(I(2),Diagonal(Σ))*[V'[:,1:p+l];V'[:,p+l+1:2(p+l)]]
    
    Qᵣ = Q□*U #dominant flop cost
    
    if verb > 0 #brief verbose output of number of discarded/retained dimensions
        print("Dimensions retained: $(size(Qᵣ,2))\n")
        print("Dimensions purged: $(size(Q□,2)-size(Qᵣ,2))\n\n")
    end
    
    return Qᵣ,Uₚ₊₁,Hₚ₊₁
end

"""
    continueBTOAR(M⁻¹::Function, D::Function, K::Function, Qᵣ::Matrix{Complex{Float64}}, Uₚ₊₁::Matrix{Complex{Float64}}, Hₚ₊₁::Matrix{Complex{Float64}}, k::Int, l::Int; deftol::Float64=1e-10)

Continue the BTOAR algorithm after a restart.

# Arguments
 -`M⁻¹::Function`: function that provides left-multiplication of `n`×`l` matrices by the inverse of M.\n
 -`D::Function`: function that provides left-multiplication of `n`×`l` matrices by the matrix D.\n
 -`K::Function`: function that provides left-multiplication of `n`×`l` matrices by the matrix K.\n
 -`Qᵣ::Matrix{Complex{Float64}}`: locked second-order Krylov basis.\n
 -`Uₚ₊₁::Matrix{Complex{Float64}}`: locked auxiliary matrix.\n
 -`Hₚ₊₁::Matrix{Complex{Float64}}`: locked auxiliary matrix.\n
 -`k::Int`: number of new blocks to add to the subspace.\n
 -`l::Int`: block size/width.\n
 -`deftol::Float64`: internal numerical tolerance for deflation detection, defaults to `1e-10`.\n
 -`verb::Int`: verbosity option. 0: no verbosity, 1: some verbosity, 2: full verbosity. Full verbosity has a large performance impact.\n

# Returns
 -`Q::Matrix`: second-order Krylov subspace basis.\n
 -`U::Matrix`: auxiliary matrix.\n
 -`H::Matrix`: auxiliary matrix.\n
"""
function continueBTOAR(M⁻¹::Function,D::Function,K::Function,Qᵣ::Matrix,Uₚ₊₁::Matrix,Hₚ₊₁::Matrix,k::Int,l::Int;deftol::Float64=1e-10,verb::Int=0)
    n = size(Qᵣ,1) #take n implicitly
    
    if deftol < eps(Float64) #warning about setting deflation tolerance too low
        @warn "deftol should not be set lower than ϵ≈2.22×10⁻¹⁶ (deftol=$deftol)\nSetting deftol to ϵ"
        deftol = eps(Float64)
    end
    if deftol ≥ 1.0 #not sure what the highest reasonable deftol would be
        @warn "deftol way too large (deftol=$deftol)\nSetting deftol to ϵl=$(l*eps(Float64))"
        deftol = l*eps(Float64)
    end
    if k*l+size(Qᵣ,2) > n #warning about setting k too large
        @warn "kl+r greater than n, expect deflation (kl+r = $(k*l+size(Qᵣ,2)), n = $n)"
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
        Rⱼ = M⁻¹(D(Q*-U[1:Int(size(U,1)/2),size(U,2)-l+1:size(U,2)]) + K(Q*-U[Int(size(U,1)/2)+1:size(U,1),size(U,2)-l+1:size(U,2)])) #take next Rⱼ block vector
        Sⱼ = zeros(Complex{Float64},size(Q,2),l) #preallocate Sⱼ
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
        if m < l #if we have deflation
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
        Uʰ = U[1:Int(size(U,1)/2),size(U,2)-l+1:size(U,2)] #we copy this because we want to modify it in the orthogonalisation without modifying Uⱼ
        H = [H zeros(size(H,1),l);zeros(l,size(H,2)) zeros(l,l)] #expand H (zeros for now)
        for i = 1:size(U,2) #second-level orthogonalisation
            H[i:i,size(H,2)-l+1:size(H,2)] = U[1:Int(size(U,1)/2),i:i]'*Sⱼ + U[Int(size(U,1)/2)+1:size(U,1),i:i]'*Uʰ #fill in new block column of Hₖ
            Sⱼ -= U[1:Int(size(U,1)/2),i:i]*H[i:i,size(H,2)-l+1:size(H,2)] #orthogonalise Sⱼ against Uⱼ,₁
            Uʰ -= U[Int(size(U,1)/2)+1:size(U,1),i:i]*H[i:i,size(H,2)-l+1:size(H,2)] #and Uʰ against Uⱼ,₂
        end
        for i = 1:size(U,2) #second-level reorthogonalisation
            Hᵣₑₛ = U[1:Int(size(U,1)/2),i:i]'*Sⱼ + U[Int(size(U,1)/2)+1:size(U,1),i:i]'*Uʰ #not the full residual
            Sⱼ -= U[1:Int(size(U,1)/2),i:i]*Hᵣₑₛ #reorthogonalise Sⱼ against Uⱼ,₁
            Uʰ -= U[Int(size(U,1)/2)+1:size(U,1),i:i]*Hᵣₑₛ #and Uʰ against Uⱼ,₂
            H[i:i,size(H,2)-l+1:size(H,2)] += Hᵣₑₛ #correct block column of Hₖ
        end
        Qᵗ,Rᵗ = qr([Sⱼ;Rʰ;Uʰ]) #standard QR factorisation, not RRQR
        Qᵗ = Matrix(Qᵗ) #we must force explicit formation of Qᵗ
        if rank(Rᵗ,rtol=deftol) < l #this means (at least partial) breakdown of the concurrent block Arnoldi procedure
            if verb == 2
                rk = rank(Rᵗ,rtol=deftol)
                print("🟥Breakdown at j=$j (rank=$rk)\n")
            elseif verb == 1
                print("🟥j=$j\n")
            end
            @warn "breakdown at j=$j, returning only Qⱼ" #warn the user, to prevent silent failure from causing problems
            return Q,nothing,nothing #return Qⱼ because it might still be useful (probably not)
        end
        H[size(H,1)-l+1:size(H,1),size(H,2)-l+1:size(H,2)] = Rᵗ #fill in bottom-right block entry of Hₖ
        U = [[[U[1:Int(size(U,1)/2),:];zeros(size(Qᵗ,1)-size(U,1),size(U,2));U[Int(size(U,1)/2)+1:size(U,1),:]] Qᵗ];zeros(size(Qᵗ,1)-size(U,1),size(U,2)+size(Qᵗ,2))] #form new columns and rows of Uⱼ
    end
    if verb > 0
        print("🟦\n\n")
    end
    return Q,U,H
end

"""
    quadEigRBTOAR(M::AbstractMatrix, D::AbstractMatrix, K::AbstractMatrix, req::Int=100, tol::Float64=1e-10, kl_max::Int, l::Int; step::Int=10, σ::Union{Float64,ComplexF64}=0.0+0.0im, smallest::Bool=true, keep::Function=every, dtol::Float64=1e-10, rrv::Int=0, flvd::Bool=true, verb::Int=0, check_singular::Bool=false, give_up::Int=10)

Compute some eigenpairs of the QEP `(λ²M + λD + K)x=0` using the restarted block TOAR algorithm.

# Arguments
 -`M::AbstractMatrix`: mass matrix from QEP.\n
 -`D::AbstractMatrix`: damping matrix from QEP.\n
 -`K::AbstractMatrix`: stiffness matrix from QEP.\n
 -`req::Int`: required number of eigenpairs. Make sure this is at most `kl_max/2` (sometimes much lower depending on the QEP). Note that the number of returned eigenpairs will often be slightly larger than `req`.\n
 -`tol::Float64`: maximum permissible backward error residual `ρ` for an eigenpair to be returned.\n
 -`kl_max::Int`: maximum subspace size before restart. Defaults to `300`, reduce this if memory consumption is an issue but set it significantly larger than `req`.\n
 -`l::Int`: block size/width. Defaults to `1`. It is not advised to set this higher than `5`.\n
 -`step::Int`: minimum number of blocks to add to the subspace between checks for convergence. Defaults to `10`, you may wish to set this lower for a higher `l`.\n
 -`σ::Union{Float64,ComplexF64}`: shift point for shift-and-invert transformation. Defaults to `0.0`. Should be set within the domain of interest.\n
 -`which::Symbol`: which eigenvalues to target. `:SM` will target eigenvalues closest to `σ`, `:LM` targest those furthest from `σ`. Defaults to `:SM`.\n
 -`keep::Function`: function that accepts a `ComplexF64` eigenvalue and returns whether it is within the domain of interest. Defaults to always true.\n
 -`dtol::Float64`: internal numerical tolerance for deflation/breakdown detection. Don't change this unless you know what you're doing.\n
 -`rrv::Int`: the number of inverse power iterations to use in Ritz vector refinement (default `0`). Not currently implemented.\n
 -`flvd::Bool`: whether to apply Fan, Lin & Van Dooren scaling to the QEP. Default `true`.\n
 -`verb::Int`: verbosity level. 0: no verbosity, 1: some verbosity, 2: full verbosity. Full verbosity has a large performance impact (not fully implemented yet).\n
 -`check_singular::Bool`: whether to check if the QEP is close to being singular, default `true`. This test can be expensive and could give false positives for some QEPs.\n
 -`give_up::Int`: how many restarts to allow before terminating in failure.\n
 -`glob::Bool`: whether to store the \"best yet\" computed eigenvalues, eigenvectors and residuals in global variables `glob_λ`, `glob_X` and `glob_ρ` during execution (default `false`). This allows manual interruption of the function without losing the results.\n
 -`extra_space::Int`: how much to allow the subspace to be slightly larger than the number of required eigenvalues. Not implemented yet, will have no effect.\n

# Returns
 -`λ::Vector`: array of Ritz values.\n
 -`X::Matrix`: array of Ritz vectors.\n
 -`ρ::Vector`: array of backward error residuals for returned eigenpairs `λ`,`X`.\n
"""
function quadEigRBTOAR(M::AbstractMatrix,D::AbstractMatrix,K::AbstractMatrix;req::Int=100,tol::Float64=1e-10,kl_max::Int=300,l::Int=1,step::Int=10,σ::Union{Float64,ComplexF64}=0.0+0.0im,which::Symbol=:SM,keep::Function=every,dtol::Float64=1e-10,rrv::Int=0,flvd::Bool=true,verb::Int=0,check_singular::Bool=false,give_up::Int=10,glob::Bool=false,extra_space::Int=0)
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
    if req > kl_max
        @warn "req should not be larger than kl_max; the algorithm may stagnate"
    end
    if l*step > 50
        @warn "l*step = $(l*step), consider setting step lower to avoid building more subspace than necessary"
    end
    if step < 1
        @error "step must be positive (got $step)"
    end
    if l > 5
        @warn "it is not reccommended to set l greater than 5 (l = $l)"
    end
    if l < 1
        @error "l must be positive (got $l)"
    end
    if (dtol > 1e-6) || (dtol < 1e-15)
        @warn "bad value for dtol (dtol = $dtol)"
    end
    if which ∉ [:SM;:LM]
        @error "valid values of which are :SM and :LM (got $which)"
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
    if rrv < 0
        @error "rrv must be positive (got $rrv)"
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
    
    transformed_keep(λ) = keep.(((λ ./ γ) .- σ) .^ -1) #to keep the right eigenvalues regardless of spectral transformation
    
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
    
    #if we are using RRVs, initialise RRV quantities
    if rrv > 0
        MQᵀMQ = zeros(ComplexF64,0,0)
        MQᵀDQ = zeros(ComplexF64,0,0)
        MQᵀKQ = zeros(ComplexF64,0,0)
        DQᵀMQ = zeros(ComplexF64,0,0)
        DQᵀDQ = zeros(ComplexF64,0,0)
        DQᵀKQ = zeros(ComplexF64,0,0)
        KQᵀMQ = zeros(ComplexF64,0,0)
        KQᵀDQ = zeros(ComplexF64,0,0)
        KQᵀKQ = zeros(ComplexF64,0,0)
    end

    Q,U,H,_,_,_ = BTOAR(M⁻¹,D¹,K¹,rand(ComplexF64,n,l),maximum([step,Int(floor(req/l))]),deftol=dtol,verb=verb) #initialise by building extra large step
    m = size(Q,2) #there might have been deflations
    good = 0 #number of acceptable eigenpairs computed
    good_ones = Bool[] #preallocate this as empty to avoid crash
    restarts = 0 #restart counter for give_up

    global glob_λ,glob_X,glob_ρ #allows the user to interrupt without losing information
    best_λ = ComplexF64[]; best_X = zeros(ComplexF64,n,0); best_ρ = ComplexF64[] #assign them empty just in case something was in them before

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
        if rrv > 0
            thing = MQ[:,1:mₗₐₛₜ]'MQ[:,mₗₐₛₜ+1:m] #couldn't think of a good variable name lol
            MQᵀMQ = [MQᵀMQ thing;thing' MQ[:,mₗₐₛₜ+1:m]'MQ[:,mₗₐₛₜ+1:m]]
            MQᵀDQ = [MQᵀDQ MQ[:,1:mₗₐₛₜ]'DQ[:,mₗₐₛₜ+1:m];MQ[:,mₗₐₛₜ+1:m]'DQ[:,1:mₗₐₛₜ] MQ[:,mₗₐₛₜ+1:m]'DQ[:,mₗₐₛₜ+1:m]]
            MQᵀKQ = [MQᵀKQ MQ[:,1:mₗₐₛₜ]'KQ[:,mₗₐₛₜ+1:m];MQ[:,mₗₐₛₜ+1:m]'KQ[:,1:mₗₐₛₜ] MQ[:,mₗₐₛₜ+1:m]'KQ[:,mₗₐₛₜ+1:m]]
            DQᵀMQ = MQᵀDQ' #take advantage of past computations (this is copied by ref, so no expensive memory copying)
            thing = DQ[:,1:mₗₐₛₜ]'DQ[:,mₗₐₛₜ+1:m]
            DQᵀDQ = [DQᵀDQ thing;thing' DQ[:,mₗₐₛₜ+1:m]'DQ[:,mₗₐₛₜ+1:m]]
            DQᵀKQ = [DQᵀKQ DQ[:,1:mₗₐₛₜ]'KQ[:,mₗₐₛₜ+1:m];DQ[:,mₗₐₛₜ+1:m]'KQ[:,1:mₗₐₛₜ] DQ[:,mₗₐₛₜ+1:m]'KQ[:,mₗₐₛₜ+1:m]]
            KQᵀMQ = MQᵀKQ'
            KQᵀDQ = DQᵀKQ'
            thing = KQ[:,1:mₗₐₛₜ]'KQ[:,mₗₐₛₜ+1:m]
            KQᵀKQ = [KQᵀKQ thing;thing' KQ[:,mₗₐₛₜ+1:m]'KQ[:,mₗₐₛₜ+1:m]]

            #display([opnorm(MQᵀMQ-MQ'MQ)/opnorm(MQᵀMQ);opnorm(MQᵀDQ-MQ'DQ)/opnorm(MQᵀDQ);opnorm(MQᵀKQ-MQ'KQ)/opnorm(MQᵀKQ);opnorm(DQᵀMQ-DQ'MQ)/opnorm(DQᵀMQ);opnorm(DQᵀDQ-DQ'DQ)/opnorm(DQᵀDQ);opnorm(DQᵀKQ-DQ'KQ)/opnorm(DQᵀKQ);opnorm(KQᵀMQ-KQ'MQ)/opnorm(KQᵀMQ);opnorm(KQᵀDQ-KQ'DQ)/opnorm(KQᵀDQ);opnorm(KQᵀKQ-KQ'KQ)/opnorm(KQᵀKQ)])
        end

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
            ρ[i] = keep(λ[i]) ? (norm(λ[i]^2*M*X[:,i]+λ[i]*D*X[:,i]+K*X[:,i]) / (abs2(λ[i])*Mₙₒᵣₘ+abs(λ[i])*Dₙₒᵣₘ+Kₙₒᵣₘ)) : 1e100 #if not in DoI, set residual to arbitrary high value
        end
        good = sum((ρ .< tol) .&& keep.(λ)) #number of acceptable residuals in domain of interest THIS LINE COULD ACTUALLY BE SIMPLIFIED NOW BECAUSE OF L.761 OPTIMISATION
        if verb == 1
            print("Subspace size: $m / $kl_max\nGood eigenpairs: $good / $req\n\n")
        elseif verb == 2
            print("Subspace size: $m / $kl_max\nTotal good eigenpairs: $(sum((ρ .< tol))) / $req\nGood eigenpairs in DoI: $good / $req\n\n")
        end
        
        #Refined Ritz vector part has to be after residual calculation
        if rrv > 0 #can't do truthiness of integers in Julia
            if good < req #no need to do anything if we already have enough good eigenpairs
                for i in 1:size(Z,2)
                    if (ρ[i] > tol) && keep(λ[i]) #we only refine ones with bad residuals in the DoI
                        θ = ((λ[i]/γ)-σ/γ) ^ (inv ? -1 : 1) #retransformed eigenvalue
                        PQᵀPQ = lu(hermitify(θ'^2*(θ^2*MQᵀMQ + θ*MQᵀDQ + MQᵀKQ) + θ'*(θ^2*DQᵀMQ + θ*DQᵀDQ + DQᵀKQ) + θ^2*KQᵀMQ + θ*KQᵀDQ + KQᵀKQ))
                        for j in 1:rrv
                            Z[:,i] = PQᵀPQ\Z[:,i] #what could be simpler than an inverse power iteration
                            Z[:,i] /= norm(Z[:,i]) #better not forget this lol
                        end
                        X[:,i] = Q*Z[:,i] #update full-size vector
                        if verb == 2
                            print("Residual before RVR: $(ρ[i])\n") #put this in verb == 2 block
                        end
                        ρ[i] = norm(λ[i]^2*M*X[:,i]+λ[i]*D*X[:,i]+K*X[:,i]) / (abs2(λ[i])*Mₙₒᵣₘ+abs(λ[i])*Dₙₒᵣₘ+Kₙₒᵣₘ) #update residual
                        if verb == 2
                            print("Residual after: $(ρ[i])\n\n")
                        end
                    end
                end
                good = sum((ρ .< tol) .&& keep.(λ)) #number of acceptable residuals in domain of interest
                if verb == 1
                    print("== REFINED RITZ VECTORS COMPUTED ==\nGood eigenpairs: $good / $req\n\n")
                elseif verb == 2
                    print("== REFINED RITZ VECTORS COMPUTED ==\nTotal good eigenpairs: $(sum((ρ .< tol))) / $req\nGood eigenpairs in DoI: $good / $req\n\n")
                end
            elseif verb == 2
                print("Refined Ritz vectors not required.\n\n")
            end
        end

        if good > size(best_λ,1)
            if verb == 2
                print("Found more good eigenvalues ($(sum(good_ones)) -> $good)\n\n")
            end
            good_ones = (ρ .< tol) .&& keep.(λ) #basically nothing to recompute this: O(kl)
            best_λ = [λ[i] for i in 1:size(λ,1) if good_ones[i]] #these are global if glob is true
            best_X = hcat([X[:,i] for i in 1:size(X,2) if good_ones[i]]...)
            best_ρ = [ρ[i] for i in 1:size(ρ,1) if good_ones[i]]
            if glob; glob_λ = best_λ; glob_X = best_X; glob_ρ = best_ρ; end
        end

        if good ≥ req #if we have found enough acceptable eigenpairs
            if verb == 2
                print("$good good eigenpairs found, returning.")
            end
            good_ones = (ρ .< tol) .&& keep.(λ) #basically nothing to recompute this: O(kl)
            λ = [λ[i] for i in 1:size(λ,1) if good_ones[i]]
            X = hcat([X[:,i] for i in 1:size(X,2) if good_ones[i]]...)
            ρ = [ρ[i] for i in 1:size(ρ,1) if good_ones[i]]
            return λ,X,ρ #we're done here
        elseif m+step*l > kl_max #if another step could expand the subspace too far
            if restarts == give_up
                @warn "restart limit exceeded, not enough eigenpairs found"
                if verb > 0
                    print("⚠️Restart limit ($give_up) exceeded, giving up. (good eigenpairs: $(size(best_λ,1)))\n\n")
                end
                return best_λ,best_X,best_ρ #return at least what we have
            end
            restarts += 1
            if verb > 0
                print("== RESTART =="*"\n"^(3-verb)) #fancy way of getting the number of newlines right
            end
            #RESTART
            if !(false in keep.(λ)) #if all computed Ritz values are inside the DoI
                #restart solely according to req
                which_eigs1(λ) = [(sum(abs(λ[i]) .≤ abs.(λ)) ≤ req) for i in 1:size(λ,1)] #we actually want to remove those closest to 0 in the transformed QEP
                Q,U,H = restartBTOAR(Q,U,H,which_eigs1,verb) #do the restart
                if verb > 0
                    print("All Ritz values inside DoI, restarted according to req (=$req).\n")
                end
            elseif sum(keep.(λ)) ≥ req #if some Ritz values are outside DoI but there are enough inside
                which_eigs2(λ) = [(sum((abs(λ[i]) .≤ abs.(λ)) .&& transformed_keep(λ)) ≤ req) && transformed_keep(λ)[i] for i in 1:size(λ,1)]
                Q,U,H = restartBTOAR(Q,U,H,which_eigs2,verb) #do the restart
                if verb > 0
                    print("Enough Ritz values inside DoI, restarted according to req (=$req) and DoI.\n")
                end
            else #if there are not enough Ritz pairs inside the DoI
                which_eigs3 = transformed_keep #just restart according to DoI only (fallback method, not quite optimal)
                Q,U,H = restartBTOAR(Q,U,H,which_eigs3,verb) #do the restart
                if verb > 0
                    print("Not enough Ritz values inside DoI, restarted according to DoI (=$req).\n")
                end
            end
            MQ = M¹(Q)
            DQ = D¹(Q)
            KQ = K¹(Q)
            QᵀMQ = Q'*M¹(Q) #no way to avoid this cost really
            QᵀDQ = Q'*D¹(Q)
            QᵀKQ = Q'*K¹(Q)
            if rrv > 0 #if we use RRVs we also have to form the eigenvector matrices again
                MQᵀMQ = MQ'MQ
                MQᵀDQ = MQ'DQ
                MQᵀKQ = MQ'KQ
                DQᵀMQ = MQᵀDQ' #take advantage of past computations
                DQᵀDQ = DQ'DQ
                DQᵀKQ = DQ'KQ
                KQᵀMQ = MQᵀKQ' #
                KQᵀDQ = DQᵀKQ' #
                KQᵀKQ = KQ'KQ
            end
        else #if we are free to grow the subspace by step
            if verb > 0
                print("== CONTINUING BTOAR ALGORITHM ==\n")
            end
            Q,U,H = continueBTOAR(M⁻¹,D¹,K¹,Q,U,H,maximum([step,minimum([Int(floor((req-good)/l)),Int(floor((kl_max-m)/l))])]),l;deftol=dtol,verb=verb) #grow the subspace by as much as we can without overflowing or overdoing it
        end
        m = size(Q,2)
    end
end
########## ADD CONDITIONAL REORTHOGONALISATION AND POST-RRQR/QR FACTORISATION REORTHOGONALISATION (maybe the latter might be hard...)
########## GIVE RESTART THE FLEXIBILITY TO ALLOW FOR A DEFLATION JUST BEFORE (this should be pretty simple, right? Just change the "+ step*l" to something)
#          I think all that is required is to make the zeroes in the W and Y matrices of conforming size (so not l but m_{something}), perhaps there are other places where we should replace l with that but it
#          should be pretty simple to find what to replace l by
########## Print out Q, U and BTOAR relation residuals after each step, not just the first lol
########## When expanding with very little space left to fill (less than step) make sure to expand some amount anyway
########## Only compute residuals for Ritz values inside DoI (set all ones outside to big number?)

########## Restarting with narrow DoI (qep7747, DoI real part between -6e4 and -5e4, tol=1e-12) doesn't seem to be working correctly.
#          When there are not enough eigenvalues in DoI, but there are some, it seems to throw them all away
#          When there are enough eigenvalues in the DoI, it still strows them all away
#          This looks an awful lot like a problem with transformed_keep()
#          Also transformed_keep() should respect inv
