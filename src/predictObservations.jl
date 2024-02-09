
using Distributions

# If we are predicting cases at a single value of θ
function sampleExpectedAndPredictiveCases(ρ, I, θ::Vector, opts)
    
    θ_matrix = repeat(θ', opts["tMax"], 1)
    return(sampleExpectedAndPredictiveCases(ρ, I, θ_matrix, opts))
    
end

# If we are predicting cases where θ changes
function sampleExpectedAndPredictiveCases(ρ, I, θ::Matrix, opts)
    
    # Calculate the expected number of cases under perfect case-ascertainment
    nrows, ncols = size(I)
    X = zeros(nrows, ncols)
    for jj = (opts["windinPeriod"]+1):ncols
        st_ind = max(1, jj-opts["maxInfToRepTime"]+1)
        en_ind = min(jj, opts["maxInfToRepTime"])
        X[:,jj] = sum(I[:,jj:-1:st_ind] .* opts["infToRepDist"][1:1:en_ind]', dims=2)
    end
    
    # and calculate the expected number of cases given the CAR at the time
    expectedCasesParticles = ρ .* X
    expectedCasesParticles[:,(1:opts["windinPeriod"])] .= 1 # Just store some blank value prior to the wind-in
    
    # and then we use expectedCasesParticles as the mean of the observation distribution and sample from this
    if minimum(θ[:,3]) > 0
        r = repeat(θ[:,3]', size(expectedCasesParticles)[1], 1)
        p = r ./(r .+ expectedCasesParticles)
        predictiveCasesParticles = rand.(NegativeBinomial.(r, p))
    else
        predictiveCasesParticles = zeros(size(expectedCasesParticles))
    end
    
    # and return our samples
    return(expectedCasesParticles, predictiveCasesParticles)
    
end



# If we are predicting wastewater at a single value of θ
function sampleExpectedAndPredictiveWastewater(I, θ::Vector, opts)
    
    θ_matrix = repeat(θ', opts["tMax"], 1)
    return(sampleExpectedAndPredictiveWastewater(I, θ_matrix, opts))
    
end

# If we are predicting wastewater at possibly varying θ
function sampleExpectedAndPredictiveWastewater(I, θ::Matrix, opts)
    
    # Calculate the expected wastewater observations
    nrows, ncols = size(I)
    X = zeros(nrows, ncols)
    for jj = 30:ncols
        st_ind = max(1, jj - opts["maxSheddingTime"]+1)
        en_ind = min(jj, opts["maxSheddingTime"])
        X[:,jj] = sum(opts["alpha"] .* I[:,jj:-1:st_ind] .* opts["infToShedDist"][1:1:en_ind]', dims=2)
    end
    expectedWastewater = X./opts["nPop"]
    expectedWastewater[:,(1:opts["windinPeriod"])] .= 1
    
    # This comes in handy if we replace "popnSizeForPredictiveWastewater" with the actual population size sampled, which can sometimes be zero.
    function safe_gamma_rand(shape, scale)
        if ismissing(shape) || ismissing(scale) || scale == 0 || shape == 0
            return(NaN)
        else
            return(rand.(Gamma(shape, scale)))
        end
    end
    
    shape = repeat(θ[:,4]', size(expectedWastewater)[1], 1) .* opts["popnSizeForPredictiveWastewater"]
    scale = expectedWastewater./shape
    predictiveWastewater = safe_gamma_rand.(shape, scale)
    
    return(expectedWastewater, predictiveWastewater)
    
    
end


