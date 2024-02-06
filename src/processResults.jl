
include("supportFuncs.jl")
include("predictObservations.jl")

# Functions to process particle filter results

# If we haven't already calculated the expected/observed cases and wastewater particles
function processAllResults(ρ, R, I, θ, opts)
    
    println("Sampling expected and predictive cases...")
    (C_expected, C_predictive) = sampleExpectedAndPredictiveCases(ρ, I, θ, opts)
    println("Sampling expected and predictive wastewater...")
    (W_expected, W_predictive) = sampleExpectedAndPredictiveWastewater(I, θ, opts)
    return(processAllResults(ρ, R, I, C_expected, C_predictive, W_expected, W_predictive, opts))
    
end


# This version uses provided predictive cases and wastewater samples
function processAllResults(ρ, R, I, C_expected::Matrix, C_predictive::Matrix, W_expected::Matrix, W_predictive::Matrix, opts::Dict)
    
    # Process results
    Tρ = calculateQuantiles(ρ, opts; var_name="ρ")
    TR = calculateQuantiles(R, opts; var_name="R")
    TI = calculateQuantiles(I, opts; var_name="I")
    Tc_expected = calculateQuantiles(C_expected, opts; var_name="Expected cases")
    Tc_predictive = calculateQuantiles(C_predictive, opts; var_name="Predictive cases")
    Tw_expected = calculateQuantiles(W_expected, opts; var_name="Expected wastewater")
    Tw_predictive = calculateQuantiles(W_predictive, opts; var_name="Predictive wastewater")
    TpRLessThanOne = calculateProbRLessThanOne(R, opts)
        
    # Tidy up and join
    select!(Tc_predictive, :lower, :upper, :date)
    rename!(Tc_predictive, :lower => :predictive_lower, :upper => :predictive_upper)
    Tc = outerjoin(Tc_expected, Tc_predictive, on=:date)
    
    select!(Tw_predictive, :lower, :upper, :date)
    rename!(Tw_predictive, :lower => :predictive_lower, :upper => :predictive_upper)
    Tw = outerjoin(Tw_expected, Tw_predictive, on=:date)
    
    # Append variable names
    Tρ.variable = "CARt"
    TR.variable = "Rt"
    TI.variable = "It"
    Tc.variable = "Cases"
    Tw.variable = "Wastewater"
    TpRLessThanOne.variable = "pRLessThanOne"
    
    # Append missing columns
    Tρ.predictive_lower = NaN
    Tρ.predictive_upper = NaN
    TR.predictive_lower = NaN
    TR.predictive_upper = NaN
    TI.predictive_lower = NaN
    TI.predictive_upper = NaN
    TpRLessThanOne.predictive_lower = NaN
    TpRLessThanOne.predictive_upper = NaN
    
    # Stack everything together
    Tout = [Tρ; TR; TI; Tc; Tw; TpRLessThanOne]
    return(Tout)
    
end



# This calculates the quantiles of our particles
function calculateQuantiles(X, opts; var_name="generic variable")
    
    Tout = DataFrame()
    Tout.mean = zeros(opts["tMax"])
    Tout.median = zeros(opts["tMax"])
    Tout.lower = zeros(opts["tMax"])
    Tout.upper = zeros(opts["tMax"])
    
    progBar = Progress(opts["tMax"]; dt=1, barlen = 50, desc="Calculating results for: " * var_name)
    for tt in 1:opts["tMax"]
        x = X[:,tt]
        if !isnan(x[1])
            Tout.mean[tt] = mean(x)
            (Tout.median[tt], Tout.lower[tt], Tout.upper[tt]) = quantile(x, [0.5, 0.025, 0.975])
        else
            Tout.mean[tt] = NaN
            Tout.median[tt] = NaN
            Tout.lower[tt] = NaN
            Tout.upper[tt] =  NaN
        end
        next!(progBar)
    end
    
    Tout.t = 1:opts["tMax"]
    Tout.date = opts["dates"]
    
    return Tout
    
end


# And also calcualtes the probability that the particles are less than 1 (useful for P(Rt<1))
function calculateProbRLessThanOne(X, opts)
    
    println("Calculating prob R less than 1...")
    
    Tout = DataFrame()    
    Tout.mean = vec(mean(X.<=1, dims=1))
    Tout.t = 1:opts["tMax"]
    Tout.date = opts["dates"]
    Tout.median = NaN
    Tout.lower = NaN
    Tout.upper = NaN

    return(Tout)
    
end




# And additional versions for when we also pass in cumulative infections
function processAllResults(ρ, R, I, CI, θ, opts)
    
    println("Sampling expected and predictive cases...")
    (C_expected, C_predictive) = sampleExpectedAndPredictiveCases(ρ, I, θ, opts)
    println("Sampling expected and predictive wastewater...")
    (W_expected, W_predictive) = sampleExpectedAndPredictiveWastewater(I, θ, opts)
    return(processAllResults(ρ, R, I, CI, C_expected, C_predictive, W_expected, W_predictive, opts))
    
end


# This version uses provided predictive cases and wastewater samples
function processAllResults(ρ, R, I, CI, C_expected::Matrix, C_predictive::Matrix, W_expected::Matrix, W_predictive::Matrix, opts::Dict)
    
    # Process results
    Tρ = calculateQuantiles(ρ, opts; var_name="ρ")
    TR = calculateQuantiles(R, opts; var_name="R")
    TI = calculateQuantiles(I, opts; var_name="I")
    TCI = calculateQuantiles(CI, opts; var_name="CI")
    Tc_expected = calculateQuantiles(C_expected, opts; var_name="Expected cases")
    Tc_predictive = calculateQuantiles(C_predictive, opts; var_name="Predictive cases")
    Tw_expected = calculateQuantiles(W_expected, opts; var_name="Expected wastewater")
    Tw_predictive = calculateQuantiles(W_predictive, opts; var_name="Predictive wastewater")
    TpRLessThanOne = calculateProbRLessThanOne(R, opts)
        
    # Tidy up and join
    select!(Tc_predictive, :lower, :upper, :date)
    rename!(Tc_predictive, :lower => :predictive_lower, :upper => :predictive_upper)
    Tc = outerjoin(Tc_expected, Tc_predictive, on=:date)
    
    select!(Tw_predictive, :lower, :upper, :date)
    rename!(Tw_predictive, :lower => :predictive_lower, :upper => :predictive_upper)
    Tw = outerjoin(Tw_expected, Tw_predictive, on=:date)
    
    # Append variable names
    Tρ.variable = "CARt"
    TR.variable = "Rt"
    TI.variable = "It"
    TCI.variable = "CI"
    Tc.variable = "Cases"
    Tw.variable = "Wastewater"
    TpRLessThanOne.variable = "pRLessThanOne"
    
    # Append missing columns
    Tρ.predictive_lower = NaN
    Tρ.predictive_upper = NaN
    TR.predictive_lower = NaN
    TR.predictive_upper = NaN
    TI.predictive_lower = NaN
    TI.predictive_upper = NaN
    TCI.predictive_lower = NaN
    TCI.predictive_upper = NaN
    TpRLessThanOne.predictive_lower = NaN
    TpRLessThanOne.predictive_upper = NaN
    
    # Stack everything together
    Tout = [Tρ; TR; TI; TCI; Tc; Tw; TpRLessThanOne]
    return(Tout)
    
end
