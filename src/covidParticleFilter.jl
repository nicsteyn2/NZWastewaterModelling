
using Distributions, ProgressMeter


# This is probably the most important script in this repo as it implements the core particle filter.


# Wrapper function for the particle filter when we are passing in a single vector for θ.
# This just repeats the θ vector T times (once for each time-step) and then calls the function below
function covidParticleFilter(ρ, R, I, Y, θ::Vector, opts)
    
    θ_matrix = repeat(θ', opts["tMax"], 1)
    return(covidParticleFilter(ρ, R, I, Y, θ_matrix, opts))
    
end


# The primary particle filter!
function covidParticleFilter(ρ, R, I, Y, θ::Matrix, opts)
    
    # Note: θ[tt,:] = [sigRho, sigR, k_c, k_w]
    
    # Pre-allocate our stepwise evidence vector (used in likelihood calculations)
    STEPWISEEV = zeros(opts["tMax"])
    
    progBar = Progress(opts["tMax"]-opts["windinPeriod"]; dt=1, barlen = 50, desc="Running SIR particle filter...", enabled=opts["showProgress"])
    for tt = (opts["windinPeriod"]+1):opts["tMax"]
        
        # Calculate convolutions of incidence with various distributions of interest. These functions are included below.
        infectionPressure = calculateInfectionPressure(tt, I, opts)
        genomePressurePerCapita = calculateGenomePressurePerCapita(tt, I, opts)
        reportingPressure = calculateReportingPressure(tt, I, opts)
        
        # Project current particles
        ρ[:,tt] = rand.(Truncated.(Normal.(ρ[:,tt-1], θ[tt,1]), 0, 1))
        R[:,tt] = rand.(Truncated.(Normal.(R[:,tt-1], θ[tt,2]*R[:,tt-1]), 0, Inf))
        try
            I[:,tt] = rand.(Poisson.(infectionPressure .* R[:,tt]))
        catch e
            println("Poisson sampling failed.")
            STEPWISEEV .= 0
            return(NaN, NaN, NaN, STEPWISEEV)
        end
        
        # Calculate case weights
        if Y.casesDataIsValid[tt] && opts["datasource"] != "W"
            μc = (ρ[:,tt] .* reportingPressure) # This is the expected value of cases
            r = θ[tt,3]
            p = r./(r .+ μc)
            W_cases = pdf.(NegativeBinomial.(r, p), Y.nCases[tt])
        else
            W_cases = ones(opts["Nx"])
        end
        
        # Calculate wastewater weights
        if Y.wwDataIsValid[tt] && opts["datasource"] != "C"
            shape = θ[tt,4] * Y.nPopInCatchment[tt]
            scale = genomePressurePerCapita./shape
            W_ww = pdf.(Gamma.(shape, scale), Y.nGenomeCopies[tt])
        else
            W_ww = ones(opts["Nx"])
        end  
        
        # Calculate overall weights
        if opts["datasource"] == "CW"
            W = W_cases .* W_ww
        elseif opts["datasource"] == "C"
            W = W_cases
        elseif opts["datasource"] == "W"
            W = W_ww
        else
            error("Invalid choice of method. Options are: cases only (C), wastewater only (W), or both (CW).")
        end
        
        # Calculate the stepwise eveidnce
        STEPWISEEV[tt] = mean(W)

        # And attempt to resample particles. If this fails then either θ is a very poor choice or Nx is too small.
        try
            resampleParticles!(tt, W, ρ, R, I, θ[tt,:], opts)
        catch e
            println("θ = " * string(θ[tt,:]))
            println("tt = " * string(tt))
            println(sum(W))
            @warn("An error occured. It is likely that the sum of particle weights is zero (or very close to zero). Returning NaN with log-likelihood = 0.")
            return(NaN, NaN, NaN, zeros(opts["Nx"]))
        end
        
        # Update progress bar
        next!(progBar)
        
    end
    
    return(ρ, R, I, STEPWISEEV)
    
end




# -------------- Function to initialise the hidden states --------------

function initialiseHiddenStates(Y, opts)
    
    # Initialise hidden-states
    ρ = zeros(opts["Nx"], opts["tMax"])
    R = zeros(opts["Nx"], opts["tMax"])
    I = zeros(opts["Nx"], opts["tMax"])
    
    # And just fill with random numbers
    for tt = 1:opts["windinPeriod"]
        ρ[:,tt] = rand(Uniform(0.05, 0.95), opts["Nx"])
        R[:,tt] = 3*rand(opts["Nx"])
        I[:,tt] = round.(rand(Uniform(Y.nCases[tt]/0.95, Y.nCases[tt]/0.05), opts["Nx"]))
    end

    # If we aren't using both datasources, then set ρ = 0.5
    if opts["datasource"] != "CW"
        ρ[:,1:opts["windinPeriod"]] .= 0.5
    end
    
    return(ρ, R, I)
    
end


# -------------- Attempt resampling of particles --------------

function resampleParticles!(tt, W, ρ, R, I, θ, opts)
    
    # Pre-allocate inds variable
    inds = wsample(1:opts["Nx"], vec(W), opts["Nx"])
    st_time = max(1, tt - opts["resamplingWindow"])
    ρ[:,st_time:tt] = ρ[inds,st_time:tt]
    R[:,st_time:tt] = R[inds,st_time:tt]
    I[:,st_time:tt] = I[inds,st_time:tt]
    
end


# -------------- Convolutions of incidence with delay/reporting/shedding distributions --------------

# Calculate infection pressure
function calculateInfectionPressure(tt, I, opts)
    st_ind = max(1, tt - opts["resamplingWindow"] + 1)
    en_ind = min(tt, opts["resamplingWindow"])
    infectionPressure = sum(I[:,tt:-1:st_ind] .* opts["genTimeDist"][1:1:en_ind]', dims=2)
    return(infectionPressure)
end


# Calculate shedding
function calculateGenomePressurePerCapita(tt, I, opts)
    st_ind = max(1, tt - opts["maxSheddingTime"] + 1)
    en_ind = min(tt, opts["maxSheddingTime"])
    genomePressure = sum(opts["alpha"] .* I[:, tt:-1:st_ind] .* opts["infToShedDist"][1:1:en_ind]', dims=2)
    genomePressurePerCapita = genomePressure/opts["nPop"]
    return(genomePressurePerCapita)
end

# Calculate expected reported cases (under perfect case ascertainment) accounting for lags
function calculateReportingPressure(tt, I, opts)
    st_ind = max(1, tt - opts["maxInfToRepTime"] + 1)
    en_ind = min(tt, opts["maxInfToRepTime"])
    reportingPressure = sum(I[:,tt:-1:st_ind] .* opts["infToRepDist"][1:1:en_ind]', dims=2)
    return(reportingPressure)
end