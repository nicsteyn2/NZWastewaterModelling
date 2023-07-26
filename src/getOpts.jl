
using Distributions, CSV, DataFrames

include("supportFuncs.jl")

# The getOpts() and getMMHOpts() specify dictionaries that contain settings for the particle filter and PMMH respectively. 
# Note that sigRho is the standard deviation of the random walk on CARt. This is called σ_CAR in the paper.

# Define particle filter options dictionary
function getOpts(Y)
    
    # Pre-allocate
    opts = Dict()
    
    # User choices
    opts["showProgress"] = true

    # Model hyperparameters and other options
    opts["Nx"] = Int(1e4) # Number of particles
    opts["tMax"] = maximum(Y.t) # Duration
    opts["resamplingWindow"] = 30 # h - the length of the resampling window
    opts["windinPeriod"] = 30 # Number of days to use in wind-in, should be ≥ opts["resamplingWindow]
    opts["dates"] = Y.date # A vector of dates
    opts["datasource"] = "CW" # Do we filter on cases only (C), wastewater only (W), or both (CW)?

    # Calculate the generation time distribution (mean and sd from Vattiatio et al DOI: 10.1016/j.epidem.2022.100657)
    gen_time_mean = 3.3
    gen_time_sd = 1.3
    gen_time_shape = (gen_time_mean^2)/(gen_time_sd^2)
    gen_time_scale = (gen_time_sd^2)/gen_time_mean
    gen_time_raw = pdf.(Gamma(gen_time_shape, gen_time_scale), 0:100)
    opts["genTimeDist"] = gen_time_raw/sum(gen_time_raw)

    # Load incubation, reporting, and shedding distributions
    shedding_raw = DataFrame(CSV.File("data/dist_ww_shedding.csv")) # Source: Hewitt et al. (2022) DOI: 10.1016/j.watres.2021.118032)
    shedding = shedding_raw.Shedding/sum(shedding_raw.Shedding)

    incubation_raw = pdf.(Weibull(1.5, 3.2), 0:1:10) # Source: Jantien et al (2022) DOI: 10.2807/1560-7917
    incubation = incubation_raw/sum(incubation_raw)

    reporting_raw = DataFrame(CSV.File("data/dist_esr_onset.csv")) # Source: NZ case data from ESR
    reporting = reporting_raw.Count/sum(reporting_raw.Count)
    
    # Specify delay distributions
    opts["infToRepDist"] = conv(incubation, reporting) # conv() is a custom function defined in src/supportFuncs.jl
    opts["maxInfToRepTime"] = length(opts["infToRepDist"])
    
    infToShedDist_tmp = conv(incubation, shedding)[4:end]
    opts["infToShedDist"] = infToShedDist_tmp/sum(infToShedDist_tmp)
    opts["maxSheddingTime"] = length(opts["infToShedDist"])
    
    # Other epidemiological parameters
    opts["alpha"] = 3e9
    opts["nPop"] = 5.15e6
    opts["popnSizeForPredictiveWastewater"] = 1e6
    
    return(opts)
    
end



# Define PMMH options dictionary
function getMMHOpts()
    
    # Initialise empty dictionary
    mmhopts = Dict()
    
    # MMH settings
    mmhopts["nSteps"] = 10000
    mmhopts["resampleLogLikFreq"] = 5
    mmhopts["mmhSave"] = 10

    # Define default parameter values and transition distributions, these must be specified at runtime
    mmhopts["proposalVariances"] = NaN
    
    # Define parameters
    mmhopts["nParams"] = 4
    mmhopts["paramNames"] = ["sigRho", "sigR", "kc", "kw"]
    mmhopts["paramNamesWithExtra"] = ["sigRho", "sigR", "kc", "kw", "did_accept", "loglik"]
    mmhopts["estimated"] = [true, true, true, true]
    
    # Define parameter bounds
    mmhopts["paramBounds"] = Dict()
    mmhopts["paramBounds"][1] = (0, 1.0) # In practice this is handled by the non-Uniform prior
    mmhopts["paramBounds"][2] = (0, 0.5) # If σ_R = 0.5 then we expect it to change by a value of 0.4 each day (around R = 1), which is very big
    mmhopts["paramBounds"][3] = (0, 400) # An upper limit of 400 places a lower bound on the observation distribution variance that is 0.25% greater than the minimum possible. Values of k_c > 400 have very little practical significance.
    mmhopts["paramBounds"][4] = (0, 0.02) # In practice values are of the order 10^-6, so this is largely irrelevant
    
    # Define prior dists on parameters
    mmhopts["paramPriorDists"] = Dict()
    mmhopts["paramPriorDists"][1] = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf) # mean 0.006, upper quantile 0.01
    # mmhopts["paramPriorDists"][1] = Truncated(Normal(0.024, 0.008160499922185493), 0, Inf) # we use this in period 1 when applying to NZ
    mmhopts["paramPriorDists"][2] = Uniform(mmhopts["paramBounds"][2][1], mmhopts["paramBounds"][2][2])
    mmhopts["paramPriorDists"][3] = Uniform(mmhopts["paramBounds"][3][1], mmhopts["paramBounds"][3][2])
    mmhopts["paramPriorDists"][4] = Uniform(mmhopts["paramBounds"][4][1], mmhopts["paramBounds"][4][2])
    
    # Define proposal dists
    mmhopts["paramProposalDist"] = Dict()
    mmhopts["paramProposalDist"][1] = (x, mmhopts) -> Truncated(Normal(x, mmhopts["proposalVariances"][1]), mmhopts["paramBounds"][1][1], mmhopts["paramBounds"][1][2])
    mmhopts["paramProposalDist"][2] = (x, mmhopts) -> Truncated(Normal(x, mmhopts["proposalVariances"][2]), mmhopts["paramBounds"][2][1], mmhopts["paramBounds"][2][2])
    mmhopts["paramProposalDist"][3] = (x, mmhopts) -> Truncated(Normal(x, mmhopts["proposalVariances"][3]), mmhopts["paramBounds"][3][1], mmhopts["paramBounds"][3][2])
    mmhopts["paramProposalDist"][4] = (x, mmhopts) -> Truncated(Normal(x, mmhopts["proposalVariances"][4]), mmhopts["paramBounds"][4][1], mmhopts["paramBounds"][4][2])
    
    # Define initial dists (so we can start by sampling from these instead of from the priors to save time)
    # if these are NOT set, the algorithm defaults to sampling from the prior distributions
    # if these are set, the algorithm samples independently uniformly in the given ranges
    mmhopts["paramInitBounds"] = Dict()
    mmhopts["paramInitBounds"][1] = missing
    mmhopts["paramInitBounds"][2] = missing
    mmhopts["paramInitBounds"][3] = missing
    mmhopts["paramInitBounds"][4] = missing

    return(mmhopts)
    
end
