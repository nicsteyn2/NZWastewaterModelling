
using ProgressMeter

include("covidParticleFilterWithCumulativeInfections.jl")
include("processResults.jl")


# This function takes in a 3D array of parameters of size (T, nParams, nParamSamples) and runs the particle filter on each parameter set of size (T, nParams).
# The function then keeps a random selection of particle trajectories from each run and stores them such that (approximately) nTargetSamples are kept in total, giving samples from our marginal posterior on the hidden states.
function runFullPosterior(nTargetSamples::Int, θ::Array{Float64, 3}, Y::DataFrame, opts::Dict, label::String; outdir="outputs/fullposterior/")
    
    # Calculate the number of samples we need to keep from each 
    nParamSamples = size(θ)[3]
    nFilterSamplesPerParam = Int(round(nTargetSamples/nParamSamples))
    nTotalSamples = nFilterSamplesPerParam * nParamSamples # This will equal nTargetSamples if nTargetSamples is a multiple of nParamSamples
    
    # Pre-allocate big results
    ρ_full = zeros(nTotalSamples, opts["tMax"])
    R_full = zeros(nTotalSamples, opts["tMax"])
    I_full = zeros(nTotalSamples, opts["tMax"])
    CI_full = zeros(nTotalSamples, opts["tMax"])
    C_expected = zeros(nTotalSamples, opts["tMax"])
    C_predictive = zeros(nTotalSamples, opts["tMax"])
    W_expected = zeros(nTotalSamples, opts["tMax"])
    W_predictive = zeros(nTotalSamples, opts["tMax"])
    
    # And run the particle filter at each set of sampled parameters
    progBar = Progress(nParamSamples; dt=1, barlen = 50, desc="Sampling from full posterior...")
    Threads.@threads for ii = 1:nParamSamples
        
        # We don't want to show progress for each indvidual filter
        opts["showProgress"] = false
        
        # Extract our current parameters
        θ_in = copy(θ[:,:,ii])
        
        # Solve the model
        (ρ, R, I, CI) = initialiseHiddenStatesWithCumulativeInfections(Y, opts)
        (ρ, R, I, CI, _) = covidParticleFilterWithCumulativeInfections(ρ, R, I, CI, Y, θ_in, opts)
        
        # Set the indices to sample from and save to
        sample_inds = sample(Int.(1:opts["Nx"]), nFilterSamplesPerParam) # Randomly select the indices of the trajectories to keep
        save_inds = (1:nFilterSamplesPerParam) .+ (nFilterSamplesPerParam * (ii - 1)) # and figure out where to save them to
        
        # Save our results in the big output files
        ρ_full[save_inds,:] = ρ[sample_inds,:]
        R_full[save_inds,:] = R[sample_inds,:]
        I_full[save_inds,:] = I[sample_inds,:]
        CI_full[save_inds,:] = CI[sample_inds,:]
        
        # We also sample our observation and predictive particles at this time-step
        (expectedCases, observedCases) = sampleExpectedAndPredictiveCases(ρ, I, θ_in, opts)
        C_expected[save_inds,:] = expectedCases[sample_inds,:]
        C_predictive[save_inds,:] = observedCases[sample_inds,:]
        
        (expectedWastewater, observedWastewater) = sampleExpectedAndPredictiveWastewater(I, θ_in, opts)
        W_expected[save_inds,:] = expectedWastewater[sample_inds,:]
        W_predictive[save_inds,:] = observedWastewater[sample_inds,:]
        
        next!(progBar)
        
    end
    
    # Calculate the results
    results = processAllResults(ρ_full, R_full, I_full, CI_full, C_expected, C_predictive, W_expected, W_predictive, opts)

    # and save
    if !isdir(outdir)
        mkdir(outdir)
    end
    CSV.write(outdir * label * ".csv", results)
    CSV.write(outdir * label * "_rawdata.csv", Y)
    
    return(results)
    
end