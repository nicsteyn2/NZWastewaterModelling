

include("src/loadNZData.jl")
include("src/getOpts.jl")
include("src/processPMMHResults.jl")
include("src/runFullPosterior.jl")



# Specify the mmhLabel to run on
function runFullPosteriorNZUpdated(nParamSamples, mmhLabel, outlabel; nChains=8, miniter=500, nTargetSamples=1e6, Nx=1e5, alpha=3e9, datasource="CW")
    
    # Load and tidy the PMMH samples
    all_samples = fetchPMMHOutputsForNZUpdated(mmhLabel, nChains, miniter)
    
    # Load data and define start and end dates for all our periods
    Y = loadNZData(ST_DATE=Date("2022-01-01")-Day(50), EN_DATE=Date("2023-09-30"))
    st_dates = [Date("2021-01-01"), Date("2022-04-01"), Date("2022-07-01"), Date("2022-10-01"), Date("2023-01-01"), Date("2023-04-01"), Date("2023-07-01")]
    en_dates = [Date("2022-03-31"), Date("2022-06-30"), Date("2022-09-30"), Date("2022-12-31"), Date("2024-01-01"), Date("2023-06-30"), Date("2023-09-30")]
    
    # Set particle filter options
    opts = getOpts(Y)
    opts["Nx"] = Int(Nx)
    opts["alpha"] = alpha
    opts["datasource"] = datasource

    # Get mmhopts, this is just to access mmhopts["nParams"] and mmhopts["paramNames"]
    mmhopts = getMMHOpts()

    # We start by sampling from each period
    θ_new = zeros(7, mmhopts["nParams"], nParamSamples)
    for period = 1:7
        θ_new[period, :, :] = sampleParamsFromMMH(all_samples[all_samples.period.==period,:], nParamSamples, mmhopts["paramNames"])'
    end

    # And then expand into the full matrix
    θ_all = zeros(opts["tMax"], mmhopts["nParams"], nParamSamples)
    for ii = 1:nParamSamples
        θ_all[:,:,ii] = makeParamMatrix(θ_new[:,:,ii], Y.date, st_dates, en_dates)
    end
    
    # And we can finally run the model!
    results = runFullPosterior(Int(nTargetSamples), θ_all, Y, opts, outlabel)

    return(results)
    
end


function fetchPMMHOutputsForNZUpdated(mmhLabel::Vector{String}, nChains, miniter)
        
    # Load PMMH samples for all periods
    all_samples = DataFrame()
    for period = 1:length(mmhLabel)
        current_mmh = mmhLabel[period]
        indir = "outputs/pmmh.nosync/" * current_mmh * "/"
        current_label = current_mmh * "_period" * string(period)
        current_samples = loadPMMHOutputs(current_label, nChains; indir=indir)
        current_samples.period .= period
        all_samples = [all_samples; current_samples]
    end
    
    # Remove missing data and cut-off the wind-in period
    all_samples = all_samples[all_samples[:,1] .!= 0, :]
    all_samples = all_samples[all_samples.iter .>= miniter, :]
    
    return(all_samples)
    
end


# Choose inputs
in_labels = ["final", "final", "final", "final", "final", "updates", "updates"]
outlabel = "updates"
nParamSamples = 100
nTargetSamples = 2e6
Nx = 1e5
miniter = 100

# Run the full posterior at our standard vaulue of α = 3e9...
runFullPosteriorNZUpdated(nParamSamples, in_labels, outlabel; alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)
