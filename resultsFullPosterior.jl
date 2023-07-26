
include("src/loadNZData.jl")
include("src/getOpts.jl")
include("src/processPMMHResults.jl")
include("src/runFullPosterior.jl")



# Specify the mmhLabel to run on
function runFullPosteriorNZ(nParamSamples, mmhLabel, outlabel; nChains=8, miniter=500, nTargetSamples=1e6, Nx=1e5, alpha=3e9, datasource="CW")
    
    # Load and tidy the PMMH samples
    all_samples = fetchPMMHOutputsForNZ(mmhLabel, nChains, miniter)
    
    # Load data and define start and end dates for all our periods
    Y = loadNZData(ST_DATE=Date("2022-01-01")-Day(50), EN_DATE=Date("2023-03-31"))
    st_dates = [Date("2021-01-01"), Date("2022-04-01"), Date("2022-07-01"), Date("2022-10-01"), Date("2023-01-01")]
    en_dates = [Date("2022-03-31"), Date("2022-06-30"), Date("2022-09-30"), Date("2022-12-31"), Date("2024-01-01")]
    
    # Set particle filter options
    opts = getOpts(Y)
    opts["Nx"] = Int(Nx)
    opts["alpha"] = alpha
    opts["datasource"] = datasource

    # Get mmhopts, this is just to access mmhopts["nParams"] and mmhopts["paramNames"]
    mmhopts = getMMHOpts()

    # We start by sampling from each period
    θ_new = zeros(5, mmhopts["nParams"], nParamSamples)
    for period = 1:5
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


# Choose inputs
in_labels = "final"
outlabel = "final"
nParamSamples = 100
nTargetSamples = 2e6
Nx = 1e5
miniter = 100

# Run the full posterior at our standard vaulue of α = 3e9...
runFullPosteriorNZ(nParamSamples, in_labels, outlabel; alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)

# ... and now at α = 2e9
runFullPosteriorNZ(nParamSamples, in_labels, outlabel * string("_alpha2e9"); alpha=2e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)

# ... and now at α = 4e9
runFullPosteriorNZ(nParamSamples, in_labels, outlabel * string("_alpha4e9"); alpha=4e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)

# # ... and now at α = 3e9 but on cases only
# runFullPosteriorNZ(nParamSamples, in_labels, outlabel * "_casesonly"; datasource="C", alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)

# # ... and now at α = 3e9 but on wastewater only
# runFullPosteriorNZ(nParamSamples, in_labels, outlabel * "_wastewateronly"; datasource="C", alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)