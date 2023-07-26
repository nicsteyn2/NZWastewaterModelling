
using Dates, DataFrames

include("src/runCurvewiseEstimates.jl")
include("src/processPMMHResults.jl")


function resultsCurvewiseEstimatesNZ(ST_DATE, EN_DATE, nParamSamples, mmhLabel, outLabel; nChains=8, miniter=500, nTargetSamples=1e6, Nx=1e5, alpha=3e9, datasource="CW", h=60, windin=60)
    
    # Load and tidy the PMMH samples
    all_samples = fetchPMMHOutputsForNZ(mmhLabel, nChains, miniter)
    
    # Load data and define start and end dates
    Y = loadNZData(ST_DATE  = ST_DATE - Day(windin), EN_DATE = EN_DATE)
    st_dates = [Date("2021-01-01"), Date("2022-04-01"), Date("2022-07-01"), Date("2022-10-01"), Date("2023-01-01")]
    en_dates = [Date("2022-03-31"), Date("2022-06-30"), Date("2022-09-30"), Date("2022-12-31"), Date("2024-01-01")]
    
    # Set particle filter opts
    opts = getOpts(Y)
    opts["Nx"] = Int(Nx)
    opts["datasource"] = datasource
    opts["alpha"] = alpha
    opts["resamplingWindow"] = h
    opts["windinPeriod"] = windin
    
    # Set MMH opts
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
    
    # Print out important info
    first_valid_date = max(Y.date[length(Y.date) - h], ST_DATE + Day(30))
    last_valid_date = Y.date[length(Y.date) - 30]
    mid_valid_date = Y.date[length(Y.date) - Int(round((h + 30)/2))]
    println()
    println("Check that the range of possible dates lie within the following values:")
    println("The first valid date is " * string(first_valid_date))
    println("The middle valid date is " * string(mid_valid_date))
    println("The last valid date is " * string(last_valid_date))
    println()
    
    # And run the model
    results = runCurvewiseEstimates(Int(nTargetSamples), θ_all, Y, opts, outLabel)
    
    return(results)
    
end


# Set options
mmhlabels = "final"
nParamSamples = 100
nTargetSamples = 2e6
miniter = 100
Nx = 1e5
windin = 60
h = 70
datasource = "CW"

# Run Feb 2022 peak
ST_DATE=Date("2022-01-07")
EN_DATE=Date("2022-04-17")
outlabel = "Feb2022PeakR" * "_" * "final"

println("Running curvewise estimator for Feb 2022 peak. We expect peak R to occur around 20 Feb, with a heavy lower tail.")
println("This code can also be used for first time R < 1, which we expect to occur around 28 Feb.")
samples = resultsCurvewiseEstimatesNZ(ST_DATE, EN_DATE, nParamSamples, mmhlabels, outlabel; nChains=8, miniter=miniter, nTargetSamples=nTargetSamples, Nx=Nx, h=h, windin=windin, datasource=datasource)

# Run July 2022 peak
ST_DATE=Date("2022-04-27")
EN_DATE=Date("2022-08-27")
outlabel = "July2022PeakCases" * "_" * "final"

samples = resultsCurvewiseEstimatesNZ(ST_DATE, EN_DATE, nParamSamples, mmhlabels, outlabel; nChains=8, miniter=miniter, nTargetSamples=nTargetSamples, Nx=Nx, h=h, windin=windin, datasource=datasource)

