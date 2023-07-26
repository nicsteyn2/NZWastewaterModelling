
include("src/loadNZData.jl")
include("src/getOpts.jl")
include("src/processPMMHResults.jl")
include("src/runFullPosterior.jl")



# Specify the mmhLabel to run on
function runFullPosteriorNZ_withShifts(OPTION, nParamSamples, mmhLabel, outlabel; nChains=8, miniter=500, nTargetSamples=1e6, Nx=1e5, alpha=3e9, datasource="CW")
    
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

    # Specify the shift we want to apply in the particle filter
    if OPTION == 1
        tmp = similar(opts["infToShedDist"])
        tmp[1] = 0
        tmp[2:end] = opts["infToShedDist"][1:end-1]
        opts["infToShedDist"] = tmp/sum(tmp)
    elseif OPTION == 2
        tmp = similar(opts["infToShedDist"])
        tmp[end] = 0
        tmp[1:end-1] = opts["infToShedDist"][2:end]
        opts["infToShedDist"] = tmp/sum(tmp)
    elseif OPTION == 3
        tmp = similar(opts["infToRepDist"])
        tmp[1] = 0
        tmp[2:end] = opts["infToRepDist"][1:end-1]
        opts["infToRepDist"] = tmp/sum(tmp)
    elseif OPTION == 4
        tmp = similar(opts["infToRepDist"])
        tmp[end] = 0
        tmp[1:end-1] = opts["infToRepDist"][2:end]
        opts["infToRepDist"] = tmp/sum(tmp)
    elseif OPTION == 5
        tmp = similar(opts["genTimeDist"])
        tmp[1] = 0
        tmp[2:end] = opts["genTimeDist"][1:end-1]
        opts["genTimeDist"] = tmp/sum(tmp)
    elseif OPTION == 6
        tmp = similar(opts["genTimeDist"])
        tmp[end] = 0
        tmp[1:end-1] = opts["genTimeDist"][2:end]
        opts["genTimeDist"] = tmp/sum(tmp)
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

runFullPosteriorNZ_withShifts(1, nParamSamples, in_labels, outlabel * "_lateshedding"; alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)
runFullPosteriorNZ_withShifts(2, nParamSamples, in_labels, outlabel * "_earlyshedding"; alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)
runFullPosteriorNZ_withShifts(3, nParamSamples, in_labels, outlabel * "_latereporting"; alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)
runFullPosteriorNZ_withShifts(4, nParamSamples, in_labels, outlabel * "_earlyreporting"; alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)
runFullPosteriorNZ_withShifts(3, nParamSamples, in_labels, outlabel * "_longgentime"; alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)
runFullPosteriorNZ_withShifts(4, nParamSamples, in_labels, outlabel * "_shortgentime"; alpha=3e9, Nx=Nx, nTargetSamples=nTargetSamples, miniter=miniter)

