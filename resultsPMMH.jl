
# This script runs the PMMH algorithm by time-period.
#
# Results are saved to the outputs/pmmh.nosync/model_name/ folder

using ArgParse, JSON

include("src/loadNZData.jl")
include("src/getOpts.jl")
include("src/PMMH.jl")
include("src/supportFuncs.jl")

function runPMMHNZ(period, nChains, model_name)
    
    # Specify extra wind-in period for parameter fitting
    normal_windin = 30
    extra_windin = 50
    
    # Pre-allocate period-specific inputs
    ST_DATE = nothing
    EN_DATE = nothing
    σrho_init = nothing
    σR_init = nothing
    kc_init = nothing
    kw_init = nothing
    propVars = nothing
    sigRho_prior = nothing

    # Set period-specific parameters
    if period == 1

        # Dates and wind-in period
        extra_windin = min(extra_windin, 30) # We don't have sufficient data for the full 50-day windin in period 1
        ST_DATE=Date("2022-01-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2022-04-01")

        # Values to initialise the MCMC at (more constrained than the priors to reduce the need for extended wind-in periods)
        σrho_init = (0.021, 0.032)
        σR_init = (0.11, 0.16)
        kc_init = (25, 34)
        kw_init = (1.3e-6, 1.8e-6)

        # Proposal variances
        propVars = [0.004, 0.024, 4.3, 1.4e-7]  # transition varianes on sigCAR, sigR, kc, kw

        # and the prior on σ_CAR
        sigRho_prior = Truncated(Normal(0.024, 0.008160499922185493), 0, Inf) # mean 0.024, upper quantile 0.04

    elseif period == 2

        ST_DATE=Date("2022-04-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2022-07-01")

        σrho_init = (0.0081, 0.011)
        σR_init = (0.067, 0.096)
        kc_init = (130, 170)
        kw_init = (4.2e-6, 5.4e-6)

        propVars = [0.0014, 0.018, 21, 5.5e-7]

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf) # mean 0.006, upper quantile 0.01

    elseif period == 3

        ST_DATE=Date("2022-07-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2022-10-01")

        σrho_init = (9.5e-3, 1.2e-2)
        σR_init = (3.5e-2, 5.5e-2)
        kc_init = (1.4e2, 1.9e2)
        kw_init = (4.0e-6, 5.2e-6)

        propVars = [0.0015, 0.01, 30, 6.6e-7]

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf)

    elseif period == 4

        ST_DATE=Date("2022-10-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2023-01-01")

        σrho_init = (9.7e-3, 1.2e-2)
        σR_init = (3.3e-2, 5.0e-2)
        kc_init = (1.4e2, 1.9e2)
        kw_init = (6.0e-6, 7.8e-6)

        propVars = [0.001, 0.0073, 22, 7.4e-7]

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf)

    elseif period == 5

        ST_DATE=Date("2023-01-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2023-04-01")

        σrho_init = (7.6e-3, 1.1e-2)
        σR_init = (3.4e-2, 5.5e-2)
        kc_init = (1.1e2, 1.6e2)
        kw_init = (5.8e-6, 7.7e-6)

        propVars = [0.0015, 0.0089, 22, 8.9e-7]

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf)

    end
    
    # Load data
    Y = loadNZData(ST_DATE=ST_DATE, EN_DATE=EN_DATE)
    
    # Set particle filter options (see src/getOpts.jl for details)
    opts = getOpts(Y)
    opts["Nx"] = Int(1e5)
    opts["windinPeriod"] = Int(normal_windin)
    
    # Set PMMH options (see src/getOpts.jl for details)
    mmhOpts = getMMHOpts()
    mmhOpts["nSteps"] = 5000
    mmhOpts["mmhSave"] = 10
    mmhOpts["period"] = period
    mmhOpts["paramInitBounds"][1] = σrho_init
    mmhOpts["paramInitBounds"][2] = σR_init
    mmhOpts["paramInitBounds"][3] = kc_init
    mmhOpts["paramInitBounds"][4] = kw_init
    mmhOpts["proposalVariances"] = propVars
    mmhOpts["paramPriorDists"][1] = sigRho_prior

    # Specify the filename & directory to write outputs to
    label = model_name * "_period" * string(period)
    outdir = "outputs/pmmh.nosync/" * model_name * "/"
    if !isdir(outdir)
        mkdir(outdir)
    end

    # Save the options to disk so we can remember our settings later
    # This file is also read by some of the .R scripts
    json_string = JSON.json(clean_opts_dict(mmhOpts))
    open(outdir * label * "_mmhopts.json", "w") do f
        write(f, json_string)
    end

    # And run the model
    runParticleMMHMultiThread(nChains, Y, mmhOpts, opts, label; extra_windin=extra_windin, outdir=outdir)
    
end


# # When testing locally
# period = 1
# nThreads = 1
# modelname = "localtest"

# When running on the cluster
period = parse(Int, ARGS[1])
nThreads = 8
modelname = "final"

runPMMHNZ(period, nThreads, modelname)