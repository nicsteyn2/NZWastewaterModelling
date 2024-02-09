
# This script runs the PMMH algorithm by time-period.
#
# Results are saved to the outputs/pmmh.nosync/model_name/ folder

using ArgParse, JSON

include("src/loadNZData.jl")
include("src/getOpts.jl")
include("src/PMMH.jl")
include("src/supportFuncs.jl")

function runPMMHNZ(period, nChains, model_name; datasource="CW")
    
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
        if datasource == "C"
            propVars = [0, 0.045, 6.2, 0]
        elseif datasource == "W"
            propVars = [0, 0.055, 0, 3e-7]
        end

        # and the prior on σ_CAR
        sigRho_prior = Truncated(Normal(0.024, 0.008160499922185493), 0, Inf) # mean 0.024, upper quantile 0.04

    elseif period == 2

        ST_DATE=Date("2022-04-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2022-07-01")

        σrho_init = (0.0081, 0.011)
        σR_init = (0.067, 0.096)
        kc_init = (130, 170)
        kw_init = (4.2e-6, 5.4e-6)

        # Proposal variances
        propVars = [0.0014, 0.018, 21, 5.5e-7]
        if datasource == "C"
            propVars = [0, 0.022, 32, 0]
        elseif datasource == "W"
            propVars = [0, 0.02, 0, 1.1e-6]
        end

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf) # mean 0.006, upper quantile 0.01

    elseif period == 3

        ST_DATE=Date("2022-07-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2022-10-01")

        σrho_init = (9.5e-3, 1.2e-2)
        σR_init = (3.5e-2, 5.5e-2)
        kc_init = (1.4e2, 1.9e2)
        kw_init = (4.0e-6, 5.2e-6)

        propVars = [0.0015, 0.01, 30, 6.6e-7]
        if datasource == "C"
            propVars = [0, 0.034, 67, 0]
        elseif datasource == "W"
            propVars = [0, 0.035, 0, 1.1e-6]
        end

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf)

    elseif period == 4

        ST_DATE=Date("2022-10-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2023-01-01")

        σrho_init = (9.7e-3, 1.2e-2)
        σR_init = (3.3e-2, 5.0e-2)
        kc_init = (1.4e2, 1.9e2)
        kw_init = (6.0e-6, 7.8e-6)

        propVars = [0.001, 0.0073, 22, 7.4e-7]
        if datasource == "C"
            propVars = [0, 0.024, 35, 0]
        elseif datasource == "W"
            propVars = [0, 0.026, 0, 1.9e-6]
        end

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf)

    elseif period == 5

        ST_DATE=Date("2023-01-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2023-04-01")

        σrho_init = (7.6e-3, 1.1e-2)
        σR_init = (3.4e-2, 5.5e-2)
        kc_init = (1.1e2, 1.6e2)
        kw_init = (5.8e-6, 7.7e-6)

        propVars = [0.0015, 0.0089, 22, 8.9e-7]
        if datasource == "C"
            propVars = [0, 0.028, 57, 0]
        elseif datasource == "W"
            propVars = [0, 0.030, 0, 1.5e-6]
        end

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf)

    elseif period == 6

        ST_DATE=Date("2023-04-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2023-07-01")

        σrho_init = (6.6e-3, 9.2e-3)
        σR_init = (1.5e-2, 3.0e-2)
        kc_init = (1.0e2, 1.3e2)
        kw_init = (4.1e-6, 5.2e-6)

        propVars = [0.0011, 0.0078, 14, 5.1e-7]
        if datasource == "C"
            propVars = [0, 0.02, 30, 0]
        elseif datasource == "W"
            propVars = [0, 0.045, 0, 3e-7]
        end

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf)

    elseif period == 7

        ST_DATE=Date("2023-07-01")-Day(normal_windin)-Day(extra_windin)
        EN_DATE=Date("2023-10-01")

        σrho_init = (8.5e-3, 1.0e-2)
        σR_init = (4.3e-2, 6.8e-2)
        kc_init = (1.7e2, 2.6e2)
        kw_init = (6.5e-6, 8.8e-6)

        propVars = [0.001, 0.01, 38, 1.0e-6]
        if datasource == "C"
            propVars = [0, 0.02, 60, 0]
        elseif datasource == "W"
            propVars = [0, 0.045, 0, 3e-7]
        end

        sigRho_prior = Truncated(Normal(0.006, 0.0020401247804969313), 0, Inf)

    end
    
    # Load data
    Y = loadNZData(ST_DATE=ST_DATE, EN_DATE=EN_DATE)
    
    # Set particle filter options (see src/getOpts.jl for details)
    opts = getOpts(Y)
    opts["Nx"] = Int(1e5)
    opts["windinPeriod"] = Int(normal_windin)
    opts["datasource"] = datasource
    
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

    if datasource == "C"
        mmhOpts["estimated"] = [false, true, true, false]
        mmhOpts["defaultValues"] = [0, NaN, NaN, 0]
        opts["Nx"] = Int(1e4)
    elseif datasource == "W"
        mmhOpts["estimated"] = [false, true, false, true]
        mmhOpts["defaultValues"] = [0, NaN, 0, NaN]
        opts["Nx"] = Int(1e4)
    end

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
# nThreads = 3
# modelname = "localtest"

# When running on the cluster
period = parse(Int, ARGS[1])
nThreads = 8
modelname = "final"

# Full model
runPMMHNZ(period, nThreads, modelname)


# Cases/wastewater only

# # Check variance of log-lik estimates
# Y = loadNZData(ST_DATE=Date("2022-01-01")-Day(30)-Day(50), EN_DATE=Date("2022-04-01"))
# opts = getOpts(Y)
# opts["Nx"] = Int(1e4)
# opts["windinPeriod"] = Int(30)
# opts["datasource"] = "W"
# (v, logev) = estimateVarianceOfLogLikEsts(100, [0, 0.13, 0, 1.5e-6], 50, Y, opts)

# # Run PMMH
# runPMMHNZ(period, nThreads, modelname; datasource="W")