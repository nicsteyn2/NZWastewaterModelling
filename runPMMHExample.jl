
using JSON

include("src/loadNZData.jl")
include("src/getOpts.jl")
include("src/PMMH.jl")
include("src/supportFuncs.jl")

function runPMMHExample(ST_DATE, EN_DATE, model_name; Nx=1e5, alpha=3e9, h=30, windin=30, nThreads=2)
    
    # Specify extra wind-in period for parameter fitting
    extra_windin = 50

    # Load data
    Y = loadNZData(ST_DATE=ST_DATE - Day(windin) - Day(extra_windin), EN_DATE=EN_DATE)

    # Specify the directory & filename to write outputs to and ensure the directory exists
    label = model_name * "/" * model_name
    outdir = "outputs/pmmh.nosync/" * model_name * "/"
    if !isdir(outdir)
        mkdir(outdir)
    end

    # Set particle filter options (see src/getOpts.jl for details)
    opts = getOpts(Y)
    opts["Nx"] = Int(Nx)
    opts["windinPeriod"] = Int(windin)
    opts["alpha"] = alpha
    opts["resamplingWindow"] = h
    
    # Set PMMH options (see src/getOpts.jl for details)
    mmhOpts = getMMHOpts()
    mmhOpts["nSteps"] = 10000
    mmhOpts["mmhSave"] = 10
    mmhOpts["paramInitBounds"][1] = (0.02, 0.07)
    mmhOpts["paramInitBounds"][2] = (0.07, 0.23)
    mmhOpts["paramInitBounds"][3] = (21, 66)
    mmhOpts["paramInitBounds"][4] = (4.3e-7, 7.6e-7)
    mmhOpts["proposalVariances"] = [0.008, 0.024, 3.7, 0.8e-7] # transition varianes on sigRho, sigR, kc, kw
    
    # Save the options to disk so we can remember our settings later
    # This file is also read by some of the .R scripts
    json_string = JSON.json(clean_opts_dict(mmhOpts))
    open("outputs/pmmh.nosync/" * label * "_mmhopts.json", "w") do f
        write(f, json_string)
    end
    
    # And run the model
    runParticleMMHMultiThread(nThreads, Y, mmhOpts, opts, label; extra_windin=extra_windin)
    
end

# When testing locally
ST_DATE = Date("2022-01-01")
EN_DATE = Date("2022-03-31")
modelname = "example"
runPMMHExample(ST_DATE, EN_DATE, modelname; nThreads=2)