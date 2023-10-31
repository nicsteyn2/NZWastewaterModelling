
# This code is a minimal working example of running the particle filter conditional on a single value of θ.
# Results from this script are not used in the paper as we always marginalise out theta rather than conditioning on it.

using Plots

include("src/loadNZData.jl")
include("src/getOpts.jl")
include("src/covidParticleFilter.jl")
include("src/processResults.jl")
include("src/supportFuncs.jl")


function runFilter(θ, ST_DATE, EN_DATE, label; Nx=1e5, datasource="CW", alpha=3e9, h=30, windin=30)
    
    # Load data. Replace this with your own data loader to run on different data.
    Y = loadNZData(ST_DATE  = ST_DATE - Day(windin), EN_DATE = EN_DATE)
    
    # Fetch options dictionary and specify inputs
    opts = getOpts(Y)
    opts["Nx"] = Int(Nx)
    opts["datasource"] = datasource
    opts["alpha"] = alpha
    opts["resamplingWindow"] = h
    opts["windinPeriod"] = windin
    
    # Initialise hidden states and run the model
    (ρ, R, I) = initialiseHiddenStates(Y, opts)
    (ρ, R, I, STEPWISEEV) = covidParticleFilter(ρ, R, I, Y, θ, opts)
    
    # Process and save results
    results = processAllResults(ρ, R, I, θ, opts)
    println("Saving results to outputs/filter/ with label = " * label)
    CSV.write("outputs/filter/"*label*".csv", results)
    CSV.write("outputs/filter/"*label*"_rawdata.csv", Y)
    return(results, Y, opts)
    
end


# Specify parameter vector, start date, and end date and let the function do the rest!
θ = [0.01, 0.038, 150, 6.8e-6] # [σ_CAR, σ_R, k_c, k_w]
ST_DATE = Date("2023-01-01")
EN_DATE = Date("2023-10-01")
(results, Y, opts) = runFilter(θ, ST_DATE, EN_DATE, "testnewdata", Nx=1e6)


# We recommend plotting the results using R/makeFigure3_results.R.
# You can copy and paste the following parameters into the labelled section:
#
# # Set parameters
# filename = "example"
# filedir = "outputs/filter/"
# st_date = as.Date("2022-07-01")
# en_date = as.Date("2022-09-30")
# CARt_basetime = 50 # This is the time-point we compare CARt to to get our relative CARt
