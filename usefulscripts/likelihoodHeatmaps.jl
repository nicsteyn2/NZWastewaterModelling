
using Plots, IterTools, Dates

include("../src/estimateLikelihoods.jl")
include("../src/loadData.jl")
include("../src/getOpts.jl")
include("../src/covidParticleFilter.jl")


function estimate1D(θbase, paramVals, paramInd, ST_DATE, EN_DATE; extra_windin=50, Nx=1e5)

    nCores = Threads.nthreads()
    ST_DATE = ST_DATE - Day(30) - Day(extra_windin)

    # Load data and setup things
    Y = loadData(ST_DATE=ST_DATE, EN_DATE=EN_DATE)
    opts = getOpts(Y)
    opts["Nx"] = Int(Nx)
    opts["showProgress"] = false
    θ = deepcopy(θbase)

    # Separate into epochs for parfor
    inds = sample(collect(1:length(paramVals)), length(paramVals); replace=false)
    inds_grouped = split_into_parts(inds, max(nCores, 16))

    # Pre-allocate
    results = repeat([NaN], length(paramVals))
    
    # Run results and save
    for current_inds_ind = 1:length(inds_grouped)
        println("Running parameter set " * string(current_inds_ind) * " of " * string(length(inds_grouped)))
        current_inds = inds_grouped[current_inds_ind]
        test = estimateLikelihoodAtPoints(paramVals[current_inds], paramInd, extra_windin, Y, θ, opts)
        results[current_inds] = test
    end
    
    return(results)
    
end


function estimate2DHeatmap(θbase, x1_vals, x2_vals, x1_ind, x2_ind, ST_DATE, EN_DATE, label; extra_windin=50, Nx=1e5)

    nCores = Threads.nthreads()
    ST_DATE = ST_DATE - Day(30) - Day(extra_windin)
    
    # Load data and setup things
    Y = loadData(ST_DATE=ST_DATE, EN_DATE=EN_DATE)
    opts = getOpts(Y)
    opts["Nx"] = Int(Nx)
    opts["showProgress"] = false
    θ = deepcopy(θbase)
    
    # Get all combinations of input parameters
    combinations = vec(collect(product(x1_vals, x2_vals)))
    paramVals = hcat([c[1] for c in combinations], [c[2] for c in combinations])
    
    # Separate into epochs for parfor
    inds = sample(collect(1:size(paramVals)[1]), size(paramVals)[1]; replace=false)
    inds_grouped = split_into_parts(inds, max(nCores, 16))
    
    # Results
    results = DataFrame(paramVals, ["x1", "x2"])
    results.loglik = NaN * zeros(nrow(results))
    
    # Run results and save
    for current_inds_ind = 1:length(inds_grouped)
        println("Running parameter set " * string(current_inds_ind) * " of " * string(length(inds_grouped)))
        current_inds = inds_grouped[current_inds_ind]
        results.loglik[current_inds] = estimateLikelihoodAtPoints(paramVals[current_inds,:], [x1_ind, x2_ind], extra_windin, Y, θ, opts)
        CSV.write("outputs/misc/2dheatmap_" * label * ".csv", results)
    end
    
    return(results)

end



# ----- Example 1D estimates ------

# Specify inputs here
θ = [0.02, 0.04, 300, 5.5e-6] # Need to specify the parameters that we don't vary
sigR_vals = 0.014:0.004:0.1 # Vector of first parameter to vary
ST_DATE = Date("2022-07-01") # Start of time-period (ignoring wind-in)
EN_DATE = Date("2022-09-30") # End of time-period

# Estimate the likelihoods
estimated_logliks = estimate1D(θ, sigR_vals, 2, ST_DATE, EN_DATE; Nx=3e5)

# Plot the results
plt = scatter(sigR_vals, estimated_logliks, label=false)
plt = xlabel!("σ_R")
plt = ylabel!("Log-likelihood")
png(plt, "outputs/results/supplementary/suppl_heatmap_1d.png")



# ----- Example 2D estimates ------

# Specify inputs here
θ = [0.02, 0.04, 300, 5.5e-6] # Need to specify the parameters that we don't vary
sigRho_vals = 0.005:0.005:0.05 # Vector of first parameter to vary
kc_vals = 10:200:2010 # Vector of second parameter to vary
ST_DATE = Date("2022-07-01") # Start of time-period (ignoring wind-in)
EN_DATE = Date("2022-09-30") # End of time-period

# Estimate the likelihoods
estimated_logliks = estimate2DHeatmap(θ, sigRho_vals, kc_vals, 1, 3, ST_DATE, EN_DATE, "sigrho_kc_period3")

# Plot the results
data_filt = estimated_logliks[estimated_logliks.loglik.>=-1900,:]
data_filt.lik = exp.(data_filt.loglik .- maximum(data_filt.loglik))
plt = plot2DHeatmap(data_filt, :x1, :x2, :loglik)
plt = xlabel!("σ_CAR")
plt = ylabel!("k_c")
png(plt, "outputs/results/supplementary/suppl_heatmap_2d.png")