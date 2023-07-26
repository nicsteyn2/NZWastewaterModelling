
include("covidParticleFilter.jl")

# This file contains a bunch of functions that are useful for estimating the likelihood of models.


# Estimate the log-likelihood at a fixed parameter value. extra_windin is used to reduce the variance of the LOGEV estimates.
# This function is called by PMMH.jl so edit with caution!
function estimateLikelihoodAtFixedParam(θ, extra_windin, Y, opts)
    
    (ρ, R, I) = initialiseHiddenStates(Y, opts)
    (_, _, _, STEPWISEEV) = covidParticleFilter(ρ, R, I, Y, θ, opts)
    LOGEV = sum(log.(STEPWISEEV[(opts["windinPeriod"]+extra_windin+1):end]))
    
    return(LOGEV)
    
end


# Wrapper function for estimateLikelihoodAtFixedParam()
# Sets the parameter value of θbase at position param_inds[ii] to param_vals[ii] and then estimates the likelihood
function estimateLikelihood(param_vals, param_inds, extra_windin, Y, θbase, opts)
    
    θ = deepcopy(θbase)
    
    # Update parameter vector
    for ii in 1:length(param_vals)
        θ[param_inds[ii]] = param_vals[ii]
    end
    
    LOGEV = estimateLikelihoodAtFixedParam(θ, extra_windin, Y, opts)    
    return(LOGEV)
    
end


# When we want to estimate the likelihood at various values of θ, where only one element is changed between estimates
function estimateLikelihoodAtPoints(param_vals::Vector, param_ind, extra_windin, Y, θ, opts)
    
    loglik = nothing
    
    n_pts = length(param_vals)
    loglik = zeros(n_pts)
    progBar = Progress(n_pts; dt=1, barlen = 50, desc="Estimating likelihoods...")
    Threads.@threads for ii = 1:n_pts
        loglik[ii] = estimateLikelihood(param_vals[ii], param_ind, extra_windin, Y, θ, opts)
        next!(progBar)
    end
    
    return(loglik)
    
end


# When we want to estimate the likelihood at various values of θ, where multiple elements are changed between estimates
function estimateLikelihoodAtPoints(param_vals::Matrix, param_inds, extra_windin, Y, θ, opts)
    
    loglik = nothing
    
    n_pts = size(param_vals)[1]
    loglik = zeros(n_pts)
    progBar = Progress(n_pts; dt=1, barlen = 50, desc="Estimating likelihoods...")
    Threads.@threads for ii = 1:n_pts
        loglik[ii] = estimateLikelihood(param_vals[ii,:], param_inds, extra_windin, Y, θ, opts)
        next!(progBar)
    end
    
    return(loglik)
    
end


# This function calls estimateLikelihoodAtFixedParam() many times to estimate the variance of our likelihood estimates
# Note the theoretical target standard deviation of log-likelihood estimates for PMMH purposes is 1.2-1.3 (DOI: 10.1214/14-STS511)
function estimateVarianceOfLogLikEsts(Ntrials, θ, extra_windin, Y, opts)
    
    logevidences = zeros(Ntrials)
    opts["showProgress"] = false # Disable individual particle filter progress bars
    
    # Run model
    progBar = Progress(Ntrials; dt=1, barlen = 50, desc="Estimating likelihoods...")
    Threads.@threads for ii = 1:Ntrials
        logevidences[ii] = estimateLikelihoodAtFixedParam(θ, extra_windin, Y, opts)
        next!(progBar)
    end
    
    return(var(logevidences), logevidences)
    
end