
using Distributions

include("covidParticleFilter.jl")
include("estimateLikelihoods.jl")

# This provides the primary functions for running particle marginal Metropolis-Hastings
# These methods should all be fairly generic.


# Run multiple threads of PMMH
function runParticleMMHMultiThread(nThreads, Y, mmhOpts, opts, label; extra_windin=50, outdir="/outputs/pmmh.nosync/")
    
    Threads.@threads for ii = 1:nThreads
        runParticleMMH(Y, mmhOpts, opts, label, chainNo=ii, extra_windin=extra_windin, outdir=outdir)
    end
    
end

# Run one thread of PMMH
function runParticleMMH(Y::DataFrame, mmhOpts::Dict, opts::Dict, label::String; chainNo=-1, extra_windin=50, outdir="outputs/pmmh.nosync/")
    
    # Ensure the individual progress bars are disabled
    opts["showProgress"] = false
    
    # Check that the output directory exists (and create it if it doesn't)
    if !isdir(outdir)
        mkdir(outdir)
    end
    
    # Pre-allocate results tables
    allParametersAndLikelihoods = zeros(mmhOpts["nSteps"], mmhOpts["nParams"]+2) #final two columns are did_accept, log_likelihood
    logevidence = zeros(mmhOpts["nSteps"])
    retainedParameters = zeros(mmhOpts["nSteps"], mmhOpts["nParams"])
    acceptanceProbs = repeat([NaN], mmhOpts["nSteps"])
    didAccept = repeat([NaN], mmhOpts["nSteps"])
    
    # When this counter hits 5 we re-sample the original likelihood so we don't get stuck on some unreasonably likely value
    rejectCounter = 0 
    
    # Sample initial parameters and logevidence
    retainedParameters[1,:] = sampleInitialParameters(mmhOpts)
    logevidence[1] = estimateLikelihoodAtFixedParam(retainedParameters[1,:], extra_windin, Y, opts)
    
    # Save initial parameters and logevidence
    allParametersAndLikelihoods[1,1:mmhOpts["nParams"]] = retainedParameters[1,:]
    allParametersAndLikelihoods[1,(mmhOpts["nParams"]+1):end] = [NaN, logevidence[1]]
    
    # Run the chain
    progBar = Progress(mmhOpts["nSteps"]-1; dt=1, barlen = 50, desc="Running particle MMH (progress for thread 1 only)...", enabled=(Threads.threadid()==1))
    for ii = 2:mmhOpts["nSteps"]
        
        # Sample new parameter vector, estimate log-likelihood, calcualte acceptance prob, and choose to accept or reject
        θ_test = sampleParamProposal(retainedParameters[ii-1,:], mmhOpts)
        logev_test = estimateLikelihoodAtFixedParam(θ_test, extra_windin, Y, opts)
        acceptance_prob = calculateAcceptanceProbability(logev_test, logevidence[ii-1], θ_test, retainedParameters[ii-1,:], mmhOpts)
        did_accept = rand() < acceptance_prob
        
        # Save results
        allParametersAndLikelihoods[ii,1:mmhOpts["nParams"]] = θ_test
        allParametersAndLikelihoods[ii,(mmhOpts["nParams"]+1):end] = [did_accept, logev_test]
        acceptanceProbs[ii] = acceptance_prob
        didAccept[ii] = did_accept
        
        if did_accept # Keep our new parameter and reset the reject counter
            retainedParameters[ii,:] = θ_test
            logevidence[ii] = logev_test
            rejectCounter = 0
        elseif rejectCounter <= mmhOpts["resampleLogLikFreq"] # otherwise we reject, keep the old parameter, and increment the reject counter
            retainedParameters[ii,:] = retainedParameters[ii-1,:]
            logevidence[ii] = logevidence[ii-1]
            rejectCounter = rejectCounter + 1
        else # and if this reject counter > mmhOpts["resampleLogLikFreq"], we also resample the old likelihood
            println("Resampling current likelihood...")
            retainedParameters[ii,:] = retainedParameters[ii-1,:]
            logevidence[ii] = estimateLikelihoodAtFixedParam(retainedParameters[ii,:], extra_windin, Y, opts)
            rejectCounter = 0
        end
        
        # Write results to disk after specified number of iterations
        if (ii % mmhOpts["mmhSave"]) == 0
            println("ii = " * string(ii) * ", saving outputs to disk")
            CSV.write(outdir * label * "_chainno" * string(chainNo) * "_allparamlogliks.csv", DataFrame(allParametersAndLikelihoods, mmhOpts["paramNamesWithExtra"]))
            CSV.write(outdir * label * "_chainno" * string(chainNo) * "_retainedparams.csv", DataFrame(retainedParameters, mmhOpts["paramNames"]))
            CSV.write(outdir * label * "_chainno" * string(chainNo) * "_logliks.csv", DataFrame(logev=logevidence), writeheader=false)
            CSV.write(outdir * label * "_chainno" * string(chainNo) * "_acceptanceprobs.csv", DataFrame(acceptprob=acceptanceProbs), writeheader=false)
        end
        
        next!(progBar)
        
    end
    
    CSV.write(outdir * label * "_chainno" * string(chainNo) * "_allparamlogliks.csv", DataFrame(allParametersAndLikelihoods, mmhOpts["paramNamesWithExtra"]))
    CSV.write(outdir * label * "_chainno" * string(chainNo) * "_retainedparams.csv", DataFrame(retainedParameters, mmhOpts["paramNames"]))
    CSV.write(outdir * label * "_chainno" * string(chainNo) * "_logliks.csv", DataFrame(logev=logevidence), writeheader=false)
    CSV.write(outdir * label * "_chainno" * string(chainNo) * "_acceptanceprobs.csv", DataFrame(acceptprob=acceptanceProbs), writeheader=false)
    
end


# This calculates the standard Metropolis-Hastings acceptance probability.
function calculateAcceptanceProbability(logev_new, logev_old, θ_new, θ_old, mmhOpts)
    
    evidence_ratio  = exp(logev_new - logev_old)
    
    prior_ratio_num = 1
    prior_ratio_den = 1
    prop_ratio_num = 1
    prop_ratio_den = 1
    for ii = 1:mmhOpts["nParams"]
        if mmhOpts["estimated"][ii]
            prior_ratio_num = prior_ratio_num * pdf(mmhOpts["paramPriorDists"][ii], θ_new[ii])
            prior_ratio_den = prior_ratio_den * pdf(mmhOpts["paramPriorDists"][ii], θ_old[ii])
            prop_ratio_num = prop_ratio_num * pdf(mmhOpts["paramProposalDist"][ii](θ_new[ii], mmhOpts), θ_old[ii])
            prop_ratio_den = prop_ratio_den * pdf(mmhOpts["paramProposalDist"][ii](θ_old[ii], mmhOpts), θ_new[ii])
        end
    end
    prior_ratio = prior_ratio_num/prior_ratio_den
    prop_ratio = prop_ratio_num/prop_ratio_den
    
    acceptanceProb = evidence_ratio * prior_ratio * prop_ratio    
    return(min(acceptanceProb, 1))
    
end


# Sample from the parameter proposal distributions
function sampleParamProposal(θ_old, mmhOpts)
    
    θ_new = zeros(mmhOpts["nParams"])
    for ii = 1:mmhOpts["nParams"]
        
        if mmhOpts["estimated"][ii]
            θ_new[ii] = rand(mmhOpts["paramProposalDist"][ii](θ_old[ii], mmhOpts))
        else
            θ_new[ii] = mmhOpts["defaultValues"][ii]
        end
        
    end
    
    return(θ_new)
    
end


# Sample initial values of parameters
function sampleInitialParameters(mmhopts)
    
    paramVector = zeros(mmhopts["nParams"])
    for ii = 1:mmhopts["nParams"]
        # If we don't specify intial bounds just sample from the prior
        if ismissing(mmhopts["paramInitBounds"][ii])
            paramVector[ii] = rand(mmhopts["paramPriorDists"][ii])
            # Otherwise we sample from some intiial bounds. This allows us to run the method with a much shorter wind-in period.
        else
            paramVector[ii] = rand(Uniform(mmhopts["paramInitBounds"][ii][1], mmhopts["paramInitBounds"][ii][2]))
        end
    end
    
    return(paramVector)
    
end