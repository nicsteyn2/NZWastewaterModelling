
using DataFrames, CSV, Distributions, Statistics


# A NZ-specific funciton to load PMMH outputs
function fetchPMMHOutputsForNZ(mmhLabel::String, nChains, miniter)
    
    # Specify the standard directory where our PMMH results are saved
    indir = "outputs/pmmh.nosync/" * mmhLabel * "/"
    
    # Load PMMH samples for all periods
    all_samples = DataFrame()
    for period = 1:5
        current_label = mmhLabel * "_period" * string(period)
        current_samples = loadPMMHOutputs(current_label, nChains; indir=indir)
        current_samples.period .= period
        all_samples = [all_samples; current_samples]
    end
    
    # Remove missing data and cut-off the wind-in period
    all_samples = all_samples[all_samples[:,1] .!= 0, :]
    all_samples = all_samples[all_samples.iter .>= miniter, :]
    
    return(all_samples)
    
end


# And if we want to run each period on a different MMH output
function fetchPMMHOutputsForNZ(mmhLabel::Vector{String}, nChains, miniter)

    if length(mmhLabel) != 5
        error("Need to specify 5 mmh labels")
    end
        
    # Load PMMH samples for all periods
    all_samples = DataFrame()
    for period = 1:5
        current_mmh = mmhLabel[period]
        indir = "outputs/pmmh.nosync/" * current_mmh * "/"
        current_label = current_mmh * "_period" * string(period)
        current_samples = loadPMMHOutputs(current_label, nChains; indir=indir)
        current_samples.period .= period
        all_samples = [all_samples; current_samples]
    end
    
    # Remove missing data and cut-off the wind-in period
    all_samples = all_samples[all_samples[:,1] .!= 0, :]
    all_samples = all_samples[all_samples.iter .>= miniter, :]
    
    return(all_samples)
    
end




# Load data for all chains with a given label (see functions below to load with a specified period)
function loadPMMHOutputs(label::String, nChains::Int; indir="outputs/pmmh.nosync/")

    mcmc_samples = DataFrame()
    for chain_no = 1:nChains
        fname = indir * label * "_chainno" * string(chain_no) * "_retainedparams.csv"
        df_in = DataFrame(CSV.File(fname))
        df_in.iter = collect(1:nrow(df_in))
        df_in.chain_no .= chain_no
        mcmc_samples = [mcmc_samples; df_in]
    end
    return(mcmc_samples)

end



# Remove missing values, cut-off the wind-in, and ensure that all particle chains have the same length
function tidyPMMHOutputs(df_pmmh_in, min_iter)

    # Make a full copy to ensure we don't change anything mutable
    df_pmmh = deepcopy(df_pmmh_in)

    # Filter out all 0 outputs
    df_pmmh = df_pmmh[df_pmmh[:,1] .!= 0, :]

    # Filter out iterations prior to min_iter
    df_pmmh = df_pmmh[df_pmmh.iter .>= min_iter, :]

    # Ensure the max iter is the same for all chains
    max_iter_by_chain = combine(groupby(df_pmmh, :chain_no), :iter => maximum)
    min_max_iter = minimum(max_iter_by_chain.iter_maximum)
    df_pmmh = df_pmmh[df_pmmh.iter .<= min_max_iter,:]

    return(df_pmmh)

end


# Sample a set of parameters
function sampleParamsFromMMH(df_pmmh::DataFrame, n_samples::Int, paramNamesToSample)

    df_subset = df_pmmh[:, paramNamesToSample]
    sample_inds = sample(1:nrow(df_subset), n_samples)
    if n_samples == 1
        θ = vec(Matrix(df_subset[sample_inds,:]))
    else
        θ = Matrix(df_subset[sample_inds,:])
    end

    return(θ)

end

