

using Plots, KernelDensity

include("loadNZData.jl")
include("getOpts.jl")
include("covidParticleFilter.jl")
include("processResults.jl")

function runCurvewiseEstimates(nTargetSamples::Int, θ::Array{Float64, 3}, Y::DataFrame, opts::Dict, label::String; outdir="outputs/curvewiseests/")
    
    
    # Calculate the number of samples we need to keep from each 
    nParamSamples = size(θ)[3]
    nFilterSamplesPerParam = Int(round(nTargetSamples/nParamSamples))
    nTotalSamples = nFilterSamplesPerParam * nParamSamples # This will equal nTargetSamples if nTargetSamples is a multiple of nParamSamples
    
    # Pre-allocate results vectors
    Rmax = zeros(nTotalSamples)
    Rmin = zeros(nTotalSamples)
    Rmax_ind = zeros(nTotalSamples)
    Rmin_ind = zeros(nTotalSamples)
    
    ρmax = zeros(nTotalSamples)
    ρmin = zeros(nTotalSamples)
    ρmax_ind = zeros(nTotalSamples)
    ρmin_ind = zeros(nTotalSamples)
    
    Imax = zeros(nTotalSamples)
    Imin = zeros(nTotalSamples)
    Imax_ind = zeros(nTotalSamples)
    Imin_ind = zeros(nTotalSamples)
    
    Rsub1_ind = zeros(nTotalSamples)

    # Get range of dates to consider
    t_ind_start = max(length(Y.date) - opts["resamplingWindow"], 30)
    t_ind_end = length(Y.date) - 30 # We use this value at the end to approximate how far we need to smooth
    dates = Y.date[t_ind_start:t_ind_end]

    progBar = Progress(nParamSamples; dt=1, barlen = 50, desc="Sampling from full posterior...")
    Threads.@threads for ii = 1:nParamSamples
        
        # We don't want to show progress for each indvidual filter
        opts["showProgress"] = false
        
        # Extract our current parameters
        θ_in = copy(θ[:,:,ii])

        # Solve the model
        (ρ, R, I) = initialiseHiddenStates(Y, opts)
        (ρ, R, I, _) = covidParticleFilter(ρ, R, I, Y, θ_in, opts)
        
        # Set the indices to sample from and save to
        sample_inds = sample(Int.(1:opts["Nx"]), nFilterSamplesPerParam)
        save_inds = (1:nFilterSamplesPerParam) .+ (nFilterSamplesPerParam * (ii - 1))
        
        # Save things
        (Rmax[save_inds], Rmax_ind[save_inds], Rmin[save_inds], Rmin_ind[save_inds]) = getExtremas(R[sample_inds, t_ind_start:t_ind_end])
        (ρmax[save_inds], ρmax_ind[save_inds], ρmin[save_inds], ρmin_ind[save_inds]) = getExtremas(ρ[sample_inds, t_ind_start:t_ind_end])
        (Imax[save_inds], Imax_ind[save_inds], Imin[save_inds], Imin_ind[save_inds]) = getExtremas(I[sample_inds, t_ind_start:t_ind_end])
        Rsub1_ind[save_inds] = getSubUnit(R[sample_inds, t_ind_start:t_ind_end])

        next!(progBar)
        
    end

    # Print to console
    println()
    println("Posterior estimates:")
    print("Max R = " * string(round(mean(Rmax), digits=2)) * " (" * string(round(quantile(Rmax, 0.025), digits=2)) * ", " * string(round(quantile(Rmax, 0.975), digits=2)) * ")")
    println(" on " * string(dates[Int(round(mean(Rmax_ind)))]) * " (" * string(dates[Int(quantile(Rmax_ind, 0.025))]) * ", " * string(dates[Int(quantile(Rmax_ind, 0.975))]) * ")")
    print("Min R = " * string(round(mean(Rmin), digits=2)) * " (" * string(round(quantile(Rmin, 0.025), digits=2)) * ", " * string(round(quantile(Rmin, 0.975), digits=2)) * ")")
    println(" on " * string(dates[Int(round(mean(Rmin_ind)))]) * " (" * string(dates[Int(quantile(Rmin_ind, 0.025))]) * ", " * string(dates[Int(quantile(Rmin_ind, 0.975))]) * ")")
    print("Max CAR = " * string(round(mean(ρmax), digits=2)) * " (" * string(round(quantile(ρmax, 0.025), digits=2)) * ", " * string(round(quantile(ρmax, 0.975), digits=2)) * ")")
    println(" on " * string(dates[Int(round(mean(ρmax_ind)))]) * " (" * string(dates[Int(quantile(ρmax_ind, 0.025))]) * ", " * string(dates[Int(quantile(ρmax_ind, 0.975))]) * ")")
    print("Min CAR = " * string(round(mean(ρmin), digits=2)) * " (" * string(round(quantile(ρmin, 0.025), digits=2)) * ", " * string(round(quantile(ρmin, 0.975), digits=2)) * ")")
    println(" on " * string(dates[Int(round(mean(ρmin_ind)))]) * " (" * string(dates[Int(quantile(ρmin_ind, 0.025))]) * ", " * string(dates[Int(quantile(ρmin_ind, 0.975))]) * ")")
    print("Max It = " * string(round(mean(Imax), digits=2)) * " (" * string(round(quantile(Imax, 0.025), digits=2)) * ", " * string(round(quantile(Imax, 0.975), digits=2)) * ")")
    println(" on " * string(dates[Int(round(mean(Imax_ind)))]) * " (" * string(dates[Int(quantile(Imax_ind, 0.025))]) * ", " * string(dates[Int(quantile(Imax_ind, 0.975))]) * ")")
    print("Min It = " * string(round(mean(Imin), digits=2)) * " (" * string(round(quantile(Imin, 0.025), digits=2)) * ", " * string(round(quantile(Imin, 0.975), digits=2)) * ")")
    println(" on " * string(dates[Int(round(mean(Imin_ind)))]) * " (" * string(dates[Int(quantile(Imin_ind, 0.025))]) * ", " * string(dates[Int(quantile(Imin_ind, 0.975))]) * ")")
    Rsub1_nonan = Rsub1_ind[.!isnan.(Rsub1_ind)]
    if length(Rsub1_nonan) >= 1
        println("R first sub 1 on " * string(dates[Int(round(mean(Rsub1_nonan)))]) * " (" * string(dates[Int(quantile(Rsub1_nonan, 0.025))]) * ", " * string(dates[Int(quantile(Rsub1_nonan, 0.975))]) * "). With " * string(sum(isnan.(Rsub1_ind))) * " NaN occurences.")
    else
        println("No R < 1 found. Non-nan Rsub1_inds are: " * string(Rsub1_nonan))
    end
    samples_out = DataFrame()

    # And print to file
    date_start = max(Y.date[length(Y.date) - opts["resamplingWindow"]], minimum(Y.date) + Day(30))
    date_end = Y.date[length(Y.date) - 30]

    if !isdir(outdir)
        mkdir(outdir)
    end
    open(outdir * label * "_curvewise.txt", "w") do file

        write(file, "Posterior estimates:\n")

        write(file, "Max R = " * string(round(mean(Rmax), digits=2)) * " (" * string(round(quantile(Rmax, 0.025), digits=2)) * ", " * string(round(quantile(Rmax, 0.975), digits=2)) * ")")
        write(file, " on " * string(dates[Int(round(mean(Rmax_ind)))]) * " (" * string(dates[Int(quantile(Rmax_ind, 0.025))]) * ", " * string(dates[Int(quantile(Rmax_ind, 0.975))]) * ")\n")

        write(file, "Min R = " * string(round(mean(Rmin), digits=2)) * " (" * string(round(quantile(Rmin, 0.025), digits=2)) * ", " * string(round(quantile(Rmin, 0.975), digits=2)) * ")")
        write(file, " on " * string(dates[Int(round(mean(Rmin_ind)))]) * " (" * string(dates[Int(quantile(Rmin_ind, 0.025))]) * ", " * string(dates[Int(quantile(Rmin_ind, 0.975))]) * ")\n")

        write(file, "Max CAR = " * string(round(mean(ρmax), digits=2)) * " (" * string(round(quantile(ρmax, 0.025), digits=2)) * ", " * string(round(quantile(ρmax, 0.975), digits=2)) * ")")
        write(file, " on " * string(dates[Int(round(mean(ρmax_ind)))]) * " (" * string(dates[Int(quantile(ρmax_ind, 0.025))]) * ", " * string(dates[Int(quantile(ρmax_ind, 0.975))]) * ")\n")

        write(file, "Min CAR = " * string(round(mean(ρmin), digits=2)) * " (" * string(round(quantile(ρmin, 0.025), digits=2)) * ", " * string(round(quantile(ρmin, 0.975), digits=2)) * ")")
        write(file, " on " * string(dates[Int(round(mean(ρmin_ind)))]) * " (" * string(dates[Int(quantile(ρmin_ind, 0.025))]) * ", " * string(dates[Int(quantile(ρmin_ind, 0.975))]) * ")\n")
        
        write(file, "Max It = " * string(round(mean(Imax), digits=2)) * " (" * string(round(quantile(Imax, 0.025), digits=2)) * ", " * string(round(quantile(Imax, 0.975), digits=2)) * ")")
        write(file, " on " * string(dates[Int(round(mean(Imax_ind)))]) * " (" * string(dates[Int(quantile(Imax_ind, 0.025))]) * ", " * string(dates[Int(quantile(Imax_ind, 0.975))]) * ")\n")

        write(file, "Min It = " * string(round(mean(Imin), digits=2)) * " (" * string(round(quantile(Imin, 0.025), digits=2)) * ", " * string(round(quantile(Imin, 0.975), digits=2)) * ")")
        write(file, " on " * string(dates[Int(round(mean(Imin_ind)))]) * " (" * string(dates[Int(quantile(Imin_ind, 0.025))]) * ", " * string(dates[Int(quantile(Imin_ind, 0.975))]) * ")\n")

        write(file, "R first sub 1 on " * string(dates[Int(round(mean(Rsub1_nonan)))]) * " (" * string(dates[Int(quantile(Rsub1_nonan, 0.025))]) * ", " * string(dates[Int(quantile(Rsub1_nonan, 0.975))]) * "). With " * string(sum(isnan.(Rsub1_ind))) * " NaN occurences.\n\n")

        write(file, "Note these estimates are constrained between " * string(date_start) * " and " * string(date_end) * ". Ensure the full confidence intervals of the relevant metric are easily contained in this range.\n")

    end

    # Finally we save the files to disk
    samples_out.Rmax = Rmax
    samples_out.Rmax_ind = Rmax_ind
    samples_out.Rmax_date = dates[Int.(Rmax_ind)]
    samples_out.Rmin = Rmin
    samples_out.Rmin_ind = Rmin_ind
    samples_out.Rmin_date = dates[Int.(Rmin_ind)]

    samples_out.CARmax = ρmax
    samples_out.CARmax_ind = ρmax_ind
    samples_out.CARmax_date = dates[Int.(ρmax_ind)]
    samples_out.CARmin = ρmin
    samples_out.CARmin_ind = ρmin_ind
    samples_out.CARmin_date = dates[Int.(ρmin_ind)]

    samples_out.Rsub1 = Rsub1_ind

    CSV.write(outdir * label * "_samples.csv", samples_out)

    
    return(samples_out)
    
end


function getExtremas(X)
    
    N = size(X)[1]
    
    max_values = zeros(N)
    min_values = zeros(N)
    max_indices = zeros(N)
    min_indices = zeros(N)
    for ii in 1:N
        max_values[ii], max_indices[ii] = findmax(X[ii, :])
        min_values[ii], min_indices[ii] = findmin(X[ii, :])
    end
    return(max_values, max_indices, min_values, min_indices)
    
end


function getSubUnit(X)

    N = size(X)[1]

    first_less_than_one_indices = zeros(N) .* NaN
    for ii in 1:N
        idx = findfirst(x -> x < 1, X[ii, :])
        if idx !== nothing
            first_less_than_one_indices[ii] = idx
        end
    end
    
    return(first_less_than_one_indices)

end