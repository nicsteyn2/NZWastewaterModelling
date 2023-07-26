
using Dates

# A collection of random functions that are useful here-and-there


# Maps a dataframe to a matrix, where the value at (i,j) is given by the value of :z_var that corresponds to x_var[i] and y_var[j]
function df_to_matrix(df, x_var, y_var, z_var)
    
    df_pivot = unstack(df, x_var, y_var, z_var)
    M = Matrix(df_pivot[!,2:end])
    return(M)
    
end


# Splits an array into equal parts of size n. The final part will have size <= n.
function split_into_parts(array, n)
    len = length(array)
    idx = 1
    result = []
    while idx <= len
        end_idx = min(idx + n - 1, len)
        push!(result, array[idx:end_idx])
        idx = end_idx + 1
    end
    return result
end


# Calculates the convolution of two discrete probability distributions
function conv(a, b)
    m = length(a)
    n = length(b)
    c = zeros(m+n-1)
    for i = 1:m
        for j = 1:n
            c[i+j-1] += a[i]*b[j]
        end
    end
    return c
end


# Plots a 2D heatmap of estimated logliks (or other things)
function plot2DHeatmap(data, var1, var2, zvar)
    
    v1 = sort(unique(data[:,var1]))
    v2 = sort(unique(data[:,var2]))
    z = df_to_matrix(data, var1, var2, zvar)
    p_loglik = heatmap(v1, v2, z', right_margin=4Plots.mm, dpi=300)
    return(p_loglik)
    
end


# This expands an NxP matrix of parameters where each row corresponds to the parametesr between st_date[i] and en_date[i] into a full TxP matrix of parameters
function makeParamMatrix(θmat::Matrix, all_dates::Vector{Date}, st_date::Vector{Date}, en_date::Vector{Date}; param_names=["sigRho", "sigR", "kc", "kw"])
    
    num_periods = size(θmat)[1]
    num_params = size(θmat)[2]
    
    # Check all the sizes match
    if length(st_date) != num_periods
        @error("Length of st_date vector doesn't match number of parameter sets.")
    elseif length(en_date) != num_periods
        @error("Length of en_date vector doesn't match number of parameter sets.")
    elseif length(param_names) != num_params
        @error("Length of param_names vector doesn't match number of parameters.")
    end
    
    # This function finds the index of all dates between st_date and en_date
    function getParamInds(date_vec, st_date, en_date)
        update_ind = (date_vec .>= st_date) .* (date_vec .<= en_date)
        return(update_ind)
    end
    
    # Construct parameter df
    df_params = DataFrame(period=collect(1:num_periods))
    for ii = 1:num_params
        df_params[:,param_names[ii]] = θmat[:,ii]
    end
    
    # Construct dates df
    df_dates = DataFrame(date=all_dates, period=zeros(length(all_dates)))
    for jj = 1:num_periods
        df_dates.period[getParamInds(df_dates.date, st_date[jj], en_date[jj])] .= jj
    end

    # Check every date has been assigned a period
    if sum(df_dates.period .== 0) > 0
        @error("At least one date in all_dates is not included in the defined periods.")
    end

    # Create final table
    df = leftjoin(df_dates, df_params, on=:period)
    θ_out = Matrix(df[:,3:end])
    return(θ_out)

end


# Create .json of mmhOpts for saving
function clean_opts_dict(D; VERBOSE=false)

    D_clean = Dict()

    for (k, v) in D
        try
            JSON.json(v)
            D_clean[k] = v
        catch e
            if VERBOSE
                println("Cannot serialize key: ", k)
            end
        end
    end

    return(D_clean)

end
    