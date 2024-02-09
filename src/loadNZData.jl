
using CSV, DataFrames, GLM, StatsModels, CategoricalArrays, Dates

# These functions are used to load the NZ wastewater and reported case data. They will need to be re-written for other jurisdictions.

function loadNZData(;ST_DATE=-Inf, EN_DATE=Inf)
            
    # Load and process wastewater data
    ww = DataFrame(CSV.File("data/wastewater.csv", missingstring="NA"))
    select!(ww, :date, :nGenomeCopies, :nPopInCatchment, :wwDataIsValid)
    
    # Process case data #TODO: DOCUMENT THE INT(ROUND()) TRANSFORMATION
    cases = DataFrame(CSV.File("data/cases.csv"))
    cases = separateAndDetrendCaseData(cases)
    cases.nCases = cases.nCasesDetrended # We want to just use our detrended cases as our cases
    cases.nCases = Int.(round.(cases.nCases)) # and cases need to be rounded as we use a discrete observation distribution
    cases[:,:casesDataIsValid] .= true
    select!(cases, :date, :nCases, :casesDataIsValid)

    # Start big table
    very_first_date = min(minimum(cases.date), minimum(ww.date))
    very_end_date = max(maximum(cases.date), maximum(ww.date))
    Y = DataFrame(date = very_first_date:Day(1):very_end_date)
    
    # Join tables together
    Y = leftjoin(Y, ww, on=:date)
    Y = leftjoin(Y, cases, on=:date)

    # Filter out data that is before our desired ST_DATE and after our desired EN_DATE
    if ST_DATE != -Inf
        Y = Y[Y.date .>= ST_DATE, :]
    end
    if EN_DATE != Inf
        Y = Y[Y.date .<= EN_DATE, :]
    end

    # Append integer time
    Y.t = Dates.value.(Y.date - minimum(Y.date)) .+ 1
    
    # Not all dates feature in each series, so we replace the missings with falses
    # On these days the model will just ignore the missing data
    Y.wwDataIsValid[ismissing.(Y.wwDataIsValid)] .= false
    Y.casesDataIsValid[ismissing.(Y.nCases)] .= false
    Y.nGenomeCopies[ismissing.(Y.nGenomeCopies)] .= NaN

    # Either remove days where nGenomeCopies == 0 or set it to a small value (since 0 is outside the support of our wastewater observation distribution)
    # Y.wwDataIsValid[Y.nGenomeCopies.==0] .= false
    tmp = Y.nGenomeCopies[.!isnan.(Y.nGenomeCopies)]
    tmp = tmp[tmp.!=0]
    Y.nGenomeCopies[Y.nGenomeCopies.==0] .= 0.1 * minimum(tmp)

    # Ensure the data is ordered by time and then return it
    sort!(Y, :t)
    return(Y)
    
end



function separateAndDetrendCaseData(Y)
    
    # We separate the NZ data into these dates
    st_dates = [Date("2000-01-01"), Date("2022-01-01"), Date("2022-03-01"), Date("2022-07-01"), Date("2022-10-01"), Date("2023-01-01"), Date("2023-04-01"), Date("2023-07-01"), Date("2023-10-01")]
    en_dates = [Date("2021-12-31"), Date("2022-02-28"), Date("2022-06-30"), Date("2022-09-30"), Date("2022-12-31"), Date("2023-03-31"), Date("2023-06-30"), Date("2023-09-30"), Date("2100-01-01")]
    
    # We break the input data into sections and then stack it into a big Yall for returning
    Yall = DataFrame()
    for ii in 1:length(st_dates)
        
        # Filter by dates and remove any missing data
        st = st_dates[ii]
        en = en_dates[ii]
        Yin = Y[(Y.date .>= st) .& (Y.date .<= en) .& .!ismissing.(Y.nCases),:]
        
        # If there is actually data to work with, detrend it!
        if nrow(Yin) > 0
            Yall = [Yall; detrendCaseData(Yin)]
        end
        
    end
    
    return(Yall)
    
end


function detrendCaseData(Y)
    
    # Process day of week
    Y[!, :DayOfWeek] = Dates.dayofweek.(Y[!, :date])
    
    daysofweek = Dict(1 => "Monday", 2 => "Tuesday", 3 => "Wednesday", 
    4 => "Thursday", 5 => "Friday", 6 => "Saturday", 
    7 => "Sunday")
    
    Y[!, :DayName] = map(d -> daysofweek[d], Y[!, :DayOfWeek])
    
    # Make sure DayName is treated as a categorical variable
    transform!(Y, :DayName => categorical => :DayName)
    
    # Run the linear regression model with log of nCases
    formula = @formula(log(nCases) ~ 1 + DayName)
    lm1 = lm(formula, Y)
    Y[!, :Prediction] = exp.(predict(lm1))
    
    # Detrend nCases multiplicatively
    Y[!, :nCasesDetrended] = Y[!, :nCases] ./ Y[!, :Prediction]
    
    # Adjust detrended data so that the total remains the same
    Y[!, :nCasesDetrended] = Y[!, :nCasesDetrended] .* sum(Y[!, :nCases]) / sum(Y[!, :nCasesDetrended])
    
    return(Y)
    
end