
using DataFrames

include("../src/loadNZData.jl")
include("../src/getOpts.jl")


function fetchPosteriorParameterMeans(df_pmmh::DataFrame, periods::Int, paramNames)

    out = zeros(length(periods), length(paramNames))
    for ii = 1:length(periods)
        period_matr = Matrix(df_pmmh[df_pmmh.period.==periods[ii], paramNames])
        out[ii,:] = mean(period_matr, dims=1)
    end
    return(out)

end


# Save distributions
function exportDistributions()
    Y = loadNZData() # Only need this so we can call getOpts()
    opts = getOpts(Y)
    CSV.write("outputs/results/misc/infToRepDist.csv", DataFrame(t=0:(opts["maxInfToRepTime"]-1), p=opts["infToRepDist"]))
    CSV.write("outputs/results/misc/infToShedDist.csv", DataFrame(t=0:(opts["maxSheddingTime"]-1), p=opts["infToShedDist"]))
    CSV.write("outputs/results/misc/genTimeDist.csv", DataFrame(t=0:(length(opts["genTimeDist"])-1), p=opts["genTimeDist"]))
end