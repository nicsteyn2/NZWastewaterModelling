
using CSV, DataFrames, Plots, Dates

# This script requires the outpouts of resultsFullPosterior.jl to already exist

label = "final"

core_results = CSV.read("outputs/fullposterior/" * label * ".csv", DataFrame)
raw_data = CSV.read("outputs/fullposterior/" * label * "_rawdata.csv", DataFrame)


# Results paragraph one
rt_jan1 = core_results[(core_results.date.==Date("2022-01-01")) .& (core_results.variable.=="Rt"),:]

println()
println("------------------------------------------")
println("Rt on 1 Jan 2022 = " * string(round(rt_jan1.mean[1], digits=2)) * " (" * string(round(rt_jan1.lower[1], digits=2)) * ", " * string(round(rt_jan1.upper[1], digits=2)) * ")")
println("------------------------------------------")
println()


# Results paragraph two
raw_data_wave_one = raw_data[raw_data.date .<= Date("2022-05-01"), :]
max_cases_wave_one = maximum(raw_data_wave_one.nCases)

raw_data_wave_two = raw_data[(raw_data.date .>= Date("2022-05-01")) .& (raw_data.date .<= Date("2022-09-01")), :]
max_cases_wave_two = maximum(raw_data_wave_two.nCases)

infections = core_results[core_results.variable.=="It",:]
infections_w1 = infections[infections.date .<= Date("2022-05-01"),:]
infections_w2 = infections[(infections.date .>= Date("2022-05-01")).&(infections.date .<= Date("2022-09-01")), :]
max_infs_w1 = maximum(infections_w1.mean)
max_infs_w2 = maximum(infections_w2.mean)

println()
println("------------------------------------------")
println("Wave one cases peak = " * string(max_cases_wave_one) * " vs wave two cases peak = " * string(max_cases_wave_two))
println("Second peak was " * string((max_cases_wave_two/max_cases_wave_one)*100) * "% of the first")
println()
println("Wave one infections peak = " * string(max_infs_w1) * " vs wave two infs peak = " * string(max_infs_w2))
println("Second peak infs were " * string((max_infs_w2/max_infs_w1)*100) * "% of the first")
println("------------------------------------------")
println()


# CARt on 1 April 2022
results_2e9 = CSV.read("outputs/fullposterior/" * label * "_alpha2e9.csv", DataFrame)
results_3e9 = CSV.read("outputs/fullposterior/" * label * ".csv", DataFrame)
results_4e9 = CSV.read("outputs/fullposterior/" * label * "_alpha4e9.csv", DataFrame)

car_2e9 = results_2e9.mean[(results_2e9.date .== Date("2022-04-01")) .& (results_2e9.variable.=="CARt"),:][1]
car_2e9_l = results_2e9.lower[(results_2e9.date .== Date("2022-04-01")) .& (results_2e9.variable.=="CARt"),:][1]
car_2e9_u = results_2e9.upper[(results_2e9.date .== Date("2022-04-01")) .& (results_2e9.variable.=="CARt"),:][1]

car_3e9 = results_3e9.mean[(results_3e9.date .== Date("2022-04-01")) .& (results_3e9.variable.=="CARt"),:][1]
car_3e9_l = results_3e9.lower[(results_3e9.date .== Date("2022-04-01")) .& (results_3e9.variable.=="CARt"),:][1]
car_3e9_u = results_3e9.upper[(results_3e9.date .== Date("2022-04-01")) .& (results_3e9.variable.=="CARt"),:][1]

car_4e9 = results_4e9.mean[(results_4e9.date .== Date("2022-04-01")) .& (results_4e9.variable.=="CARt"),:][1]
car_4e9_l = results_4e9.lower[(results_4e9.date .== Date("2022-04-01")) .& (results_4e9.variable.=="CARt"),:][1]
car_4e9_u = results_4e9.upper[(results_4e9.date .== Date("2022-04-01")) .& (results_4e9.variable.=="CARt"),:][1]

println()
println("-------------------CAR ON 1 APRIL 2022-----------------------")
println("α = 2e9 implies CAR = " * string(round(car_2e9, digits=2)) * " (" * string(round(car_2e9_l, digits=2)) * ", " * string(round(car_2e9_u, digits=2)) * ")")
println("α = 3e9 implies CAR = " * string(round(car_3e9, digits=2)) * " (" * string(round(car_3e9_l, digits=2)) * ", " * string(round(car_3e9_u, digits=2)) * ")")
println("α = 4e9 implies CAR = " * string(round(car_4e9, digits=2)) * " (" * string(round(car_4e9_l, digits=2)) * ", " * string(round(car_4e9_u, digits=2)) * ")")
println("-------------------------------------------------------------")
println()




# CARt on 1 Jan 2022 (for older comparisons)
results_2e9 = CSV.read("outputs/fullposterior/" * label * "_alpha2e9.csv", DataFrame)
results_3e9 = CSV.read("outputs/fullposterior/" * label * ".csv", DataFrame)
results_4e9 = CSV.read("outputs/fullposterior/" * label * "_alpha4e9.csv", DataFrame)

car_2e9 = results_2e9.mean[(results_2e9.date .== Date("2022-01-01")) .& (results_2e9.variable.=="CARt"),:][1]
car_2e9_l = results_2e9.lower[(results_2e9.date .== Date("2022-01-01")) .& (results_2e9.variable.=="CARt"),:][1]
car_2e9_u = results_2e9.upper[(results_2e9.date .== Date("2022-01-01")) .& (results_2e9.variable.=="CARt"),:][1]

car_3e9 = results_3e9.mean[(results_3e9.date .== Date("2022-01-01")) .& (results_3e9.variable.=="CARt"),:][1]
car_3e9_l = results_3e9.lower[(results_3e9.date .== Date("2022-01-01")) .& (results_3e9.variable.=="CARt"),:][1]
car_3e9_u = results_3e9.upper[(results_3e9.date .== Date("2022-01-01")) .& (results_3e9.variable.=="CARt"),:][1]

car_4e9 = results_4e9.mean[(results_4e9.date .== Date("2022-01-01")) .& (results_4e9.variable.=="CARt"),:][1]
car_4e9_l = results_4e9.lower[(results_4e9.date .== Date("2022-01-01")) .& (results_4e9.variable.=="CARt"),:][1]
car_4e9_u = results_4e9.upper[(results_4e9.date .== Date("2022-01-01")) .& (results_4e9.variable.=="CARt"),:][1]

println()
println("----------------- CAR ON 1 JAN 2022 -------------------------")
println("α = 2e9 implies CAR = " * string(round(car_2e9, digits=2)) * " (" * string(round(car_2e9_l, digits=2)) * ", " * string(round(car_2e9_u, digits=2)) * ")")
println("α = 3e9 implies CAR = " * string(round(car_3e9, digits=2)) * " (" * string(round(car_3e9_l, digits=2)) * ", " * string(round(car_3e9_u, digits=2)) * ")")
println("α = 4e9 implies CAR = " * string(round(car_4e9, digits=2)) * " (" * string(round(car_4e9_l, digits=2)) * ", " * string(round(car_4e9_u, digits=2)) * ")")
println("-------------------------------------------------------------")
println()

