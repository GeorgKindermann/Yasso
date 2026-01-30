using Downloads
eval(Meta.parse(String(take!(Downloads.download("https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/julia/y20_subroutine.jl", IOBuffer())))))
#include("y20_subroutine.jl") #In case you have it local on your computer

using StaticArrays

const theta = yasso20.get_theta()

init = zero(SVector{5, Float64})
infall = SVector(0.5, 0.1, 0.1, 0.2, 0.0)  #0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
time = 1.   #Time to run
avgT = SVector(-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.)   #Temp month average [deg C]
sumP = 600. #Precip annual summ [mm]
diam = 2.   #size [cm]
leach = 0.  #Leaching

A = yasso20.get_A(theta, avgT, sumP, diam, leach) 

for year in 0:9
    init = yasso20.get_next_timestep(A, init, infall, time)
    println("$year $(round.(init, digits=3))")
end

yasso20.get_spin(A, infall)

function tt(n)
    result = zero(SVector{5, Float64})
    for i in 1:n
        infall2 = infall * i / n
        avgT2 = avgT .+ i / n
        A = yasso20.get_A(theta, avgT2, sumP, diam, leach)
        result = yasso20.get_next_timestep(A, result, infall, time, 6)
    end
    result
end

tt(1)
@time result = tt(1000000)
# 0.600205 seconds (7.00 M allocations: 350.945 MiB, 3.32% gc time)

result'
# 1.40538  0.149233  0.223337  3.02026  6.67175
