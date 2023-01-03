using Downloads
<<<<<<< HEAD:15/julia/simpleExample.jl
eval(Meta.parse(String(take!(Downloads.download("https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/15/julia/y15_subroutine.jl", IOBuffer())))))
=======
eval(Meta.parse(String(take!(Downloads.download("https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/15/y15_subroutine.jl", IOBuffer())))))
>>>>>>> 140064075045444e2d0fe1d14ac8a9b3060cbc11:15/simpleExample.jl
#include("y15_subroutine.jl") #In case you have it local on your computer


theta = [0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26]
init = zeros(5)
infall = [0.5, 0.1, 0.1, 0.2, 0] #0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
time = 1.   #Time to run
avgT = 10.  #Temp annual average [C]
sumP = 600. #Precip annual summ [mm]
ampT = 12.  #Amplitude (max. difference of month averages / 2) [C]
diam = 2.   #size [cm]
leach = 0.  #Leaching

yasso15.set_theta(theta)
yasso15.set_A(avgT, sumP, ampT, diam, leach)

for year in 0:9
    init = yasso15.get_next_timestep(init, infall, time)
    println("$year $(round.(init, digits=3))")
end

yasso15.get_spin(infall)


@time begin
    n = 1000000
    infall2 = copy(infall)
    result = zeros(5)
    yasso15.set_theta(theta)
    for i in 1:n
        infall2 .= infall * i / n
        avgT2 = avgT + i / n
        yasso15.set_A(avgT2, sumP, ampT, diam, leach)
        result .= yasso15.get_next_timestep(result, infall, time)
    end
end
#12.989973 seconds (160.00 M allocations: 8.702 GiB, 4.68% gc time)
result
# 2.61325  0.275456  0.392325  8.35314  10.5658

