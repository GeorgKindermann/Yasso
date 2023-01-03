using Downloads
eval(Meta.parse(String(take!(Downloads.download("https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/julia/y20_subroutine.jl", IOBuffer())))))
#include("y20_subroutine.jl") #In case you have it local on your computer


theta = [0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25]
init = zeros(5)
infall = [0.5, 0.1, 0.1, 0.2, 0] #0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
time = 1.   #Time to run
avgT = [-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.]  #Temp month average [C]
sumP = 600. #Precip annual summ [mm]
diam = 2.   #size [cm]
leach = 0.  #Leaching

yasso20.set_theta(theta)
yasso20.set_A(avgT, sumP, diam, leach)

for year in 0:9
    init = yasso20.get_next_timestep(init, infall, time)
    println("$year $(round.(init, digits=3))")
end

yasso20.get_spin(infall)


@time begin
    n = 1000000
    infall2 = copy(infall)
    result = zeros(5)
    yasso20.set_theta(theta)
    for i in 1:n
        infall2 .= infall * i / n
        avgT2 = avgT .+ i / n
        yasso20.set_A(avgT2, sumP, diam, leach)
        result .= yasso20.get_next_timestep(result, infall, time)
    end
end
# 13.324012 seconds (153.37 M allocations: 7.812 GiB, 3.84% gc time, 0.55% compilation time)

result'
# 1.40538  0.149233  0.223337  3.02026  6.67175
