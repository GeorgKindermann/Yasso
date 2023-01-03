module yasso20

export set_theta, set_A, get_spin, get_next_timestep

theta = Vector{Float64}(undef, 35)
A = Matrix{Float64}(undef, 5, 5)

function set_theta(theta_)
  global theta = copy(theta_)
  tabs = [1:4; 32; 35]
  theta[tabs] .= -abs.(theta[tabs])
  nothing
end

function set_A(avgT, sumP, diam, leach)
  global A[1:4, 5] .= 0.0

  m3 = sumP / 1000.0

  #Average temperature dependence
  tem = sum(exp.(theta[22] * avgT + theta[23] * avgT .^ 2)) / 12.
  temN = sum(exp.(theta[24] * avgT + theta[25] * avgT .^ 2)) / 12.
  temH = sum(exp.(theta[26] * avgT + theta[27] * avgT .^ 2)) / 12.

  #Precipitation dependence
  tem = tem * (1.0 - exp(theta[28] * m3))
  temN = temN * (1.0 - exp(theta[29] * m3))
  temH = temH * (1.0 - exp(theta[30] * m3))

  #Size class dependence -- no effect if d == 0.0
  size_dep = 1
  if (diam > 0.0)
    size_dep = min(1.0, (1.0 + theta[33] * diam + theta[34] * diam^2)^theta[35])
  end
  #Calculating matrix a (will work ok despite the sign of alphas)
  A[1:6:13] = theta[1:3] * tem * size_dep
  A[4, 4] = theta[4] * temN * size_dep
  A[5, 5] = theta[32] * temH #no size effect in humus
  dAbs = abs.(A[1:6:19])
  idx = 5
  for i in 0:3
    for j in 0:3
      if (i != j)
        A[1+j*5+i] = theta[idx] * dAbs[1+j]
        idx += 1
      end
    end
  end
  #mass flows AWEN -> H (size effect is present here)
  A[5:5:20] .= theta[31] * dAbs
  #Leaching (no leaching for humus) 
  if (leach < 0.0)
    A[1:6:19] .+= leach * m3
  end
  nothing
end

function get_spin(infall)
  A \ infall * -1
end

function get_next_timestep(init, infall, time)
  A \ (exp(A * time) * (A * init + infall) - infall)
end

#import Expokit #Does not improve mutch
#function get_next_timestepB(init, infall, time)
#  A \ (Expokit.chbv(A * time, A * init + infall) - infall)
#end

end
