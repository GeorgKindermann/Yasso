module yasso15

export set_theta, set_A, get_spin, get_next_timestep

theta = Vector{Float64}(undef, 35)
A = Matrix{Float64}(undef, 5, 5)

function set_theta(theta_)
  global theta = copy(theta_)
  tabs = [1:4; 32; 35]
  theta[tabs] .= -abs.(theta[tabs])
  nothing
end

function set_A(avgT, sumP, ampT, diam, leach)
  global A[1:4, 5] .= 0.0

  m3 = sumP / 1000.0
  #temperature annual cycle approximation
  m0 = (1.0 / sqrt(2.0) - 1.0) / pi
  m1 = 1.0 / sqrt(2.0) / pi
  m2 = (1.0 - 1.0 / sqrt(2.0)) / pi
  clim = 4.0 * ampT
  te = Vector{Float64}(undef, 4)
  te[1] = avgT + clim * m0
  te[2] = avgT - clim * m1
  te[3] = avgT + clim * m2
  te[4] = avgT + clim * m1

  #Average temperature dependence
  tem = sum(exp.(theta[22] * te + theta[23] * te .^ 2)) / 4
  temN = sum(exp.(theta[24] * te + theta[25] * te .^ 2)) / 4
  temH = sum(exp.(theta[26] * te + theta[27] * te .^ 2)) / 4

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

end
