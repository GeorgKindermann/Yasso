module yasso20

export get_theta, get_A, get_spin, get_next_timestep

using LinearAlgebra, StaticArrays

@inline @fastmath function exp5FIX(A::SMatrix{5, 5, Float64, 25})
    s = 0.0
    @inbounds for i in 1:25; s += A[i]^2; end
    p = sqrt(s)
    
    normiter = 2.0; j = 2
    while p > normiter; normiter *= 2.0; j += 1; end

    C = A / normiter
    res = one(SMatrix{5,5,Float64,25}) + C
    
    D = C * C; res += D * 0.5
    D = D * C; res += D * 0.16666666666666666
    D = D * C; res += D * 0.041666666666666664
    D = D * C; res += D * 0.008333333333333333
    D = D * C; res += D * 0.001388888888888889
    D = D * C; res += D * 0.0001984126984126984

    @inbounds for _ in 2:j; res = res * res; end
    res
end

@inline @fastmath function exp5(A::SMatrix{5, 5, Float64, 25}, q::Int=11)
    # Norm-Berechnung (Frobenius-Norm)
    s = 0.0
    @inbounds for i in 1:25; s += A[i]^2; end
    p = sqrt(s)
    
    # Scaling
    normiter = 2.0; j = 2
    while p > normiter
        normiter *= 2.0
        j += 1
    end

    C = A / normiter
    res = one(SMatrix{5,5,Float64,25})
    res += C
    
    # Taylor-Reihe für q Terme
    D = C
    for i in 2:(q-1)
        D = D * C
        D = D / Float64(i)
        res += D
    end

    # Squaring
    @inbounds for _ in 1:(j-1)
        res = res * res
    end
    
    return res
end

@inline @fastmath function get_A(theta::SVector{35, Float64}, avgT::SVector{12, Float64}, sumP::Float64, diam::Float64, leach::Float64)
  m3 = sumP / 1000

  #Average temperature dependence
  avgT2 = avgT .^ 2
  tem = sum(exp.(theta[22] * avgT + theta[23] * avgT2)) / 12.
  temN = sum(exp.(theta[24] * avgT + theta[25] * avgT2)) / 12.
  temH = sum(exp.(theta[26] * avgT + theta[27] * avgT2)) / 12.

  #Precipitation dependence
  tem = tem * (1.0 - exp(theta[28] * m3))
  temN = temN * (1.0 - exp(theta[29] * m3))
  temH = temH * (1.0 - exp(theta[30] * m3))

  #Size class dependence -- no effect if d == 0.0
  size_dep = (diam > 0.0) ? min(1.0, (1.0 + theta[33]*diam + theta[34]*diam^2)^theta[35]) : 1.0

  #Calculating matrix a (will work ok despite the sign of alphas)
  A = MMatrix{5, 5, Float64, 25}(undef)
  @inbounds begin
    k = tem * size_dep
    A[1, 1] = theta[1] * k
    A[2, 2] = theta[2] * k
    A[3, 3] = theta[3] * k
    A[4, 4] = theta[4] * temN * size_dep
    A[5, 5] = theta[32] * temH
    
    dAbs = abs.(@SVector [A[1], A[7], A[13], A[19]])
    idx = 5
    for i in 1:4, j in 1:4
      if i != j
        A[i, j] = theta[idx] * dAbs[j]
        idx += 1
      end
    end
    #mass flows AWEN -> H (size effect is present here)
    for j in 1:4; A[5, j] = theta[31] * dAbs[j]; end
    #Leaching (no leaching for humus) 
    if (leach < 0.0)
      lm3 = leach * m3
      for j in 1:4; A[j, j] += lm3; end
    end
  end
  SMatrix(A)
end

@inline @fastmath function get_next_timestep(A::SMatrix{5, 5, Float64}, init::SVector{5, Float64}, infall::SVector{5, Float64}, t::Float64 = 1., q::Int = 11)
  #A \ (exp(A * t) * (A * init + infall) - infall)
  A \ (exp5(A * t, q) * (A * init + infall) - infall)
end

function get_theta(theta = [0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25])
  @assert length(theta) == 35
  tabs = [1:4; 32; 35]
  theta[tabs] .= -abs.(theta[tabs])
  SVector{35, Float64}(theta)
end

function get_spin(A::SMatrix{5, 5, Float64}, infall::SVector{5, Float64})
  A \ infall * -1
end
end
