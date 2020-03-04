# In this file a light R-function of Yasso07 is created

# Copyright (C) <2017>  <Finnish Meteorological Institute>
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Based on Yasso07 description Tuomi & Liski 17.3.2008
# Created by Taru Palosuo in December 2011


#  Instructions  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

# 1) first run the source code for the function with "source(....r)"
# 2) then you can use the yasso07-function by just calling it yasso07.light(..)
# 3) Input needed for the function is the same as for the matrix version:
#        1. MeanTemperature - mean annual temperatures [C]
#        2. TemperatureAmplitude - temperature amplitudes [C] 
#        3. Precipitation - annual precipiations [mm]
#        4. InitialCPool - initial C pools of model compartments, length 5, [any unit]
#        5. LitterInput - mean litter input, 5 columns AWENH, [any unit]
#        6. WoodySize - size of woody litter (for non-woody litter this is 0)
#        7. Yasso07Parameters - these in the format applied in the fortran version, length 44
#        8. Time - simulation time 

# NOTE that this function eats only one type of material at the time. So, non-woody and different woody litter
# materials needs to be calculated separately.

# The output of the function is the matrix AWENH compartments at the given time since the simulation start


# Basics  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

# additional R libraries (as needed)
library(Matrix)


# Function definition   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

# The input for this function is provided as vectors, row = year

yasso07.light =  function(MeanTemperature, TemperatureAmplitude,Precipitation,InitialCPool,LitterInput,WoodySize,Yasso07Parameters,SimulationTime) {

#    MeanTemperature,...                % Mean annual temperature, C
#    TemperatureAmplitude,...           % (T_max-T_min)/2, C
#    Precipitation,...                  % Annual rainfall, mm
#    InitialCPool,...                   % AWENH, kg
#    LitterInput,...                    % AWENH, kg/a

MT=MeanTemperature
TA=TemperatureAmplitude
PR=Precipitation
PR=PR/1000;               # conversion from mm to meters

LC=InitialCPool
LI=LitterInput

PA=Yasso07Parameters
WS=WoodySize;              

TI = SimulationTime

alfa=c(-PA[1], -PA[2], -PA[3], -PA[4], -PA[35])   # Vector of decomposition rates

# Creating the matrix A_p (here called p)

row1 = c(-1, PA[5], PA[6], PA[7], 0)
row2 = c(PA[8], -1, PA[9], PA[10], 0)
row3 = c(PA[11], PA[12], -1, PA[13], 0)
row4 = c(PA[14], PA[15], PA[16], -1, 0)
row5 = c(PA[36], PA[36], PA[36], PA[36], -1)

p = matrix(c(row1, row2, row3, row4, row5), 5, 5, byrow=T)

# temperature dependence parameters
beta1 = PA[17]
beta2 = PA[18]
gamma = PA[26]

# Woody litter size dependence parameters
delta1 = PA[39]
delta2 = PA[40]
r = PA[41]

T1=MT+4*TA/pi*(1/sqrt(2)-1)          # Eq. 2.4 in model description
T2=MT-4*TA/(sqrt(2)*pi)              # Eq. 2.5 in model description
T3=MT+4*TA/pi*(1-1/sqrt(2))          # Eq. 2.6 in model description
T4=MT+4*TA/(sqrt(2)*pi)              # Eq. 2.7 in model description 

# k following Eq. 3 in Tuomi et al. 2009. Eco.Mod. 220: 3362-3371
k=alfa*mean(exp(beta1*c(T1,T2,T3,T4)+beta2*(c(T1,T2,T3,T4)^2))*(1-exp(gamma*PR)))     

# the effect of wl size as in Eq. 3.1 in model description
k = c(k[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),k[5])
    
A=p%*%diag(k)                             # Matrix multiplication in R: %*%

#analytical solution as in Eq. 1.3 in model description
LC = as.array(solve(A)%*% (expm(A*TI)%*%(A%*%InitialCPool+LI)-LI))

LC  # prints this as a result of a function
    
}  # end of yasso07 function

