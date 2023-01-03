# Yasso15
Yasso soil carbon model

Includes code and coefficients for the versions 15.

The core calculaton, originaly written in fortran, is rewirtten in C++, R and Julia.

# For C++ Users

The linear solver, based originaly on gaussian elimination, is repaced with Crout LU decomposition.

The originaly fixed number of Taylor terms is selectable.

To the original method, calculating matrix exponential with Taylor series, a method directly computing exp(A * t) * v using Chebyshev from Expokit was added.

A simple example shows the usage inside C++.


# For R Users

Simply download [simpleExample.r](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/15/r/simpleExample.r) showing how to use it and if you want a local copy of the function also download [y15_subroutine.r](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/15/r/y15_subroutine.r).


# For Julia Users

Simply download [simpleExample.jl](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/julia/simpleExample.jl) showing how to use it and if you want a local copy of the function also download [y15_subroutine.jl](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/15/julia/y15_subroutine.jl).
