# Yasso20
Yasso soil carbon model

Includes code and coefficients for the versions 20.

The core calculaton, originaly written in fortran, is rewirtten in C++, R and Julia.

# For C++ Users

The linear solver, based originaly on gaussian elimination, is repaced with Crout LU decomposition.

The originaly fixed number of Taylor terms is selectable.

To the original method, calculating matrix exponential with Taylor series, a method directly computing exp(A * t) * v using Chebyshev from Expokit was added.

A simple example shows the usage inside C++.

# For Rust Users

Download the source and use >RUSTFLAGS="-C target-cpu=native" cargo run --release< to start the example.

# For R Users

Simply download [simpleExample.r](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/r/simpleExample.r) showing how to use it and if you want a local copy of the function also download [y20_subroutine.r](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/r/y20_subroutine.r).


# For Julia Users

Simply download [simpleExample.jl](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/julia/simpleExample.jl) showing how to use it and if you want a local copy of the function also download [y20_subroutine.jl](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/julia/y20_subroutine.jl).


# For Python Users

Download [simpleExample.jl](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/python/simpleExample.py) showing how to use it. A local copy of the function needs to be downloades in addition to the same folder [y20_subroutine.py](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/python/y20_subroutine.py).
