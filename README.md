# Yasso
Yasso Soil Carbon Model

Includes code and coefficients for the versions 07, 15 and 20.

Code and data collected from [Official releases of the soil carbon model YASSO](https://github.com/YASSOmodel), [Jari Liski](https://github.com/JariLiski), [Jussi Kollin](https://github.com/JussiKollin), [Adobi Elizabeth Okwuosa](https://github.com/IKWENZI/Yasso15) and [Siagosheva](https://github.com/siagosheva/yasso07ui).

The core calculaton, originaly written in fortran, is rewirtten in C++, R, Julia and Python.

# For C++ Users

The linear solver, based originaly on gaussian elimination, is repaced with Crout LU decomposition.

The originaly fixed number of Taylor terms is selectable.

To the original method, calculating matrix exponential with Taylor series, a method directly computing exp(A * t) * v using Chebyshev from Expokit was added.

A simple example shows the usage inside C++.

# For Rust Users

Download the source and use "cargo run --release" to start the example.

# For R Users

Simply download [simpleExample.r](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/r/simpleExample.r) showing how to use it and if you want a local copy of the function also download [y20_subroutine.r](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/r/y20_subroutine.r).


# For Julia Users

Simply download [simpleExample.jl](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/julia/simpleExample.jl) showing how to use it and if you want a local copy of the function also download [y20_subroutine.jl](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/julia/y20_subroutine.jl).

# For Python Users

Download [simpleExample.jl](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/python/simpleExample.py) showing how to use it. A local copy of the function needs to be downloades in addition to the same folder [y20_subroutine.py](https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/python/y20_subroutine.py).


# Speed comparison

The here provided implementations for YASSO20 can make n simulations steps per second in the different languages on an AMD Ryzen 7 PRO 5850U:

1. C++: 1'821'407 iterations per Second (2'593'899 when using fast-math)
3. Julia: 1'744'245 iterations per Second
2. Fortran: 1'642'435 iterations per Second
3. Rust: 1'236'705 iterations per Second
4. R: 4'462 iterations per Second
5. Python: 4'030 iterations per Second

C++, Julia, Fortran or Rust are about 300-600 times faster than R or Python.