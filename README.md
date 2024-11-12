# Enclosure-of-Lyapunov-exponents

This repository contains the Julia code asoociated with the paper ["Rigorous enclosure of Lyapunov exponents of stochastic flows"](https://arxiv.org/abs/2411.07064).

You will need to install the relevant packages from the julia REPL with the command

```julia
    using Pkg; Pkg.add(" BandedMatrices, BlockArrays, BlockBandedMatrices, Combinatorics, IntervalArithmetic, LinearAlgebra, LaTeXStrings, Plots, Polynomials, PolynomialRoots, Serialization, SparseArrays")
```

This implementation crucially relies on the package [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) for rigorously controlling rounding errors.
