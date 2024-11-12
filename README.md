# Enclosure-of-Lyapunov-exponents

This repository contains the Julia code asoociated with the paper ["Constructive proofs for some semilinear PDEs on $H^2(e^{|x|^2/4},\mathbb{R}^d)$"]([https://arxiv.org/abs/2404.04054](https://arxiv.org/abs/2411.07064)).

You will need to install the relevant packages from the julia REPL with the command

```julia
    using Pkg; Pkg.add(" BandedMatrices, BlockArrays, BlockBandedMatrices, Combinatorics, IntervalArithmetic, LinearAlgebra, LaTeXStrings, Plots, Polynomials, PolynomialRoots, Serialization, SparseArrays")
```

This implementation crucially relies on the package [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) for rigorously controlling rounding errors.
