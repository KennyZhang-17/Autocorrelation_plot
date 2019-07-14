# AutoCorrelation.jl

![](https://github.com/KennyZhang-17/Autocorrelation_plot/blob/master/theoretical_vs_sample_autocorrelation.png)

## Installation
To install Julia, go to https://julialang.org/

* The package is **not** yet registered in `METADATA.jl` and can only be installed with `using Pkg; Pkg.add("https://github.com/KennyZhang-17/Autocorrelation_plot")`.

* Another way is to `clone https://github.com/KennyZhang-17/Autocorrelation_plot` first, and add the following in `~/.julia/config/startup.jl`
(create if doesn't exist):
```
push!(LOAD_PATH, "<path to dir of this cloned repo>")
```

## Basic Usage
```julia
julia> using AutoCorrelation

julia> kenny(0, 500, 10, 1, 0.1)
1.0622498624093526 #this is (expected - mean^2) / variance
```

Checkout `demo` for commented usage.
