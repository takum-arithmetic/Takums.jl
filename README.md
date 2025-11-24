# Takums.jl

[![Build Status](https://github.com/takum-arithmetic/Takums.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/takum-arithmetic/Takums.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This package implements [takum arithmetic](https://arxiv.org/abs/2404.18603), a new tapered
precision machine number format. Eight new data types are defined, namely
`Takum8`, `Takum16`, `Takum32` and `Takum64` for the (floating-point) takums,
and `LogTakum8`, `LogTakum16`, `LogTakum32` and `LogTakum64` for the
logarithmic takums, a logarithmic number system. Internally the reference
C99 library [libtakum](https://github.com/takum-arithmetic/libtakum) is used
for the actual computations for the most part.

Using this package one is able to evaluate takums for real-world applications
in terms of precision. Given it is a software implementation the performance
is overall worse than the respective usual IEEE 754 floating-point hardware
implementations, but it is sufficient to evaluate the formats for reasonably
sized numerical benchmarks.

## Usage

The eight takum number types `Takum8`, `Takum16`, `Takum32`, `Takum64`
`LogTakum8`, `LogTakum16`, `LogTakum32`, `LogTakum64` have been implemented
to behave as much as any built-in floating-point type. They are subtypes
of `AbstractFloat` even though logarithmic takums themselves are strictly
speaking not a floating-point number format, but a logarithmic number system.
However, as the majority of numerical code is written to accept `AbstractFloat`s
rather than `Real`s and takums share many properties of a typical floating-point
number system, this decision was made deliberately and for good reasons.

```julia
julia> using Takums

julia> x = Takum8(8.4)
Takum8(8)

julia> sqrt(x)
Takum8(3)
```

A wide range of functions are implemented for the number types. If you find
a floating-point function you need that is not yet implemented, please raise
an issue.

## Author and License

Takums.jl is developed by Laslo Hunhold and licensed under the ISC
license. See LICENSE for copyright and license details.
