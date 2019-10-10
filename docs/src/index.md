# Flows.jl
A flow-like API to solve differential equations. 

## Rationale
Many numerical algorithms in dynamical systems theory require the action of a flow operator associated to a dynamical system, rather than a solve-an-initial-value-problem-and-store-its-solution approach. This package provides an API geared towards this needs, without the ambition to be a general purpose package.

## Installation
This package is not registerd in the Julia's `METADATA.jl` but can be installed in the Julia REPL. First, hit `]` to enter package mode, then type
```julia
add https://github.com/gasagna/Flows.jl
```

## Table of Contents
```@contents
Pages = ["quickstart.md", "monitors.md", "available-methods.md", "coupled.md", "quadrature.md", "examples.md", "advanced.md", "api.md",]
Depth = 1
```

