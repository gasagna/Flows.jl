# Flows.jl API

# The Flow operator

The basic building block of this package is the [`Flows.Flow`](@ref) object, a discrete approximation of the flow of a dynamical system. Here is a list of possible constructors.
```@docs
flow
```

Objects of type `Flow` satisfy a callable interface, with additional arguments possible.
```@docs
Flows.Flow
```