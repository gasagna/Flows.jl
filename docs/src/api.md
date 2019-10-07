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

# Monitor objects

```@docs
Monitor
reset!
times
samples
```

# Coupled objects
```@docs
Coupled
couple
couplecopy
getindex
similar
size
```

# Storage objects
```@docs
RAMStorage
period
isperiodic
timespan
```

# Semi-implicit methods
```@docs
ImcA!
```

# Time Stepping
```@docs
TimeStepConstant
TimeStepFromStorage
AbstractTimeStepFromHook
```