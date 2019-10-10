# Full public API

## Flow operator API
The basic building block of this package is the [`Flows.Flow`](@ref) object, a discrete approximation of the flow of a dynamical system. Here is a list of possible constructors.
```@docs
flow
```

Objects of type `Flow` satisfy a callable interface, with additional arguments possible.
```@docs
Flows.Flow
```

## Monitor API
```@docs
Monitor
reset!
times
samples
```

## Coupled API
```@docs
Coupled
couple
couplecopy
getindex
similar
size
```

## Storage API
```@docs
RAMStorage
period
isperiodic
timespan
```

## Integration methods API
```@docs
RK4
CNRK2
CB3R2R2
CB3R2R3c
CB3R2R3e
CB4R3R4
ImcA!
```

## Time Stepping API
```@docs
TimeStepConstant
TimeStepFromStorage
AbstractTimeStepFromHook
```