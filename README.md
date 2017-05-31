# IMEXRKCB.jl

This repository contains an efficient Julia implementation of three of the eight low-storage implicit-explicit Runge-Kutta time integration schemes recently developed by Daniele Cavaglieri & Tom Bewley ([JCP, 286, pp. 172-193, 2015](http://www.sciencedirect.com/science/article/pii/S0021999115000352)). These schemes are particularly suited for high-dimensional systems with stiff behaviour, such as systems arising from the discretisations of the Navier-Stokes equations for turbulent flow simulation. 

More specifically, the third-order accurate schemes __IMEXRKCB3e__ and __IMEXRKCB3c__ are available using the three-register implementation of [2R] IMEXRK schemes, while the fourth-order accurate __IMEXRKCB4__ scheme is available using the four-register implementation for [3R] schemes. These three also have embedded schemes for local error estimation that can be enabled if required. However, no step-size control feature is currently implemented. 

## Usage ##

### Definition of the state ###
We use duck-typing throughout to impose as little restriction as possible on user types. The state type can be anything, not simply a vector, but it must implement a basic indexing interface. For example, assume the state is represented by
```julia
struct FooBar
    data        # actual state variables, e.g. a vector, matrix, shared array...
    otherfield1 # some other stuff
    otherfield2 # more stuff
end
``` 
Then `FooBar` should implement
```julia
Base.getindex( f::FooBar,      i::Int)
Base.setindex!(f::FooBar, val, i::Int)
Base.similar(f::FooBar)
Base.eachindex(f::FooBar)
```
where the first two are used to loop over the state variables, using the information from `eachindex`, and `similar` is used 
to allocate internal storage. If the state is an `Array{T, N}`, e.g. a `Vector` than all of this is not needed.

### Defining the operators ###
The dynamical system to be integrated is defined by its stiff and non stiff components, where the former is assumed to be linear. To use the integrators, users needs to define appropriate behaviour for the types used to represent the two components. The stiff part, denoted by `A` (of type, say, `Atype`) should extend two methods, one from `Base`, the other from the `IMEXRKCB` package:
```
Base.A_mul_B!(out::FooBar, A::Atype, x::FooBar)
IMEXRKCB.ImcA!(A::Atype, c::Real, y::FooBar, z::FooBar)
```
The first calculates the action of the discretised linear operator `A` on `x` and stores the result in `out`. The latter solves
for `z` such that `(I-cA)z = y`, where `I` is the identity operator and `c` is a real constant. The result is written in `z`. The variable `A` can contain all the storage required to perform these two operations, e.g. storage for LU factorisation, etc..

The non stiff part `g` should be callable, with signature
```julia
g(t::Real, x::FooBar, dxdt::FooBar)
```
where `dxdt` is the time derivative of the state associated with the non stiff terms, `x` is the current state, and `t` is time. `dxdt` will be overwritten.

### Scheme definition ###
To use the integrator, one needs first to define the scheme that will be used. The two third order integration schemes are constructed using 
```julia
scheme = IMEXRK3R2R(IMEXRKCB3e, x, false),
```
where the first argument is a tableau (IMEXRKCB3e/IMEXRKCB3e for three register implementations of 2R scheme), `x` is an object of type `FooBar` (used to allocate internal storage for the scheme) and the last argument is a boolean flag to activate the embedded scheme. When this is activated, the state calculate by the embedded scheme at the end of a time step can be read from `scheme.storage[5]`. Note that in the construction of the scheme, `Base.similar(::FooBar)` is called to allocate storage.

The fourth order scheme (using a four register implementation of the 3R scheme) is constructed using 
```julia
scheme = IMEXRK4R3R(IMEXRKCB4, x, false),
```
*Do not mix implementations and schemes*, as it leads to wrong integration result.

### Time stepping ###
A time step is calculated with
```julia
step!(scheme, g, A, t, Δt, x)
```
where `x` will be overwritten with the state at the end of the time step, from time `t` to time `t+Δt`. Note that 
this function does not allocate extra memory. All the required storage for the step is allocated at the time of 
constructing `scheme`. This is the only building block provided in this package. More advanced integrators can be built using this low-level primitive.
