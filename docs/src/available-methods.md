# Available integration schemes
Currently only a handful of integration schemes are supported. These are:

  * [`RK4`](@ref), a classical fourth order Runge-Kutta method, for non-stiff problems
  * a set of low-storage IMEX method developed by Daniele Cavaglieri and Thomas Bewley at UCSD [^1] for stiff problems, where stiffness arises from the linear term
    * [`CB3R2R2`](@ref), a second order accurate method
    * [`CB3R2R3e`](@ref) and [`CB3R2R3c`](@ref), two third order accurate methods 
    * [`CB4R3R`](@ref), a fourth order accurate method
  * [`CNRK2`](@ref), a classical second order Crank-Nicholson-Runge-Kutta method for stiff problems

## Usage
All methods have constructors with similar signatures.

### Standard problems and coupled systems
For standard problems including coupled systems, the constructor accepts an object of the type used to represent the state (see [Flows.jl Defining the integration scheme](@ref)). For instance, to construct an [`RK4`](@ref) object for a system defined by a $4\times 4$ Julia `Matrix` type
```julia
m = RK4(zeros(4, 4))
```

For [Flows.jl Coupled dynamical systems](@ref), the object passed two the constructor should be a [`Coupled`](@ref) object, consistent with the state type.

### Linearised equations 
For linearised equations marched over an `AbstractStorage` the constructor accepts an additional argument, depending on whether the equations correspond to a forward or adjoint problem. An `RK4` method for the forward problem with the same type as discussed before is constructed as
```julia
m = RK4(zeros(4, 4), ContinuousMode(false))
```
while for the adjoint problem
```julia
m = RK4(zeros(4, 4), ContinuousMode(true))
```
The object `ContinuousMode` signals that we are solving a continuous approximation of the linearised equations. Discrete adjoint solvers are planned, but not implemented yet.

## References
[^1] Cavaglieri, D. and Bewley, T., 2015. Low-storage implicit/explicit Runge–Kutta schemes for the simulation of stiff high-dimensional ODE systems. Journal of Computational Physics, 286, pp.172-193.