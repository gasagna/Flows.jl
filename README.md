# Flows.jl
A convenient flow-like API for numerical integration of differential equations.

# A brief explanation of why this package exists: a disclaimer!
There are now plenty of state-of-the-art differential equations Julia packages out there, but none satisfied by exact needs when I first looked into this. What I originally required was a package offering a flow-like API to a differential equation, as many numerical algorithms in dynamical systems theory require the `action` of the flow operator associated to a dynamical system, rather than the more usual solve-an-initial-value-problem-and-store-its-solution approach. I also needed solvers for stiff ODEs arising from discretisations of dissipative PDE (i.e. the Navier-Stokes equations governing fluid flow), and a way to solve the linearised tangent equation and their discretely-consistent adjoint version. 

# Quick start 
## Available numerical methods
We currently support a handful of numerical methods, much less than other serious differential equations packages. More precisely we support:

- a classical fourth order Runge-Kutta method, for non-stiff equations
- second, third and fourth order accurate, low-storage IMEX methods developed by Cavalieri and Bewley (JCP 2015) for stiff problems, where stiffnes arises from the linear term.
- a classical second order Crank-Nicholson-Runge-Kutta method

## Fundamentals
Let a dynamical system be defined by some evolution rule $f$ as
```math
\dot{\mathbf{x}}(t) = \mathbf{f}(\mathbf{x}(t))
```
where $t$ is time and $x$ is the state, generically an element of a vector space $\mathcal{X}$. The flow of $f$ is the nonlinear operator $\mathbf{\Phi}^t : \mathcal{X} \mapsto \mathcal{X}$ that maps a point $\mathbf{x}_0$ at $t=t_0$ to its forward-in-time image $\mathbf{x}(t)$, i.e. $\mathbf{x}(t) = \mathbf{\Phi}^{t-t_0}(\mathbf{x}_0)$. 

This package is based on satisfying this interface. 

## Interfaces
### State objects
The package works transparently on state vectors `x` of arbitrary type, say objects of some type `X`. The requirement of expressing computations with object of type `Vector` is quite restrictive in practice and would not leverage the multiple dispatch capabilities of Julia. The only restriction on the type `X` is that it should be mutable and should support the dot-broadcasting notation for algebraic expressions including scalars and objects of the same type, e.g.
```julia
x .= 3.0.*y .+ 4.0.*z
```
for three elements `x`, `y` and `z` of type `X`. Internally, the package makes extensive use of the dot notation, e.g. when updating cache objects during a time step. In addition, the type `X` should provide methods for the Base functions
    * `Base.similar(x)`
    * `Base.copy(x)`
These are required internally, e.g. to pre-allocate create cache objects used for the rest of the time integration. The mutability restriction prevents having, e.g., immutable `StaticArrays` objects as state vectors, which would probably lead to faster code on small problems. This is less important for the typical applications of this package, e.g. turbulent flow simulations where the state is generally given by multiple large three-dimensional arrays.

### Governing equations
Given state vectors `x` as objects of some Julia type `X`, specifying dynamical system depends on whether an explicit or semi-implicit integration method is used. For explicit methods, a dynamical system is specified by a Julia function or callable object that has a method with signature
```julia
f(t::Real, x::X, dxdt::X).
```
This function is supposed to operate in place, and should modify the content of the third argument `dxdt` with the time derivative at point `x`, and additionally, time `t` for non-autonomous systems. Modifying the content of the second argument in the body of this function is undefined behaviour and generally discouraged. For implicit methods, an additional term is required for the integration of the linear stiff term. This should be an object `A` that has a method for the `LinearAlgebra` function `mul!` with signature
```julia
LinearAlgebra.mul!(out::X, A, x::X).
```
User should define a method for this function to computes the action of `A` on `x` and stores it in `out`, without modifying the content of `x`. In addition, the object `A` should define a method for the function `Flows.ImcA`, provided by this package, with signature
```julia
Flows.ImcA(A, c::Real, y::X, z::X)
```
that solves the linear system 
```math
    (A - c I)z = y
```
for some right hand side `y` of type `X` and a scalar `c` and stores the result in `z`.

## Cache objects
Performing a time step typically requires temporary objects for the intermediate computations, e.g. the internal stages of a Runge-Kutta method, where the number of temporaries depends on the integration method. In this package, these temporaries get preallocated by defining an instance of one of the subtypes of `AbstractMethod` provided by the package. These are


```julia
x0 = MyState()
F = flow(f, RK4(x0), TimeStepConstant(0.1))
```

Here, the flow operator `F` is now a Julia function that propagates a point `x0` in state space forward in time by some amount of time. For instance
```julia
F(x0, (0, 10))
```
propagates `x0` from t=0 to t=10. Note that for performance reasons, the flow operator operates in place, i.e. it overwrites its first input argument.

## Calculating integrals
It is very often necessary to calculate the integral of some function $f$ along a trajectory of the system.

To calculate such integrals we use a pure quadrature equation (), i.e. augment the system by a further equations that gets integrated in time jointly with the main problem. 

```julia
function quadfun(t, xq::Coupled, dqdt)
    x = first(xq)
    dqdt[1] = x[1]^2
    dqdt[2] = x[1]*x[2]
    return dqdt
end
```

```julia
# define an initial condition for the integration
x0 = Float64[1.0, 3.0, 4.0]
# initialise the quadrature component to zero
q0 = Float64[0.0, 0.0]

# define the augmented flow operator, by coupling the right hand side of differential equation and the additional quadrature equation
F = flow(couple(f, quadfun), RK4(couple(x0, q0)), TimeStepConstant(0.1))
```

This flow operator can be now called like
```julia
F(couple(x0, q0), (0, 10))
```
when we now propagate `x0` forward in time and simoultaneously compute the required integrals, whose value will be contained in the vector `q0` after the integration.

## Monitoring the solution
Often, one need to store a 

To do so, this package provides the `Monitor` object, which can be defined by
```julia
m = Monitor(x0, somefunction)
```
where `x0` is of the same type, and `somefunction` is a function providing the quantity we want to monitor. For instance, we might need to monitor the norm of the solution as time progresses, in which case
```julia
m = Monitor(x0, norm)
```
would suffice. The second argument can be an anonymous function. For instance, if we want to monitor the 25 state, we could define
```julia
m = Monitor(x0, x->x[25])
```

When a monitor object is defined, 

```julia
F(x0, (0, T), m)
```

and then access 