# Flows.jl
A convenient flow-like API for numerical integration of differential equations.

# A brief explanation of why this package exists: a disclaimer!
There are now plenty of state-of-the-art differential equations Julia packages out there, but none satisfied by exact needs when I first looked into this. What I originally required was a package offering a flow-like API to a differential equation, as many numerical algorithms in dynamical systems theory require the `action` of the flow operator associated to a dynamical system, rather than the more usual solve-an-initial-value-problem-and-store-its-solution approach. I also needed solvers for stiff ODEs arising from discretisations of dissipative PDE (i.e. the Navier-Stokes equations governing fluid flow), and a way to solve the linearised tangent equation and their discretely-consiste adjoint version. 


# Using the package
## Available numerical methods
We currently support a handful of numerical methods, much less than other serious differential equations packages. More precisely we support:

- a classical fourth order Runge-Kutta method, for non-stiff equations
- second, third and fourth order accurate, low-storage IMEX methods developed by Cavaliery and Bewley (JCP) for stiff problems, where stiffnes arises from the linear term.

## Defining the flow of a dynamical systems
Let a dynamical system be defined by some evolution rule $f$


```julia
x0 = MyState()
F = flow(f, RK4(x0, :NORMAL), TimeStepConstant(0.1))
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
F = flow(couple(f, quadfun), RK4(couple(x0, q0), :NORMAL), TimeStepConstant(0.1))
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