# Calculating integrals

## Fundamentals
It is sometimes necessary to calculate the integral of some function $g : \mathcal{X} \rightarrow \mathbb{R}$ along a trajectory of the system, i.e. the integral
```math
I = \int_0^T g(\mathbf{x}(t))\mathrm{d}t
```
where $\mathbf{x}(t)\in\mathcal{X}$ is the system's state.

This can be achieved by adding a pure quadrature equation, i.e. augmenting the dynamical system with a further equations that gets integrated in time jointly with the main problem. In other words, the coupled system
```math
\left\{
\begin{array}{c}
  \dot{\mathbf{x}}(t) &= \mathbf{f}(t, \mathbf{x}(t))\\
  \dot{I}(t) &= g(\mathbf{x}(t))
\end{array}
\right.
```
is integrated from $0$ to $T$ with initial condition $I(0) = 0$. It is easy to see that the final condition $I(T)$  is the numerical approximation of the integral above. 

The usefulness of this approach is that the integral is obtained without doing much, at the end of the solution. In addition, the order of accuracy of the integration methods used to advance the coupled state is the order of accuracy of the estimation of the integral above.

## Example
As an example let us compute the integral of the function 
```math
g(\mathbf{x}) = x_1^2
```
for a system with three degrees of freedom arranged in a standard Julia `Vector`. First we define the quadrature equation 
```julia
g(t, x, dxdt, I, dIdt) = (dIdt[1] = x[1]^2; return dIdt)
```
Note the structure of the signature.

We then define some initial conditions, making sure `I[1]` is initially set to 
zero.
```julia
x = Float64[1.0, 3.0, 4.0]
I = Float64[0.0]
```
Note that because the augmented state must be mutable, the extra quadrature variable must be stored in a mutable container, a `Vector`. In fact, one can calculate the integral of as many functions as desired, just define the function `g` above to fill as many components of `dIdt` as required.

We then define a flow operator, by coupling the right hand side of the differential equation (assume it is defined by a julia function `f`) and the additional quadrature equation
```julia
F = flow(couple(f, g), RK4(couple(x, I)), TimeStepConstant(0.1))
```
Note how the state and the quadrature variables are coupled toghether and passed to the `RK4` constructor. This makes sure that the temporary objects in the Runge-Kutta scheme are consistent with the extra degree of freedom.

The flow operator can be now called with
```julia
F(couple(x, I), (0, 10))
```
where `x` is advanced forward in time from $t=0$ to $t=10$. Upon return `I[1]` contains an approximation of the integral
```math
\int_0^{10} x_1^2(t)\mathrm{d}t
```

!!! example
    It is possible to monitor the value of the quadrature variable along a simulation. For instance, let us calculate the cumulative average of the function $g(\mathbf{x}(t))$, the quantity
    ```math
    \bar{g}(\tau) = \frac{1}{\tau} \int_0^\tau g(\mathbf{x}(t))\mathrm{d}t
    ```

    To do so, simply define a `Monitor` object that observes the quadrature component
    ```julia
    mon = Monitor(couple(x, I), xq->xq[2])
    ```

    After numerical integration, samples of the function $\bar{g}$ can be obtained as
    ```julia
    g_bar = samples(mon)./times(mon)
    ``` 
    where the first element is undefined.
