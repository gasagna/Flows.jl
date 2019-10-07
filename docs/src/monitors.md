# Monitor objects
A `Flow` operator simply maps a state vector forward in time by some specified amount. It operates in place, and does not store or record anything during the trajectory. However, it is sometimes useful to record some quantity, for instance one of the degrees of freedom, or maybe some integral quantity along a trajectory. This can be achieved by using a `Monitor` object.

## Basic usage
The constructor of the `Monitor` type has the signature
```julia
Monitor(x::X, g::Union{Callable, Function})
```
The first argument is an object of some user defined type, say `X`, the same type used to represent the system's state. The second argument is a function or callable object that we use to 'observe' the state along the simulation. It must accept a single argument of type `X`, must have the signature
```julia
g(::X)
```
and can return anything. 

!!! example
    This example demonstrates how to define an object to monitor the first state of a dynamical system with three degrees of freedom.
    ```julia
    mon = Monitor(zeros(3), x->x[1])
    ```

    Note how the second argument is simply an anonymous function that extracts the first element. A more elegant approach is also
    ```julia
    mon = Monitor(zeros(3), first)
    ```

In practice, in the constructor, the function is called on the first argument, the type of the output is analysed and storage to hold more elements of the same type is allocated.

To monitor an observable during a trajectory, the monitor object can be passed as a third argument to a `Flow` object. During the integration, a sample of the observable is taken t the end of every time step, including one sample at the beginning of the trajectory.

!!! example
    For instance, assume `F` is a `Flow` object for the Lorenz equations and we want to monitor the norm of the state vector over a short trajectory from $t=0$ to $t=1$. This can be achieved by
    ```julia
    mon = Monitor(zeros(3), norm)
    F(x, (0, 1), mon)
    ```
    
At the end of the integration, the content of the `Monitor` object `mon` can be accessed by two helper functions. The first
```julia
samples(mon)
```
returns a Julia `Vector` with samples of the observed function, while 
```julia
times(mon)
```
returns a `Vector` containing the times whan the samples are taken. This is quite useful for plotting the some observable of system as a function of time.

!!! note
    The observable function can really return anything. For instance, if we want to observe the state, we can define a monitor with the `copy` function.
    ```julia
    mon = Monitor(zeros(3), copy)
    F(x, (0, 1), mon)
    # samples(mon) now contains the state at the end of every time step
    ```

    If we want to monitor more quantities, we can pass a function that returns a `Tuple`, like so
    ```julia
    mon = Monitor(zeros(3), x->(x[1], x[2]^2))
    F(x, (0, 1), mon)
    ```
    so that `samples(mon)` returns a vector of `Tuple`s.


## Advanced usage
The behaviour of `Monitor` object can be customised more finely. Consult the API page for more details.