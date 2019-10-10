# Advanced features

## Call dependencies
When solving several coupled systems, the default settings are such that the signature of the first, second, third and so on functions should be
```julia
f1(t, x1, dx1dt)
f2(t, x1, dx1dt, x2, dx2dt)
f3(t, x1, dx1dt, x2, dx2dt, x3, dx3dt)
f4(t, x1, dx1dt, x2, dx2dt, x3, dx3dt, x4, dx4dt)
```
This is not always desirable. For instance, assume we need to solve the problem
```math
\tag{1}\label{eq}
\left\{
\begin{array}{cc}
  \dot{\mathbf{x}}(t)   =& \mathbf{f}(t, \mathbf{x}(t))\\
  \dot{\mathbf{y}}_1(t) =& \mathbf{g}_1(t, \mathbf{x}(t), \mathbf{y}_1(t))\\
  \dot{\mathbf{y}}_2(t) =& \mathbf{g}_2(t, \mathbf{x}(t), \mathbf{y}_2(t))\\
  \dot{\mathbf{y}}_3(t) =& \mathbf{g}_3(t, \mathbf{x}(t), \mathbf{y}_3(t))\\
\end{array}
\right.
```
where we might wish to express the structure for the function calls. 

This can be achieved in this package by using a [`CallDependency`](@ref) object. The constructor accepts as many tuples of integers as many equations we need to solve. Each tuple specifies what states are required in the signature. For instance, for the example ($\ref{eq}$) the correct call dependency specification is 
```julia
deps = CallDependency((1,), (1, 2), (1, 3), (1, 4))
```

Clearly, the signatures of the functions `f`, `g1`, `g2`, and `g3` must be consistent with this specification, i.e. should be
```julia
f(t, x, dxdt)
g1(t, x, dxdt, y1, dy1dt)
g2(t, x, dxdt, y2, dy2dt)
g3(t, x, dxdt, y3, dy3dt)
```

The `deps` object is then passed as an additional argument to the constructor of the [`Flows.Flow`](@ref) object. For instance, for an explicit method
```julia
F = flow(couple(f, g1, g2, g3), deps, otherargs...)
```
