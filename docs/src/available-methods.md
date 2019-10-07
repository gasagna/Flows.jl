# Available numerical methods
We currently support a handful of numerical methods, much less than other serious differential equations packages. More precisely we support:

- a classical fourth order Runge-Kutta method, for non-stiff equations
- second, third and fourth order accurate, low-storage IMEX methods developed by Cavalieri and Bewley (JCP 2015) for stiff problems, where stiffnes arises from the linear term.
- a classical second order Crank-Nicholson-Runge-Kutta method