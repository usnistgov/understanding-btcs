Refined convergence
-------------------

Surprisingly, we get similar approximate convergence rates for numerical experiments
for $\kappa=1$ and $\kappa=-1/2$.

![Convergence kappa nonzero](./refined-convergence-small-ks.png)

How come the error in the slope of the line is much larger?
Could it be because we are comparing a simulation with $n$
points to a simulation with $2n$ points (or so).
Would the error in the slope decrease if we compared only
to a simulation with a very large number of points?

When we consider larger (and positive) $\kappa$,
we still get the same convergence rates.
However, when $\kappa$ is very large, (25.6 in our case),
we find it challenging to observe the convergence rate
numerically.

![Convergence kappa large](./refined-convergence-large-positiveks.png)

TODO
----

1.  Change description of $\tau_{\star}$ in paper

2.  Write up math behind mesh refinement analysis

3.  Consider questions above

4.  Figure out what is going on with $\kappa=25.6$. 
    Consider trying to perform numerical experiments with some intermediate $\kappa$ first, say $\kappa = 12$.
