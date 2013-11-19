TODO
====

* Change the DE optimizers and RandomSearcher over to column-major order for individuals like NES uses. This is what we will use form now on since it is Julia's default and might give some speed benefits.

* Implement search space bounding via fitness penalty instead of explicit bounding. The latter is messy and not clear what is the right way to bound. Instead we use fitness(x) = fitness(feasible(x)) + alpha * norm(x - feasible(x), 2) where alpha is 1e-3 or similar, i.e. we penalize going outside the box. This is more general and simpler. Update DE code accordingly.