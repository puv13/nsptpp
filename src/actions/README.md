# Actions #

This folder contains implementations of different "theories", e.g., QCD or CP(N) models.

Every theory is defined in terms of an implementation of
the left Lie derivative w.r.t. the group, i.e. $\Nabla S$, where $S$ is the action.

Conventions:
  * Each "theory" has its own namespace.
  * Functor semantics, i.e. classes with call operators.
  * In each namespace, a class `Force` is defined which possesses a implementation of
    `Lattice<Expansion<G, order>> operator()(Lattice<Expansion<G, order>> const & gauge) const`