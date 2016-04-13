# sapo
Reachability Computaion and Parameter Sythesis for Polynomial Systems

Sapo is a C++ prototype tool designed to solve the reachability and parameter synthesis problem
for polynomial dynamical systems. The dynamics of the considered system are defined as discrete-time polynomials,
the sets reachable by the system can be represented with boxes,
parallelotopes (the n-dimensional generalization of parallelograms),
or parallelotope bundles (sets of parallelotopes whose intersections symbolically represent a polytope),
while the parameter sets are represented by polytopes.

The reachability computation and parameters refinement are carried out representing
the polynomials in Bernstein form and transforming the problems into Linear Programs.
