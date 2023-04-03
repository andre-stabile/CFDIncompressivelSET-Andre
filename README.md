2D Finite Element solver for fluid-structure interaction with incompressible flows and rigid body structures. Uses T2T2 elements and SUPG, PSPG and LSIC stabilizations and an Arbitrary Lagrangian Eulerian (ALE) formulation. It follows from an initial code for incompressible flows by raksanches.

Three branches are available:
* Uncoupled-problem: the structure is excited by the fluid flow, but no coupling takes place, in a way that the fluid is not affected by the structure
* Weakly-coupled-problem: a block-diagonal weak coupling approach is adopted, in such a way that the structure's position and mesh update are done outside the Newton-Raphson loop of the fluid-flow
* Strongly-coupled-problem: a block-diagonal strong coupling approach is adopted, in such a way that the structure's position and mesh update are done inside the Newton-Raphson loop of the fluid-flow
