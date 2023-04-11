# PolynomialSystems
Example input-output systems of polynomial-type.

## Lorenz63 (_Lorenz63_)

The classic 3-dimensional quadratic system used to study chaotic phenomena. 
The default parameters produce the classic butterfly simulation.  To 
generate inputs to the system, we consider applied to the first coordinate.

## Lorenz96 (_Lorenz96_)

A generalization of the Lorenz63 system to higher state dimensions.  The 
quadratic system is written as a ring, where each node has connections to
the node in front and the two nodes behind it.  A uniform control is applied 
to every state.

## Rossler (_Rossler_)

Another classic 3-dimensional quadratic system that exhibits an interesting
attractor.

## van der Pol ring (_vanderPolRing_)

A ring of van der Pol oscillators (oscillators that are weakly coupled to
their neighbors).  This cubic model is the subject of numerous synchronization
studies.  Individual controls are applied to the velocity components of selected
oscillators by in index parameter.  Each oscillator can optionally have 
different strength nonlinearities and coupling parameters.
