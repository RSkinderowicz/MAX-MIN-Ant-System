# MAX-MIN Ant System (MMAS)

This is a single-file implementation of the MAX-MIN Ant System (MMAS)
as described in:
*Stützle, Thomas, and Holger H. Hoos. "MAX-MIN Ant System." Future generation
computer systems 16.8 (2000): 889-914.*

It is intended for educational purposes and lacks generality. It uses by
default so-called _candidate lists_ to speed up the search process.  The
parameters are set to values suggested in the paper. No pheromone resetting is
implemented, and only the iteration best solution is used when depositing
pheromone. For the deeper understanding of how the MMAS works, please, refer to
the paper.

It supports (partially) loading of (A)TSP instances from the well-known [TSPLIB] repository.


## Building

Run `make` from the command line -- by default, it requires GCC version with
C++11 support but any other C++ compiler with C++11 support should be fine.
The default target is `release,` i.e., optimized version but it can be changed
to `debug.`


## Running

By default, it tries to load `kroA100.tsp.` TSP instance from the current directory
however, the instance path can be given as a parameter.

As the algorithm is being executed, the program outputs information about
each improved solution found.
By default, the algorithm is terminated after ants have constructed 1 million
solutions.


Example output:

    Read line: NAME: kroA100
    Read line: TYPE: TSP
    Read line: COMMENT: 100-city problem A (Krolak/Felts/Nelson)
    Read line: DIMENSION: 100
    Read line: EDGE_WEIGHT_TYPE : EUC_2D
    Read line: NODE_COORD_SECTION
    Greedy solution cost: 27807
    initial_pheromone: 3.59622e-07
    New best solution found with cost: 27696 at iteration 0
    New best solution found with cost: 26399 at iteration 0
    New best solution found with cost: 25311 at iteration 0
    New best solution found with cost: 25126 at iteration 1
    New best solution found with cost: 25065 at iteration 1
    New best solution found with cost: 24807 at iteration 2
    New best solution found with cost: 24694 at iteration 6
    ...
    New best solution found with cost: 21368 at iteration 58240
    New best solution found with cost: 21282 at iteration 58241



[TSPLIB]: <http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/>
