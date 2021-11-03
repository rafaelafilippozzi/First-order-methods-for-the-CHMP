# First-order-methods-for-the-CHMP
Routine used in the article "First-order methods for the convex hull membership problem.

# Folders and Files

- Algorithms: In this folder are the algorithms used in the experiments.
    * `TriangleAlgorithm.m`: Triangle Algorithm with random pivot
    * `GreedyTriangleAlgorithm.m`: Greedy Triangle Algorithm
    * `AwayStepFrankWolfeAlgorithm.m`: Frank Wolfe Algorithm with Away Step
    * `SpectralProjectedGradient.m`: Spectral Projected Gradient (SPG)
    * `simplex_proj.m`: Projection onto unit simplex for SPG
    * `nonmonotoneArmijocriterion.m`: Nonmonotone Armijo Criterion for SPG.

- Functions: Functions used to create the routine from sections 6.1 and 6.2.
    * `artificialchmp.m`: Generates artificial instances of CHMP
    * `artificialfeasibilityproblem.m`: Generates artificial instances of linear programming feasibility problem
    * `generateArandom.m`: Generates A randomly according to a uniform distribution on the unit ball of Rm
    * `generatesgraphics.m`: Generates the graphics in section 6.1

- Experiment-data-6-1: Data included in section 6.1
- Experiment-data-6-2: Data included in section 6.2

    * `ExperimentsSection6_1.m`: Routine used to generate the results of section 6.1
    * `ExperimentsSection6_2.m`: Routine used to generate the results of section 6.2
