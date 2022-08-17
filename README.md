# First-order-methods-for-the-CHMP
Routine used in the article "First-order methods for the convex hull membership problem.

# Folders and Files

* `ExperimentsSection6_1.m`: Routine used to generate the results of section 6.1

* `ExperimentsSection6_2.m`: Routine used to generate the results of section 6.2

* `mnistchmp.m`: Routine used to generate the results of section 6.3

- Algorithms: In this folder are the algorithms used in the experiments.
    * `TriangleAlgorithm.m`: Triangle Algorithm with random pivot
    * `GreedyTriangleAlgorithm.m`: Greedy Triangle Algorithm
    * `AwayStepFrankWolfeAlgorithm.m`: Frank Wolfe Algorithm with Away Step
    * `SpectralProjectedGradient.m`: Spectral Projected Gradient (SPG)
    * `simplex_proj.m`: Projection onto unit simplex for SPG
    * `nonmonotoneArmijocriterion.m`: Nonmonotone Armijo Criterion for SPG.

- Functions: auxiliary functions used by ExperimentsSection6_x.m
    * `artificialchmp.m`: Generates artificial instances of CHMP
    * `artificialfeasibilityproblem.m`: Generates artificial instances of linear programming feasibility problem
    * `generateArandom.m`: Generates A randomly according to a uniform distribution on the unit ball of Rm
    * `generatesgraphics.m`: Generates the graphics in section 6.1

- Experiment-data-6-1: Data included in section 6.1
- Experiment-data-6-2: Data included in section 6.2
- Experiment-data-6-3: Data included in section 6.3

To reproduce the experiments from Section 6.1 as they are in the article, add the folders 'Functions' and 'Algorithms' to MATLAB path and use the routine ExperimentsSection6_1.m. A folder called "visitordata-6-1" will be created with the new data/results. 
The same applies to the experiments of Section 6.2 and Section 6.3. 
Note: For the experiments in section 6.3 it is necessary to extract the MNIST data set.

