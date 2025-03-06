# MorphoTraj
This repository contains functions to aid in the comparison of morphogenetic shape data. Currently, it specifically contains two functions.

The first ('pairwiseDistances') calculates pairwise Euclidean distances betweeen all input individuals, generates a heatmap summary, and calculates a ranked list of pairwise comparisons in descending order from those that are closest to each other in Euclidean space.

The second ('draw.traj') performs a PCA on an input matrix (ideally representing a developmental series, though this is not required) to identify the primary axis of shape variation, calculates sum of squares distances between each point in the matrix and an equivalent equidistantly-separated point on a line representing the primary axis of shape variation, determines the direction of the trajectory by minimizing the sum of squares across all points, and draws this tractory as an arrow. It also reports the mean sum of squares distance across every point on the line as an indiciation of dispersal around the trajectory.
