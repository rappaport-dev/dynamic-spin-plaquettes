# Optimized Kinetic Monte Carlo for Spin Glass Systems

A high-performance **Julia** simulation of spin glass dynamics using **Glauber Dynamics**. This project implements an algorithmic optimization using a **Fenwick Tree** (Binary Indexed Tree) to reduce the time complexity of spin selection from $O(N)$ to $O(\log N)$.

## Engineering Features
* **Algorithmic Optimization:** Replaced linear search with a Fenwick Tree structure to handle spin-flip rate updates efficiently.
* **HPC Integration:** Includes **Slurm** submission scripts (`run.sl`) for execution on the Boston College cluster.
* **Reproducibility:** Scripts automatically log git commit hashes and diffs to ensure experimental reproducibility.

##  Files
* `simulation.jl`: Main simulation logic implementing the Fenwick Tree and MCMC steps.
* `run.sl`: Batch submission script for Slurm workloads.
* `plotenergy.jl`: Visualization utility for energy decay analysis.

## Example Usage
To run a simulation locally:
```bash
julia simulation.jl --L 32 --p 0.1 --beta 6.0 --T_max 10000
```

To submit a simulation to the cluster you can use Slurm:

```batch
sbatch run.sl
```
