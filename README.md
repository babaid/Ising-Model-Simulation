# Ising Model Simulation
A project cointaining a simulation of the Ising modell (known from theoretical physics), using advanced concepts such as Markov Chains and Monte Carlo simulations.
______________
## Introduction
The Ising model is a rather simple model in statistical physics describing ferromagnetism. In the model ferromagnetism arrises from the interaction of the spins of electrons. These interactions are described by the following Hamiltonian:

$$
H = -J \sum_{\langle i, j \rangle}s_i s_j - B\sum_i s_i
$$

where in the sum $ \sum_{\langle i, j \rangle}$ only the interaction between nearest neighbors is taken into account. To deal with the finite dimension of the system simulated, periodic boundary conditions are introduced, so one usually speaks of Ising chains:

$$
s_{i+N} = s_i
$$

Firstly we will look at situations with B = 0, since most of the questions we are interested in can still be answered and calculations are a bit simplier.

## Simulation? What and how? - Metropolis-Hatings algorithm
Since we are speaking in terms of quantum mechanics, statistics is going to be involved and the concept that systems try to occupy the lowest energetic states. The whole procedure involves calculating the energy for a current state, then flipping a random spin, calculating the new energy of the system. If the energy difference is negative then the probability that the spin-flip happens is 1, if its positive one calculates a transition probability for the spin-flip and samples from the according probability distribution, to determine wether the spin will be flipped or not.
This method is called Metropolis-Hastings algorithm, which is a pretty combinatorial method, mostly used if one wants to sample from an unknown, very complicated distribution.
Which spin will be flipped, is going to be determined by a random walk, wich is probably the simplest, and fastest and computationally inexpensivest way of doing this.

The observables of interest are going to be calculated using the Monte Carlo technique:

$$
A = \int_V dx f(x) = \lim_{n\to\infty}\frac{1}{N}\sum_{i=1}^{N}f(x_i)
$$


______________________________________________________________________________


This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Ising-Model-Simulation

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
