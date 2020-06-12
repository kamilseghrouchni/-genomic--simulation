# Genetics Project- Team 4 (EPFL)

### Aim
Simulate population genetics.

The program developed aims to simulate the Wright Fisher model. This consists
of having random alleles that have the same probability of integrating
the next generation. The simulation is a model of the Genetic drift within
a population.

An additional feature of the simulation is to have a more selective allele
distribution by taking into consideration specific factors such as: migration,
bottleneck, mutations and selection that will create a more realistic simulation of population
genetics.

As such, the program can be launched by using various modes. Each mode corresponds to a
specific model that will take into consideration one of the following factors :
migration, selection, mutation and time dependent population size (bottleneck
effect). Note that the mutation mode requires additional precision from the user
to decide on a specific mutation model (Cantor, Felsenstein, Kimura) and the
migration mode requires the pattern followed by the sub-groups as well as the allele
repartition (random or input by the user).

The program can be launched using two different methods:
  1- using only an input file (.txt) with the desired parameters for the simulation
  2- using the fasta and input.txt file. The input file is modified by the user to
     decide on the various simulation parameters (population size, number of generations,
     mutation marker sites, etc.).
  Note: for a simulation using mutation models, a fasta file is mandatory

Finally, the results.txt file generated for each simulation and can be used to generate
various graphs using jupyterNotebook and plotly library.

Special Feature:
The program is coded using a multithreading process. Each thread executes a single
simulation, allowing for replicas to run simultaneously. This allows for faster simulation.

How to run a simulation:
The user must use the input.txt file to choose:
- The number of generations
- The population size
- The execution mode (bottleneck, selection, migration, mutation)
- The marker sites (if using the mutation model)

If working with the mutation model, the user must precise:
- The mutation model with the appropriate probabilities

If working with the migration mode, the user must precise:
- The mutation pattern
- The migration rates

If working with the selection model, the user must precise:
- the selection rates for each allele

Presentation of the project (PPTX): <https://prezi.com/p/jfcpou3kqbcd/>
