# nn-verification 

The companion repository to the paper 

**Speeding Up Neural Network Robustness Verification via Algorithm Configuration and an Optimised Mixed Integer Linear Programming Solver Portfolio**, Matthias König, Holger H Hoos, Jan N van Rijn, to be published.

This repository provides

- the benchmark instances used in our configuration experiments
- the wrapper for configuring the commercial MIP solver Gurobi
- the wrapper for configuring the Gurobi solver directly through Venus.

# Usage

Required tools:

- Gurobi v9.0.1 (www.gurobi.com)
- SMAC v2.10.3 (http://www.cs.ubc.ca/labs/beta/Projects/SMAC/)
- MIPVerify v0.2.3 (https://github.com/vtjeng/MIPVerify.jl)
- Venus v1.01 (https://vas.doc.ic.ac.uk/software/neural/)
- Hydra v1.1 (http://www.cs.ubc.ca/labs/beta/Projects/Hydra/)
 
For our configuration experiments, we used the framework provided in the 

- Algorithm Configuration Library 2.0 (https://bitbucket.org/mlindauer/aclib2/src/master/)

Once you have installed all required tools, you can create add use the files from this repository within the aclib framework and create a scenario file. 
Hydra runs on top of this framework, which means that it simply uses the scenario file specified within aclib.
