# nn-verification 

The companion repository to the paper 

**Speeding Up Neural Network Robustness Verification via Algorithm Configuration and an Optimised Mixed Integer Linear Programming Solver Portfolio**, Matthias König, Holger H Hoos, Jan N van Rijn, to be published.

This repository provides

- an example for the configuration of the Venus verifier [1]
- an example for the configuration of the MIPVerify verifier [2]

# Usage

Required tools:

- Gurobi v9.0.1 (www.gurobi.com)
- SMAC v2.10.3 (http://www.cs.ubc.ca/labs/beta/Projects/SMAC/)
- Hydra v1.1 (http://www.cs.ubc.ca/labs/beta/Projects/Hydra/)
 
For our configuration experiments, we used the framework provided in the 

- Algorithm Configuration Library 2.0 (https://bitbucket.org/mlindauer/aclib2/src/master/)

## Install Gurobi
```
wget https://packages.gurobi.com/9.0/gurobi9.0.1_linux64.tar.gz
tar xvfz gurobi9.0.1_linux64.tar.gz
cd gurobi901/linux64
python3 setup.py install
```
Next, add the following to the .bashrc file:
```
export GUROBI_HOME="your/path/here/gurobi901/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
```

# References
[1] Botoeva E, Kouvaros P, Kronqvist J, Lomuscio A, Misener R (2020) Efficient Verification of ReLU-based Neural Networks via Dependency Analysis. In: Proceedings of The Thirty-Fourth AAAI Conference on Artificial Intelligence (AAAI20), pp 3291–3299
[2] Tjeng V, Xiao K, Tedrake R (2019) Evaluating Robustness of Neural Networks with Mixed Integer Programming. In: Proceedings of the 7th International Conference on Learning Representations (ICLR 2019)
