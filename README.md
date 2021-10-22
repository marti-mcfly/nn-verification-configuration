# nn-verification 

The companion repository to the paper 

**Speeding Up Neural Network Robustness Verification via Algorithm Configuration and an Optimised Mixed Integer Linear Programming Solver Portfolio**, Matthias König, Holger H Hoos, Jan N van Rijn, to be published. 

**Abstract:** Despite their great success in recent years, neural networks have been found to be vulnerable to adversarial attacks. These attacks are often based on slight perturbations of given inputs that cause them to be misclassified. Several methods have been proposed to formally prove robustness of a given network against such attacks. However, these methods typically give rise to high computational demands, which severely limit their scalability. Recent state-of-the-art approaches state the verification task as a minimisation problem, which is formulated and solved as a mixed-integer linear programming (MIP) problem. We extend this approach by leveraging automated algorithm configuration techniques and, more specifically, construct a portfolio of MIP solver configurations optimised for the neural network verification task. We test this approach on two recent, state-of-the-art MIP-based verification engines, MIPVerify and Venus, and achieve substantial improvements in CPU time by average factors of up to 4.7 and 10.3, respectively.

See https://aisecure-workshop.github.io/aml-iclr2021/papers/36.pdf for an extended abstract.

This repository provides

- an example for the configuration of the Venus verifier [1];
- an example for the configuration of the MIPVerify verifier [2].

In both cases, the considered network is an MNIST classifier designed for robustness taken from [3].

# Usage

## Required tools:

- Gurobi v9.0.1 (www.gurobi.com)
- SMAC v2.10.3 (http://www.cs.ubc.ca/labs/beta/Projects/SMAC/)
- Hydra v1.1 (http://www.cs.ubc.ca/labs/beta/Projects/Hydra/)

Verification engines:

- Venus v1.10 (https://vas.doc.ic.ac.uk/software/neural/)
- MIPVerify v0.2.3 (https://github.com/vtjeng/MIPVerify.jl)
 
For our configuration experiments, we used the framework provided in the 

- Algorithm Configuration Library 2.0 (https://bitbucket.org/mlindauer/aclib2/src/master/)

## Install Gurobi
```
wget https://packages.gurobi.com/9.0/gurobi9.0.1_linux64.tar.gz
tar xvfz gurobi9.0.1_linux64.tar.gz
cd gurobi901/linux64
python3 setup.py install
```
Next, add the following to the ```.bashrc``` file:
```
export GUROBI_HOME="your/path/here/gurobi901/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
```

To run Gurobi you need to obtain a (free academic) license.

## Install Venus requirements

```
pip3 install -r venus-requirements.txt
```

## General setup

Once you have downloaded the aclib library and installed its requirements, you should

1. place the contents from ```instances``` in the similary named folder in aclib and adapt the instance paths in training and test set files;
2. place the contents from ```target_algorithms``` in the similary named folder in aclib.

## Configure Venus

To configure Venus, we use a wrapper that directly accesses the Gurobi parameters through the Venus interface. To start the configuration procedure, you should

1.  adapt the network path in ```target_algorithms/venus/run_venus.py```;
2.  adapt the binary path in ```target_algorithms/venus/wrapper.py```;
3.  specify a scenario for SMAC (see ```example_scenario_venus.txt``` for reference) and place it in the respective aclib folder;
4.  run Hydra with the following commands: 

```/your/path/here/hydra-1.1-development-cae8151/bin/hydra --num-iterations 4 --num-smac-runs 2 --num-configs-per-iter 1 --rungroup Hydra_Venus --num-run 1 --smacOptions /your/path/here/aclib2/scenarios/.../your-scenario-file.txt --smac-execution-options /your/path/here/hydra-1.1-development-cae8151/smac-execution-options-local.txt```

The resulting Gurobi configuration(s) can be set directly in Venus to perform evaluation on the full MNIST dataset (we provide a full version of the MNIST dataset in ```.pkl``` format, see ```MNIST-full```).

## Configure MIPVerify

To configure MIPVerify, we extracted MIP problem formulations and configure Gurobi independently. To start the configuration procedure, you should

1. adapt the binary path in ```target_algorithms/gurobi902/wrapper.py```;
2. specify a scenario for SMAC (see ```example_scenario_gurobi.txt``` for reference) and place it in the respective aclib folder
3. run Hydra as outlined above.

The resulting Gurobi configuration(s) can be set directly in MIPVerify to perform evaluation on the full MNIST dataset. MIPVerify contains functions that can automatically import the network from [3] (```get_example_network_params("MNIST.RSL18a_linf0.1_authors")```) and load the full MNIST dataset (```read_datasets("MNIST")```).

# References
[1] Botoeva E, Kouvaros P, Kronqvist J, Lomuscio A, Misener R (2020) Efficient Verification of ReLU-based Neural Networks via Dependency Analysis. In: Proceedings of The Thirty-Fourth AAAI Conference on Artificial Intelligence (AAAI20), pp 3291–3299

[2] Tjeng V, Xiao K, Tedrake R (2019) Evaluating Robustness of Neural Networks with Mixed Integer Programming. In: Proceedings of the 7th International Conference on Learning Representations (ICLR 2019)

[3] Raghunathan A, Steinhardt J, Liang P (2018) Certified Defenses against Adversarial Examples. arXiv preprint arXiv:180109344
