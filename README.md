# QuantumFrankWolfe
Official implementation of Q-FW: A Hybrid Classical-Quantum Frank-Wolfe for Quadratic Binary Optimization (ECCV 2022)

Created by <a href="https://alpyurtsever.github.io/" target="_blank">Alp Yurtsever</a>, <a href="https://tolgabirdal.github.io/" target="_blank">Tolga Birdal</a>, <a href="https://people.mpi-inf.mpg.de/~golyanik/" target="_blank">Vladislav Golyanik</a>.

This repository contains the implementation of our [ECCV 2022 paper *Quantum Frank Wolfe*](https://arxiv.org/abs/2203.12633), a Hybrid Classical-Quantum Frank-Wolfe for Quadratic Binary Optimization. In particular, we release code for the main algorithm and provide an example on solving graph matching.

## Abstract
We present a hybrid classical-quantum framework based on the Frank-Wolfe algorithm, Q-FW, for solving quadratic, linearly-constrained, binary optimization problems on quantum annealers (QA). The computational premise of quantum computers has cultivated the re-design of various existing vision problems into quantum-friendly forms. Experimental QA realizations can solve a particular non-convex problem known as the quadratic unconstrained binary optimization (QUBO). Yet a naive-QUBO cannot take into account the restrictions on the parameters. To introduce additional structure in the parameter space, researchers have crafted ad-hoc solutions incorporating (linear) constraints in the form of regularizers. However, this comes at the expense of a hyper-parameter, balancing the impact of regularization. To date, a true constrained solver of quadratic binary optimization (QBO) problems has lacked. Q-FW first reformulates constrained-QBO as a copositive program (CP), then employs Frank-Wolfe iterations to solve CP while satisfying linear (in)equality constraints. This procedure unrolls the original constrained-QBO into a set of unconstrained QUBOs all of which are solved, in a sequel, on a QA. We use D-Wave Advantage QA to conduct synthetic and real experiments on two important computer vision problems, graph matching and permutation synchronization, which demonstrate that our approach is effective in alleviating the need for an explicit regularization coefficient.

## Contents

*gen_synth_multi_graph_2d*: generates synthetic data

*visualize_synth_match_data*: visualizes the data genrated in gen_synth_multi_graph_2d

*test_synth_data*: initial tests on the synthetic data

License: TBA. 
