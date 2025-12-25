This repository contains experimental code for evaluating computationally expensive and offline optimization algorithms applied to airfoil shape optimization problems. The airfoil shape optimization framework is adapted from the work: [Optimization of an airfoil shape using Genetic algorithm](https://www.mathworks.com/matlabcentral/fileexchange/69095-optimization-of-an-airfoil-shape-using-genetic-algorithm-ga).

The algorithms implemented and compared in this study include:

**Expensive Optimization Algorithms:**
- **SHPSO**: Surrogate-assisted Hierarchical Particle Swarm Optimization  
  Yu H, Tan Y, Zeng J, et al. Surrogate-assisted hierarchical particle swarm optimization[J]. Information Sciences, 2018, 454: 59-72.
- **ESAO**: Evolutionary Sampling Assisted Optimization  
  Wang X, Wang G G, Song B, et al. A novel evolutionary sampling assisted optimization method for high-dimensional expensive problems[J]. IEEE Transactions on Evolutionary Computation, 2019, 23(5): 815-827.
- **DESO**: Data-driven Evolutionary Sampling Optimization  
  Huixiang Z, Wenyin G, Ling W. Data-driven evolutionary sampling optimization for expensive problems[J]. Journal of Systems Engineering and Electronics, 2021, 32(2): 318-330.
- **TS-DDEO**: Two-stage Data-driven Evolutionary Optimization  
  Zhen H, Gong W, Wang L, et al. Two-stage data-driven evolutionary optimization for high-dimensional expensive problems[J]. IEEE transactions on cybernetics, 2021, 53(4): 2368-2379.
- **ESA**: Evolutionary Sampling Agent  
  Zhen H, Gong W, Wang L. Evolutionary sampling agent for expensive problems[J]. IEEE Transactions on Evolutionary Computation, 2022, 27(3): 716-727.

**Offline Data-driven Evolutionary Optimization Algorithms:**
- **DDEA-SE**: Offline Data-driven EA with Selective Surrogate Ensembles  
  Wang H, Jin Y, Sun C, et al. Offline data-driven evolutionary optimization using selective surrogate ensembles[J]. IEEE Transactions on Evolutionary Computation, 2018, 23(2): 203-216.
- **DDEA-LDG**: Data-driven EA with Localized Data Generation  
  Li J Y, Zhan Z H, Wang C, et al. Boosting data-driven evolutionary algorithm with localized data generation[J]. IEEE Transactions on Evolutionary Computation, 2020, 24(5): 923-937.
- **MS-DDEO**: Offline Data-driven EO based on Model Selection  
  Zhen H, Gong W, Wang L. Offline data‐driven evolutionary optimization based on model selection[J]. Swarm and Evolutionary Computation, 2022, 71: 101080.
- **MSEA**: Offline Evolutionary Optimization with Problem-driven Model Pool & Weighted Model Selection  
  Zhen H, Xue B, Gong W, et al. Offline evolutionary optimization with problem-driven model pool design and weighted model selection indicator[J]. Swarm and Evolutionary Computation, 2025, 97: 102034.

This project serves as a testbed for benchmarking advanced optimization strategies in the context of aerodynamic shape design, where objective function evaluations are costly or only historical/offline data is available.