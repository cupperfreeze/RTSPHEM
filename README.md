## Overview

**RTSPHEM** — **Reactive Transport Solver in Porous Homogenized Evolving Media** includes software and data sets developed in the porous media research conducted at the Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU), Germany, within the Research Training Group RTG 2339 funded by DFG. More precisely, the code provided is related to the following publications:

- [-1] S. Gärttner, F. Frank, F. Woller, A. Meier, N. Ray:
**Estimating relative diffusion from 3D micro-CT images using CNNs**, arXiv, 2022.
https://arxiv.org/abs/2208.03337
- [0] M. Kelm, S. Gärttner, C. Bringedal, B. Flemisch, P. Knabner, N. Ray:
**Comparison study of phase-field and level-set method for three-phase systems including two minerals**, Computational Geosciences 26, pp. 545–570, 2022.
https://link.springer.com/article/10.1007/s10596-022-10142-w
- [1] S. Gärttner, F. O. Alpak, A. Meier, N. Ray, F. Frank:
**Estimating permeability of 3D micro-CT images by physics-informed CNNs based on DNS**, Computational Geosciences, 2023.
https://link.springer.com/article/10.1007/s10596-022-10184-0
- [2] S. Gärttner, P. Frolkovič, P. Knabner, N. Ray: 
**Efficiency of micro-macro models for ractive two-mineral systems**, SIAM MMS 20(1), pp. 433-461, 2022.
https://doi.org/10.1137/20M1380648
- [3] P. Frolkovič, N. Gajdošová, S. Gärttner, N. Ray:
**Voronoi implicit interface method for geometry evolution of two minerals with applications in reactive porous media**, Proceedings Of The Conference Algoritmy, pp. 121–130, 2020.
http://www.iam.fmph.uniba.sk/amuc/ojs/index.php/algoritmy/article/view/1566
- [4] S. Gärttner, P. Frolkovič, P. Knabner, N. Ray:
**Efficiency and accuracy of micro-macro models for mineral dissolution**, Water Resources Research 56 (8), 2020.
https://doi.org/10.1029/2020WR027585
- [5] N. Ray, J. Oberlander, P. Frolkovič:
**Numerical investigation of a fully coupled micro-macro model for mineral dissolution and precipitation**, Computational Geosciences 23, pp. 1173–1192, 2019.
https://doi.org/10.1007/s10596-019-09876-x
- [6] F. Frank; P. Knabner:
**Convergence analysis of a BDF2/mixed finite element discretization of a Darcy–Nernst–Planck–Poisson system**, ESAIM: Mathematical Modelling and Numerical Analysis (ESAIM: M2AN), 51 (5), pp. 1883–1902, 2017.
https://dx.doi.org/10.1051/m2an/2017002
- [7] R. Schulz, N. Ray, F. Frank, H. Mahato, P. Knabner:
**Strong solvability up to clogging of an effective diffusion–precipitation model in an evolving porous medium**, European Journal of Applied Mathematics, 28 (2), pp. 179–207, 2017.
https://dx.doi.org/10.1017/S0956792516000164
- [8] N. Ray, T. van Noorden, F. Frank, P. Knabner:
**Multiscale modeling of colloid and fluid dynamics in porous media including an evolving microstructure**, Transport in Porous Media 95 (3), pp. 669–696, 2012.
https://dx.doi.org/10.1007/s11242-012-0068-z
- [9] F. Frank, N. Ray, P. Knabner:
**Numerical investigation of homogenized Stokes–Nernst–Planck–Poisson systems**, Computing and Visualization in Science, 14 (8), pp. 385–400, 2011.
https://dx.doi.org/10.1007/s00791-013-0189-0


Besides finite element transport and flow solvers on fixed and variable geometries, this package includes an implementation of first and second order upwind finite difference level-set solvers and the Voronoi implicit interface method [2] for normal advection. For examination of such geometries, finite element solvers for diffusion and permeability computation are included (cf. cell problems in [2]). Furthermore, pretrained convolutional neural networks for permeability prediction as well as the training data acquisition and training script are contained. Scripts related to the simulations performed in the above publications are enclosed. The authors do not guarantee the exactness of results obtained using this software.

## Software environment

Our software is written in Matlab and builds upon the PDE toolbox HyPHM by Florian Frank (FAU) (https://www1.am.uni-erlangen.de/HyPHM/).
Therefore, execution requires a recent version of Matlab. In order use the full capacity of our implementation, the following Matlab toolboxes need to be installed as well:
- Deep Learning Toolbox,
- Parallel Computing Toolbox,
- Partial Differential Equation Toolbox.

Furthermore, some functions used in this project are also present in C/C++ code or optimized to work with Matlab Compiler and need to be compiled before usage (mex):
- `/LevelSetTriplePoint/src/FMM_iteration.m` (accelerate VIIM steps),
- `/LevelSetTriplePoint/src/FMM_iteration_WithS.m` (accelerate VIIM steps),
- `/ReactiveTransport/ExternalRoutines/Substitution/backsubs.cpp` (accelerate `Stokes`/`Perm` solves).

If the compiled routines are available to Matlab, they are automatically used to improve performance.

This software is compatible with several different linear solvers, preconditioners and matrix assembly methods that are not native to Matlab. To use all implemented features, make sure the following packages are installed:

- (1) **ILU(k) Preconditioner**: Killian Miller, MATLAB Central File Exchange, 2020 https://www.mathworks.com/matlabcentral/fileexchange/48320-ilu-k-preconditioner
- (2) **ilupack**: Bollhöfer Matthias and Yousef Saad: Multilevel preconditioners constructed from inverse-based ILUs, SIAM Journal on Scientific Computing, 27, 2006. http://ilupack.tu-bs.de
- (3) **fsparse**: Stefan Engblom, Dimitar Lukarski: Fast Matlab compatible sparse assembly on multicore computers, Parallel Computing, 2016. https://github.com/stefanengblom/stenglib/blob/master/Fast/fsparse.m

## Content

### `Data/`

Contains files in `.mat` format including data for simulation scripts:

- `ChenRelease.mat`:           Simulation output corresponding to [2], microscale simulation (due to their large size, these files are available on request).
- `CNNPermGeneral.mat`:        Contains a CNN to predict (1,1) component of permeability tensor on general geometries.
- `CNNPermRectangle.mat`:      Contains a CNN to predict (1,1) component of permeability tensor on rectangular geometries.
- `LevelConv2eps.mat`:         Contains VIIM convergence test corresponding to [2].
- `PermDataSPE.mat`:           Permeability data taken from SPE10 project.
- `PermLookupRect.mat`:        Permeability values corresponding rectangles (0.02,0.04,...,1) x (0.02,0.04,...,1).
- `SPE10.mat`:                 Simulation results corresponding to micro-macro simulation in [2] (due to their large size, these files are available on request).
- `TrainCNNDataGeneral.mat`:   Training data for above CNN.
- `TrainCNNDataRectangle.mat`: Training data for above CNN.
- `Further Data`:              Further simulation Data regarding [4] are available here: https://data.4tu.nl/articles/dataset/Efficiency_and_Accuracy_of_Micro-Macro_Models_for_Mineral_Dissolution_-_Simulation_Data/12776063.		

### `LevelSetSolverOld/`

Discontinued version of a standard level-set solver; still needed for some applications. The folder `TestScripts` contains files for testing and visualization.

### `LevelSetTriplePoint/`

Level-set method.

- `src/`: Source code for standard level-set solver and VIIM for an arbitrary number of phases, cf. [3]. Includes Eikonal solver, interface velocity extender, interface length/position evaluation and wrappers for compact application in micro-macro simulations.
- Scripts and visualization for different test cases (T-shape, circle) and convergence studies.

### `ReactiveTransport/`

Reactive transport simulations and source code.

- `HyPHM/`: Adapted version of HyPHM by Florian Frank; changes in short:
  - massive performance optimization for the classes `Grid`, `Transport` and `Stokes` by vectorization (possibly using Matlab's automated multithreading capabilities),
  - new classes `StokesLevel`, `TransportLevel` for simulations on geometry defined by a level-set function,
  - new class `Perm` for highly optimized calculation of permeability tensors, builds upon Stokes, with additional P1-bubble/P1 discretization option (compared to more expensive P2/P1),
  - Robin boundary conditions available for class `TransportLevel`.
  - *Warning*: new routines only tested for target application. For general use, double checking with original HyPHM implementation is strongly recommended.
- `littleScripts/`: scripts for main simulation (evaluation, visualization, ...).
- `Machine Learning/`: scripts for training data acquisition and network training for permeability prediction. Usage of Cuda GPUs is recommended.
- `scripts/`: Simulation scripts:
	- `BenchPart1Newton.m` : Benchmark test, cf. [4] Chapter 5.3.
	- `Chen.m`, cf. [2] Chapter 4.2. (similar to `ChenNeumann.m` and `ChenRobin.m` which use different outer boundary conditions)
	- `Dolomite_CircleDetection.m`, cf. [4] Chapter 3.2.2.
	- `DolomiteMicroMacro.m`, cf. [4] Chapter 5.3.
	- `DolomiteMicroStokesNewton.m`, cf. [4] Chapter 5.3.
        - `LS_ODE_Simulation.m`: Reactive two-mineral simulation without spatial resolution, cf. [0] Section 4.3
        - `LS_PDE_Simulation.m`: Reactive two-mineral simulation with spatial resolution (diffusion-advection-reaction equation), cf. [0] Section 4.4 / 4.5
	- `RandomPorosityField.m`, cf. [4] Chapter 5.2.
	- `Two_Min_Triple.m`: Working example micro-macro model with 2 mineral phases.
	- `Two_Min_TripleDataSet.m`, cf. [2]  Chapter 4.3.
	- `Two_Min_TripleDataSet_CNN.m`, cf. [2]  Chapter 4.3.
- `src/`: Source code:
	- `CellProblems/`: 
		- Extended finite element solver for diffusion calculation on geometries given by a level-set (P1 discretization), simultaneously calculating porosity and fluid-solid interface surface.
		- FEM solver for permeability calculation on geometries given by a level-set function, including local mesh refinement along interface and Uzawa's algorithm for saddle point problems.
	- `CircleDetection/`: Adaptivity feature
		- check if solid geometry resembles a circle by principle component analysis and estimate radius,
		- use radius to interpolate hydrodynamic parameters from lookup table. 
	- `LevelSetSolver2ndOrder/`: Finite difference level-set solver using a first or second order upwind scheme on quadratic Cartesian grids. Optimized for narrow band level-set functions. (Computations restricted to non-isnan-nodes).
	- `MacroscaleNewton/`: Newton solver to handle nonlinearly coupled transport equations (implicit Euler).
		- Uses HyPHM assembly routines and Armijo step size control (cf. [5]).
		- Different solver options for linear systems: (defined by global variable `Solver`):
			- `StandardDirect`: Matlab's backslash operator (default)
			- `ilupack`: use ilupack solver (needs (2) installed)
			- `BlockPrec`: Take advantage of system's structure by specialized solver, cf. [2] Chapter 3.2. Block diagonal ilu(k) preconditioner with resorting/rescaling (needs (1), (2) installed).
			- `ilukIterative`: Matlab's bicgstab with ilu(k) preconditioner (needs (1) installed).
			- `StandardIterative`: Matlab's bicgstab with ILU preconditioner.

### `3DRelDiffCNN/`

Code related to the training of a CNN on 3D segmented binary image data for scalar relative diffusion prediction, cf. [-1].
The underlying data-set is derived from the voxel-data of the `3DPermCNN/` experiment and can be derived using /Miscellaneous/ApplyMorph.m
This project shares a lot of code with `3DPermCNN/`.

- `CNN/`: Contains script to train CNN and pretrained network itself. Usage of Cuda GPUs is recommended.
- `ComputedData/`: Computed diffusion values (DiffData) and sample characteristics (surface area, porosity).
- `TauFactor/`: Scripts for performing diffusion computation using TauFactor (doi: 10.1016/j.softx.2016.09.002)
- `Miscellaneous/`: Scripts for plotting, morphological model and heuristic law fitting

### `3DPermCNN/`

Code related to the training of a physics-informed CNN on 3D segmented binary image data for scalar permeability prediction, cf. [1].

- `CNN/`: Contains script to train physics-informed CNN and pretrained network itself. Usage of Cuda GPUs is recommended.
- `ComputedData/`: Computed permeability values (PermData...) and sample characteristics (SampleData...) for different types of sandstone.
- `Miscellaneous/`: Scripts for pore space manipulation (using level-set methods), generation of samples with rectangular channel, inner surface area computation and evaluation of CNN computational performance. 
- `Postprocessing/`: Scripts for CNN validation and plots; Performance statistics of forward solver.
- `Preprocessing/`: Sampling operations, removal of disconnected pore space, computation of maximum flow values (physics input).
- `RawConversion/`: Reading raw files into Matlab and vice versa.
- `VoxelData/`: Preprocessed subsample voxel images. (Will be uploaded separatly and linked as soon as paper [1] is finally published.)

## Acknowledgements
The authors like to thank Jens Oberlander whose implementations provided the basis for many routines as well as Peter Frolkovič for his expertise in level-set algorithms.

Please report any bugs or improvement suggestions to Stephan Gärttner, Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany, gaerttner@math.fau.de.

## How to cite
```bibtex
@Online{RTSPHEM,
  Title                    = {RTSPHEM -- Reactive Transport Solver in Porous Homogenized Evolving Media},
  Author                   = {Stephan G\"arttner and Florian Frank},
  Url                      = {https://github.com/cupperfreeze/RTSPHEM/},
  DOI                      = {10.5281/zenodo.5166669},
  Year                     = {2021},
  Note                     = {Accessed: XXX},
  Organization             = {Department Mathematik, Friedrich-Alexander-Universit\"at Erlangen-N\"urnberg}
}
```
