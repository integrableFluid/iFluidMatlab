# iFluid: A MATLAB Library for Thermodynamic Bethe Ansatz and Generalized Hydrodynamics

`iFluid` is a MATLAB open-source library designed to solve equations related to the **Thermodynamic Bethe Ansatz (TBA)** and **Generalized Hydrodynamics (GHD)**. This framework is built to be highly flexible and extensible, making it easy for users to study integrable and near-integrable systems.

## Features

### 🛠 **Flexible, Extensible Framework**
- The core TBA and GHD equations are implemented in the **abstract class** `iFluidCore`, providing a universal backbone.
- Users can extend this class to implement **custom models** by defining model-specific properties such as scattering kernels.
- Pre-built support for popular models, including Hard rods, Lieb-Liniger, sine-Gordon, and XXZ spin chain.

### 🧮 **Custom `fluidcell` Data Structure**
- A specialized wrapper for MATLAB N-D arrays with overloaded operators for efficient kernel operations, including:
  - **Matrix multiplication (`*`)**: Simplifies convolutions over rapidity and quasi-particle types.
  - **Matrix "division" (`\`)**: Solves linear systems (e.g., for dressing operations) by inverting kernel.
  - **Example:** The discrete (implicit) dressing equation

    $$q_a^{\mathrm{dr}} (\theta_i) = \sum_{b,j} \left[I_{a,b} (\theta_i - \theta_j)  q_b (\theta_j) - d\theta_j T_{a,b} (\theta_i - \theta_j)  \vartheta_b(\theta_j)  q_b^{\mathrm{dr}}(\theta_j) \right] $$

    Is solved using the `fluidcell` operators as 
     ```matlab
     q_dr = (I - T.*transpose(dt.*v) / q;
     ```

### 🔍 **Advanced Solvers**
- A modular implementation of the **GHD partial differential equation solver**:
  - Abstract solver design allows users to implement custom solvers.
  - Pre-included solvers are based on **Backwards Semi-Lagrangian (BSL)** methods.
  - Supports higher-order integration (Runge-Kutta up to 4th order) for solving characteristic equations.

### ✨ **Beyond Euler-scale**
- Built-in functionality for:
  - **Position and time dependent couplings** (and associated effective accelerations).  
  - Higher-order GHD with **diffusion**.
  - Collision integrals for systems with **non-integrable perturbations**.

### 🚀 **GPU Support**
- Accelerate calculations using MATLAB's built-in **GPU arrays**.
- Particularly beneficial for simulations with large rapidity grids.

### 🌊 **Zero-Temperature & Whitham Theory**
- A dedicated sub-library for **zero-temperature/ground state calculations**:
  - State parameterization using **Fermi contours**, i.e. water-bag simulations.
  - Includes tools for **Whitham modulation theory**, enabling the study of **dispersive shock waves**.

## Installation

Clone this repository to your local machine:
```bash
git clone https://github.com/your-username/iFluid.git
```
Then add the library to your MATLAB path:
```matlab
addpath('path/to/iFluid');
```

## Citation
When publishing calculations performed with iFluid, please cite the [iFluid paper](https://arxiv.org/abs/2001.02547). 

Details and benchmarks of the Backward Semi-Lagragian schemes are documented [in this paper](https://arxiv.org/abs/2212.12349).

## Documentation
Documentation for adding new models and solvers can be found [here!](https://integrablefluid.github.io/iFluidDocumentation/)

For ussage of iFluid, see example scripts located in the \Examples folder.

## References
iFluid uses the following scripts for utility purposes. They can be found in the \utils folder.

#### Gauss-Legendre Quadrature (legzo)

[https://github.com/Pazus/Legendre-Gauss-Quadrature](https://github.com/Pazus/Legendre-Gauss-Quadrature)

#### ConsoleProgressBar
Evgeny Pr (2019). ConsoleProgressBar [https://www.mathworks.com/matlabcentral/fileexchange/30297-consoleprogressbar](https://www.mathworks.com/matlabcentral/fileexchange/30297-consoleprogressbar), MATLAB Central File Exchange. 
