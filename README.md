## Arbovirus Modeling and Uncertainty Quantification Toolbox

**ARBO: Arbovirus Modeling and Uncertainty Quantification Toolbox** is a comprehensive Matlab/C++ package designed for the simulation and analysis of arbovirus nonlinear dynamics. Developed with an educational approach, **ARBO** is intuitive and user-friendly, making it accessible for researchers and students alike.

### Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Usage](#usage)
- [Documentation](#documentation)
- [Reproducibility](#reproducibility)
- [Authors](#authors)
- [Citing ARBO](#citing-arbo)
- [License](#license)
- [Institutional support](#institutional-support)
- [Funding](#funding)

### Overview
**ARBO** was developed to simulate the nonlinear dynamics of an epidemic model to describe the Zika Virus outbreak in Brazil. It includes modules for solving initial value problems, calibration problems, model enrichment, and uncertainty quantification.

This code was developed to simulate the nonlinear dynamics of a epidemic model to describe Zika Virus outbreak in Brazil. It also solves an inverse problem to calibrate the underlying dynamic model parameters using real data as reference. These results are reported in the following paper:
- **E. Dantas, M. Tosin, A. Cunha Jr**, *Calibration of a SEIR–SEI epidemic model to describe the Zika virus outbreak in Brazil*, Applied Mathematics and Computation, vol. 338, pp. 249-259, 2018 <a href="https://doi.org/10.1016/j.amc.2018.06.024" target="_blank">DOI</a>

A third module includes a model enrichment approach, that uses discrepancy operator calibrated with data to compensate epidemic uncertainties in the epidemic model structure. The framework and some results are reported in:
- **R. E. Morrison, A. Cunha Jr**, *Embedded model discrepancy: A case study of Zika modeling*, Chaos, v. 30, pp. 051103, 2020 <a href="https://doi.org/10.1063/5.0005204" target="_blank">DOI</a>

The code also includes an uncertainty quantification module, that uses a probabilistic model to deal with the model parameters uncertainties. This framework and some results of the stochastic simulations are reported in:
- **E. Dantas, M. Tosin, A. Cunha Jr**, *An uncertainty quantification framework for a Zika virus epidemic model*, Journal of Computational Interdisciplinary Sciences, v. 10, pp. 91-96, 2019 <a href="http://dx.doi.org/10.6062/jcis.2019.10.02.0163" target="_blank">DOI</a>

### Features
- **Data Sets:** Preprocessed data for simulations
- **Initial Value Problem:** Solves forward problems using ODE solvers
- **Calibration Problem:** Calibrates model parameters using real data
- **Model Enrichment:** Enhances models with discrepancy operators
- **Uncertainty Quantification:** Propagates uncertainties via Monte Carlo method

### Usage
To get started with **ARBO**, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/americocunhajr/ARBO.git
   ```
2. Navigate to the code directory:
   ```bash
   cd ARBO/ARBO-1.0
   ```

The Matlab main routines and functions of the code are described below:
- main_SEIR_SEI_IVP_XX.m - Defines parameters and IC for the forward problem; solves the IVP with ode45; plots the time series
- rhs_SEIR_SEI.m - System of diferential equations for the IVP
- main_SEIR_SEI_TRR_XX.m - Sets up scenarios for the inverse problem and options for the TRR solver; plots time series
- ObjFun_SEIR_SEI.m - Sets up the objective function used in the main file for the inverse problem
- main_SEIR_SEI_MC_XX.m - Compute the propagation of uncertainties via Monte Carlo method

A description C++ program can be seen inside model_enrichment directory, where you can find a README file with instructions.

### Documentation
The routines in **ARBO** are well-commented to explain their functionality. Each routine includes a description of its purpose, as well as inputs and outputs. 

### Reproducibility
Simulations done with **ARBO** are fully reproducible, as can be seen on this <a href="https://codeocean.com/capsule/8169007/tree/v4" target="_blank">CodeOcean capsule</a>

### Authors
- Michel Tosin
- Eber Dantas
- Americo Cunha
- Rebecca E. Morrison

### Citing ARBO
If you use **ARBO** in your research, please cite the following publications:
- *M. Tosin, E. Dantas, A. Cunha Jr, R. E. Morrison, ARBO: Arbovirus modeling and uncertainty quantification toolbox, Software Impacts, vol. 12, pp. 100252, 2022 https://doi.org/10.1016/j.simpa.2022.100252*
- *E. Dantas, M. Tosin, A. Cunha Jr, Calibration of a SEIR–SEI epidemic model to describe the Zika virus outbreak in Brazil,  Applied Mathematics and Computation, vol. 338, pp. 249-259, 2018 https://doi.org/10.1016/j.amc.2018.06.024*
- *E. Dantas, M. Tosin, A. Cunha Jr, An uncertainty quantification framework for a Zika virus epidemic model, Journal of Computational Interdisciplinary Sciences, v. 10, pp. 91-96, 2019 http://dx.doi.org/10.6062/jcis.2019.10.02.0163*
- *R. E. Morrison, A. Cunha Jr, Embedded model discrepancy: A case study of Zika modeling, Chaos, v. 30, pp. 051103, 2020 https://doi.org/10.1063/5.0005204*

```
@article{Tosin2022ARBO,
   author  = {M. Tosin and E. Dantas and A. {Cunha~Jr} and R. E. Morrison},
   title   = {{ARBO: A}rbovirus modeling and uncertainty quantification toolbox},
   journal = {Software Impacts},
   year    = {2022},
   volume  = {12},
   pages   = {100252},
   doi     = {10.1016/j.simpa.2022.100252},
}
```

```
@article{Dantas2018p249,
   author  = {E. Dantas and M. Tosin and A. {Cunha~Jr}},
   title   = {Calibration of a {SEIR–SEI} epidemic model to describe the {Z}ika virus outbreak in {B}razil},
   journal = {Applied Mathematics and Computation},
   year    = {2018},
   volume  = {338},
   pages   = {249-259},
   doi     = {10.1016/j.amc.2018.06.024},
}
```

```
@article{Dantas2019p91,
   author  = {E. Dantas and M. Tosin and A. {Cunha~Jr}},
   title   = {An uncertainty quantification framework for a {Z}ika virus epidemic model},
   journal = {Journal of Computational Interdisciplinary Sciences},
   year    = {2019},
   volume  = {10},
   pages   = {91-96},
   doi     = {10.6062/jcis.2019.10.02.0163},
}
```

```
@article{Morrison2020p051103,
   author  = {R. E. Morrison and A. {Cunha~Jr}},
   title   = {Embedded model discrepancy: {A} case study of {Z}ika modeling},
   journal = {Chaos},
   year    = {2020},
   volume  = {30},
   pages   = {051103},
   doi     = {10.1063/5.0005204},
}
```

### License

**ARBO** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.

<img src="logo/mit_license_red.png" width="10%"> 

### Institutional support

<img src="logo/logo_uerj_color.jpeg" width="10%"> &nbsp; &nbsp; <img src="logo/logo_uc-boulder_color.png" width="35%">

### Funding

<img src="logo/faperj.jpg" width="20%"> &nbsp; &nbsp; <img src="logo/cnpq.png" width="20%"> &nbsp; &nbsp; <img src="logo/capes.png" width="10%"> &nbsp; &nbsp; <img src="logo/stem2d_logo.png" width="30%">
