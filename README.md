<img src="logo/ARBO.png" width="40%">

**ARBO: Arbovirus Modeling and Uncertainty Quantification Toolbox** is a package for simulation and analysis of arbovirus nonlinear dynamics. The implementation follows an educational style, to make its use very intuitive. 

This package includes the following modules:
- Data Sets
- Initial Value Problem
- Calibration Problem
- Model Enrichment
- Uncertainty Quantification

The Matlab main routines and functions of the code are described below:
- main_SEIR_SEI_IVP_XX.m - Defines parameters and IC for the forward problem; solves the IVP with ode45; plots the time series
- rhs_SEIR_SEI.m - System of diferential equations for the IVP
- main_SEIR_SEI_TRR_XX.m - Sets up scenarios for the inverse problem and options for the TRR solver; plots time series
- ObjFun_SEIR_SEI.m - Sets up the objective function used in the main file for the inverse problem
- main_SEIR_SEI_MC_XX.m - Compute the propagation of uncertainties via Monte Carlo method

A description C++ program can be seen inside model_enrichment directory, where you can find a README file with instructions.

## Software History

This code was developed to simulate the nonlinear dynamics of a epidemic model to describe Zika Virus outbreak in Brazil. It also solves an inverse problem to calibrate the underlying dynamic model parameters using real data as reference. These results are reported in the following paper:
- *E. Dantas, M. Tosin, A. Cunha Jr, Calibration of a SEIR–SEI epidemic model to describe the Zika virus outbreak in Brazil,  Applied Mathematics and Computation, vol. 338, pp. 249-259, 2018 https://doi.org/10.1016/j.amc.2018.06.024*

A third module includes a model enrichment approach, that uses discrepancy operator calibrated with data to compensate epidemic uncertainties in the epidemic model structure. The framework and some results are reported in:
- *R. E. Morrison, A. Cunha Jr, Embedded model discrepancy: A case study of Zika modeling, Chaos, v. 30, pp. 051103, 2020 https://doi.org/10.1063/5.0005204*

The code also includes an uncertainty quantification module, that uses a probabilistic model to deal with the model parameters uncertainties. This framework and some results of the stochastic simulations are reported in:
- *E. Dantas, M. Tosin, A. Cunha Jr, An uncertainty quantification framework for a Zika virus epidemic model, Journal of Computational Interdisciplinary Sciences, v. 10, pp. 91-96, 2019 http://dx.doi.org/10.6062/jcis.2019.10.02.0163*

## Reproducibility

Simulations done with **ARBO** are fully reproducible, as can be seen on this <a href="https://codeocean.com/capsule/2674283/tree" target="_blank">CodeOcean capsule</a>

## Authors
- Michel Tosin
- Eber Dantas
- Americo Cunha
- Rebecca E. Morrison

## Citing ARBO

We kindly ask users to cite the following references in any publications reporting work done with **ARBO**:
- *E. Dantas, M. Tosin, A. Cunha Jr, Calibration of a SEIR–SEI epidemic model to describe the Zika virus outbreak in Brazil,  Applied Mathematics and Computation, vol. 338, pp. 249-259, 2018 https://doi.org/10.1016/j.amc.2018.06.024*
- *E. Dantas, M. Tosin, A. Cunha Jr, An uncertainty quantification framework for a Zika virus epidemic model, Journal of Computational Interdisciplinary Sciences, v. 10, pp. 91-96, 2019 http://dx.doi.org/10.6062/jcis.2019.10.02.0163*
- *R. E. Morrison, A. Cunha Jr, Embedded model discrepancy: A case study of Zika modeling, Chaos, v. 30, pp. 051103, 2020 https://doi.org/10.1063/5.0005204*

```
@article{Dantas2018p249,
   author  = {E. Dantas and M. Tosin and A. {Cunha~Jr}},
   title   = {Calibration of a {SEIR–SEI} epidemic model to describe the {Z}ika virus outbreak in {B}razil},
   journal = {Applied Mathematics and Computation},
   year    = {2018},
   volume  = {338},
   pages   = {249-259},
   doi     = {https://doi.org/10.1016/j.amc.2018.06.024},
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
   doi     = {http://dx.doi.org/10.6062/jcis.2019.10.02.0163},
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
   doi     = {https://doi.org/10.1063/5.0005204},
}
```

## License

**ARBO** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.

## Institutional support

<img src="logo/logo_uerj_color.jpeg" width="10%"> &nbsp; &nbsp; <img src="logo/logo_uc-boulder_color.png" width="25%">

## Funding

<img src="logo/faperj.jpg" width="20%"> &nbsp; &nbsp; <img src="logo/cnpq.png" width="20%"> &nbsp; &nbsp; <img src="logo/capes.png" width="10%"> &nbsp; &nbsp; <img src="logo/stem2d_logo.png" width="30%">

