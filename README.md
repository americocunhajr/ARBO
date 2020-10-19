<img src="logo/ZikaVD.png" width="40%">

**ZikaVD - Zika Virus Dynamics** is an easy to run Matlab code to simulate the nonlinear dynamics of the Zika virus. The implementation follows an educational style, to make its use very intuitive. 

This package includes the following modules:
- Initial Value Problem
- Calibration Problem
- Uncertainty Quantification

The most important routines and functions of the code are described below:
- main_SEIR_SEI.m - Defines parameters and IC for the forward problem; solves the IVP with ode45; plots the time series
- rhs_SEIR_SEI.m - System of diferential equations for the IVP
- TRR_main_SEIR_SEI.m - Sets up scenarios for the inverse problem and options for the TRR solver; plots time series
- TRR_rhs_SEIR_SEI.m - Adjusts the system of ODE according to the scenario chosen in the main file
- TRR_FunctionOutput_SEIR_SEI.m - Sets up the objective function used in the main file for the inverse problem
- main_SEIR_SEI_MonteCarlo.m - Compute the propagation of uncertainties via Monte Carlo method
- post_SEIR_SEI_MonteCarlo.m - Post processing for Monte Carlo simulation

## Software History

This code was developed to simulate the nonlinear dynamics of a epidemic model to describe Zika Virus outbreak in Brazil. It also solves an inverse problem to calibrate the underlying dynamic model parameters using real data as reference. These results are reported in the following paper:
- *E. Dantas, M. Tosin, A. Cunha Jr, Calibration of a SEIR–SEI epidemic model to describe the Zika virus outbreak in Brazil,  Applied Mathematics and Computation, vol. 338, pp. 249-259, 2018*
https://doi.org/10.1016/j.amc.2018.06.024

The code also includes an uncertainty quantification module, that uses a probabilistic model to deal with the model parameters uncertainties. This framework and some results of the stochastic simulations are reported in:
- *E. Dantas, M. Tosin, A. Cunha Jr, An uncertainty quantification framework for a Zika virus epidemic model, Journal of Computational Interdisciplinary Sciences, v. 10, pp. 91, 2019*
http://dx.doi.org/10.6062/jcis.2019.10.02.0163

## Authors
- Eber Dantas
- Michel Tosin
- Americo Cunha

## Citing ZikaVD

We kindly ask users to cite the following references in any publications reporting work done with **ZikaVD**:
- *E. Dantas, M. Tosin, A. Cunha Jr, Calibration of a SEIR–SEI epidemic model to describe the Zika virus outbreak in Brazil,  Applied Mathematics and Computation, vol. 338, pp. 249-259, 2018*
https://doi.org/10.1016/j.amc.2018.06.024
- *E. Dantas, M. Tosin, A. Cunha Jr, An uncertainty quantification framework for a Zika virus epidemic model, Journal of Computational Interdisciplinary Sciences, v. 10, pp. 91-96, 2019*
http://dx.doi.org/10.6062/jcis.2019.10.02.0163

```
@article{Dantas2018p249,
  author = {E. Dantas and M. Tosin and A. {Cunha Jr}},
  title = {Calibration of a {SEIR–SEI} epidemic model to describe the {Z}ika virus outbreak in {B}razil},
  journal = {Applied Mathematics and Computation},
  year = {2018},
  volume = {338},
  pages = {249-259},
  doi = {https://doi.org/10.1016/j.amc.2018.06.024},
}

@article{Dantas2019p91,
  author = {E. Dantas and M. Tosin and A. {Cunha Jr}},
  title = {An uncertainty quantification framework for a {Z}ika virus epidemic model},
  journal = {Journal of Computational Interdisciplinary Sciences},
  year = {2019},
  volume = {10},
  pages = {91-96},
  doi = {http://dx.doi.org/10.6062/jcis.2019.10.02.0163},
}
```

## License

**ZikaVD** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.
