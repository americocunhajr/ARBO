<img src="logo/ZikaVD_logo.png" width="40%">

**ZikaVD - Zika Virus Dynamics** is an easy to run Matlab code to simulate the nonlinear dynamics of the Zika virus. The implementation follows an educational style, to make its use very intuitive. This package includes the following modules:

- main_SEIR_SEI.m - Defines parameters and IC for the forward problem; solves the IVP with ode45; plots the time series.
- rhs_SEIR_SEI.m - System of diferential equations for the IVP.
- TRR_main_SEIR_SEI.m - Sets up scenarios for the inverse problem and options for the TRR solver; plots time series.
- TRR_rhs_SEIR_SEI.m - Adjusts the system of ODE according to the scenario chosen in the main file.
- TRR_FunctionOutput_SEIR_SEI.m - Sets up the objective function used in the main file for the inverse problem.


Further details about can be seen in:
- *E. Dantas, M. Tosin, A. Cunha Jr, Calibration of a SEIR–SEI epidemic model to describe the Zika virus outbreak in Brazil,  Applied Mathematics and Computation, vol. 338, pp. 249-259, 2018*
https://doi.org/10.1016/j.amc.2018.06.024

## Authors:
- Eber Dantas
- Michel Tosin
- Americo Cunha

## Citing ZikaVD:

We kindly ask users to cite the following reference in any publications reporting work done with **ZikaVD**:
- *E. Dantas, M. Tosin, A. Cunha Jr, Calibration of a SEIR–SEI epidemic model to describe the Zika virus outbreak in Brazil,  Applied Mathematics and Computation, vol. 338, pp. 249-259, 2018*
https://doi.org/10.1016/j.amc.2018.06.024


## License

**ZikaVD** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.
