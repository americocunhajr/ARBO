 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This file contains the code for the forward model.
 *-----------------------------------------------------------------*/

#include "model.h"
/* #include "dynamics_info.h" */
#include <cmath>
#include <vector>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <assert.h>
#include "eigen3/Eigen/Dense"
/* //antioch */
/* #include <antioch/kinetics_evaluator.h> */
/* #include <antioch/cea_evaluator.h> */

#ifndef __EPS_ABS
#define __EPS_ABS 1e-8
#endif
#ifndef __EPS_REL
#define __EPS_REL 1e-8
#endif

//first define the function for the ODE solve
int zikaFunction( double t,
                 const double Y[],
                 double dYdt[],
                 void* params)
{
  //here, params is sending the function all the reaction info
  dynamics_info dyn = *(dynamics_info *) params;

  //reduced model parameters:
  double bh = 1/11.3; //\beta_h
  double ah = 1/5.9;  //\alpha_h
  double g = 1/7.9;   //\gamma
  double d = 1/11.;   //\delta
  double bv = 1/8.6;  //\beta_v
  double av = 1/9.1;  //\alpha_v
  double nv = 1.;     //N_v
  double nh = 206 * pow(10,6);

  const unsigned int n_s = dyn.N_s;
  const unsigned int inad_type = dyn.Inad_type;
  const unsigned int pf = dyn.Params_factor;
  const unsigned int n_deltas = pf * n_s;
  std::vector<double> delta = dyn.Deltas;

  //use pops to copy ``populations'' of state variables
  std::vector<double> pops(n_s + 1,0);
  for (unsigned int i = 0; i < n_s + 1; i++){
    pops[i] = Y[i];
    if(pops[i] <= 0){
      pops[i] = 0;
    }
  }

  //SEIR-SEI model
  dYdt[0] = -bh * pops[0] * pops[6] / nv;
  dYdt[1] = bh * pops[0] * pops[6] / nv - ah * pops[1];
  dYdt[2] = ah * pops[1]  - g * pops[2];
  dYdt[3] = g * pops[2];
  dYdt[4] = d * nv - bv * pops[4] * pops[2] / nh - d * pops[4];
  dYdt[5] = bv * pops[4] * pops[2] / nh - (av + d) * pops[5];
  dYdt[6] = av * pops[5] - d * pops[6];
  dYdt[7] = ah * pops[1];

//inadequacy formulation
  if ( inad_type == 0) {
    for (unsigned int i = 0; i < n_s; i++){
      dYdt[i] += 0.;
    }
  }
  else if ( inad_type == 1) {
    for (unsigned int i = 0; i < n_s; i++){
      dYdt[i] += delta[pf*i+0]*pops[i] + delta[pf*i+1]*std::abs(dYdt[i]);
      /* dYdt[i] += delta[pf*i+0] + delta[pf*i+1]*pops[i] + delta[pf*i+2]*std::abs(dYdt[i]); */
    }
  }
  else if ( inad_type == 2) {
    for (unsigned int i = 0; i < n_s; i++){
      dYdt[i] += delta[pf*i+0]*pops[i] + delta[pf*i+1]*std::abs(dYdt[i]) + 
          delta[pf*i+2]*std::pow(pops[i],2) + delta[pf*i+3]*std::pow(dYdt[i],2);
    }
  }
  else if ( inad_type == 3) {
    for (unsigned int i = 0; i < n_s; i++){
        for (unsigned int j = 0; j < n_s; j++){
            dYdt[i] += delta[pf*i + 2*j +0]*pops[j] + delta[pf*i + 2*j + 1]*std::abs(dYdt[j]);
            dYdt[i] += delta[pf*i + 2*j +2]*std::pow(pops[j],2) + delta[pf*i + 2*j + 3]*std::pow(dYdt[j],2);
        }
    }
  }

  if (Y[7] > 3.0e9) {
      //std::cout << "she's gonna blow!\n\n";
      //throw 20;
  }
  return GSL_SUCCESS;
}

//jacobian for ode solve---------------------------------------------
int zikaJacobian( double t, 
				const double Y[],
				double *dfdY,
				double dfdt[],
				void* params )
{
  return GSL_SUCCESS;
}

void zikaComputeModel(
  std::vector<double>&  initialValues,
  std::vector<double>&  timePoints,
  dynamics_info*        dyn,
  std::vector<double>&  returnValues)
{  
  // Compute model
  // GSL prep
  unsigned int dim = initialValues.size();
  unsigned int n_s = dim - 1;
  gsl_odeiv2_system sys = { zikaFunction, 
			   zikaJacobian, 
			   dim, dyn };
  
  double h = 1e-10;    //initial step-size
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rkf45,h,1e-8,1e-4);   
  // initialize values
  double Y[dim];
  for (unsigned int i = 0; i < dim; ++i){
    Y[i] = initialValues[i];
    returnValues[i] = initialValues[i];
  };
  double t = 7.0;
  double prevt;
  double sumY;
//  std::cout << "Starting Integration..." << std::endl;
  // GSL integration
  double finalTime;
  for (unsigned int i = 1; i < timePoints.size(); i++){
    finalTime = timePoints[i];
    while (t < finalTime)
      {
     //   prevt = t;
     //   sumY = 0;
        // t and y are updated and placed back into those variables
        int status = gsl_odeiv2_driver_apply( d, // necessary gsl vars
               &t,    // current time
               finalTime,    // maY time
               Y );   // current solution values (at time t)
//        std::cout<<"Time = "<<t<<"\nAfter integration Y values : \n";
//        for (unsigned int i = 0; i <= n_species; i++) std::cout<<Y[i]<<"\n";
        /* std::cout<<"t = "<<t<<"\n"; */
        /* std::cout<<"Y[0] = "<<Y[0]<<"\n"; */
        /* std::cout<<"Y[1] = "<<Y[1]<<"\n"; */
        //std::cout<<"T = "<<Y[n_species]<<"\n";
//        for( i = 0; i<7;i++) sumY+=Y[i];
//        std::cout<<"N = "<<sumY<<"\n";
        // check that the evolution was successful
        #ifdef UQ_FATAL_TEST_MACRO
          UQ_FATAL_TEST_MACRO( status != GSL_SUCCESS,
             0,
             "ZIKA",
             "The status of GSL integration != GSL_SUCCESS" );
        #else 
          if ( status != GSL_SUCCESS )
          {
     //       std::cout<< "h approx = " << t-prevt<<"\n\n";
            std::cout << "ERROR: status of GSL integration != GSL_SUCCESS" <<
              std::endl;
            assert( status == GSL_SUCCESS );
          }
        #endif
      }
    //  std::cout << " h is " << h << std::endl;
    // save results, right now return values are all the species of the reduced
    // model
    for (unsigned int j = 0; j < dim; j++){
      returnValues[dim*i +j] = Y[j];}
  }
  // std::cout << "C = " << Y[7] << std::endl;
  // std::cout << "O2 = " << Y[1] << std::endl;
  // deallocate memory
  gsl_odeiv2_driver_free( d );
}
