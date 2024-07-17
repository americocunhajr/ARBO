 /*-------------------------------------------------------------------
 * ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
 *-----------------------------------------------------------------*/

/*-------------------------------------------------------------------
 * This file contains the code for the user defined likelihood data 
 * class and the user defined likelihood routine.
 *-----------------------------------------------------------------*/

#include "likelihood.h"
#include "dynamics_info.h"
#include "model.h"
#include <cmath>
#include <stdio.h>
#include <fstream>

// Constructor
likelihoodRoutine_Data::likelihoodRoutine_Data(
    const QUESO::BaseEnvironment& env,
    const std::vector<double> & times, 
    std::vector<double> & ics,
    const std::vector<double> & csc,
    double & var,
    dynamics_info * dynInfo)
: m_env(&env),
  m_times(times),
  m_ics(ics),
  m_csc(csc),
  m_var(var),
  m_dynMain(dynInfo)
{
}

// Destructor
likelihoodRoutine_Data::~likelihoodRoutine_Data()
{
}

//------------------------------------------------------
// The user defined likelihood routine
//------------------------------------------------------

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);
    
  if (paramDirection && functionDataPtr && gradVector && hessianMatrix && hessianEffect) 
  {
    // Just to eliminate INTEL compiler warnings
  }
  
  env.subComm().Barrier(); 
  
  // Compute likelihood 
  // get data from likelihood data structure
  const std::vector<double>&  times
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_times;
  std::vector<double>& ics
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_ics;
  const std::vector<double>& csc 
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_csc;
  double & var
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_var;
  dynamics_info *        dyn
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_dynMain;

  const unsigned int n_s = dyn->N_s;          //the number of species included in the model
  const unsigned int n_times = dyn->N_times;  //the number of time points in every time series of data
  const unsigned int n_params = dyn->Params_factor * n_s;      //the number of parameters to be calibrated

  unsigned int dim = n_s + 1;
  //set up lambda vector for loop, right now just one
  /* std::vector<double> phiPoints(n_phis,0.); */

  /* for (unsigned int j = 0; j < n_phis; j++){ */
  /*   phiPoints[j] = phis[n_times * j]; */
  /* } */

  std::vector<double> timePoints(n_times,0.);
  for (unsigned int i = 0; i < n_times; i++){
    timePoints[i] = times[i];
    // std::cout << "timePoints[i] = " << timePoints[i];
  }

  //return all times of C (cumulative cases)
  std::vector<double> returnValues(n_times * dim, 0.);

  double misfitValue = 0.;
  double diff = 0.;

  /* for (unsigned int i = 0; i < n_params; i++){  dyn->Deltas[i] = -std::exp(paramValues[i]); } */
  for (unsigned int i = 0; i < n_params; i++){  
      /* dyn->Deltas[i] = -std::abs(paramValues[i]); */ 
      dyn->Deltas[i] = paramValues[i]; 
      /* std::cout << "dyn->Deltas[i] = " << dyn->Deltas[i]<< "\n"; */
  }
  /* for (unsigned int i = 0; i < n_params; i++){  dyn->Deltas[i] = 0.; } */

  try
     {
      zikaComputeModel(ics,timePoints,dyn,returnValues);
//      std::cout << "Finished compute model" << std::endl;
      for (unsigned int j = 0; j < n_times; j++){
          //only have data for Y[7]
        diff = (returnValues[dim * j + 7] - csc[j]);
//          std::cout<<"dim*j+ 7" << dim*j+7 << "\n";
  //        std::cout<<"likelihood: " << returnValues[dim*j+7] << "compared to " << csc[j] << "\n";
          misfitValue += diff * diff / var;
          /* std::cout << "misfit = " << misfitValue << "\n"; */
          /* count += 1; */
        }
     } catch( int exception )
     {
      misfitValue = 1000000;
   }

  /* std::cout << "likelihood count = " << count << "\n"; */
  /* std::cout << " the misfit is " << misfitValue << std::endl; */
  return (-0.5 * misfitValue);
}
