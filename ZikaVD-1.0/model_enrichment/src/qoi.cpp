/*-------------------------------------------------------------------
 * ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
 *-----------------------------------------------------------------*/

/*-------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This file contains the code for the user defined qoi routine.
 *-----------------------------------------------------------------*/

#include "qoi.h"
#include "model.h"
#include "dynamics_info.h"
#include <cmath>
//------------------------------------------------------
/// The actual (user-defined) qoi routine
//------------------------------------------------------
// Constructor
qoiRoutine_Data::qoiRoutine_Data(
    const QUESO::BaseEnvironment& env,
    const std::vector<double> & times,
    std::vector<double> & ics,
    dynamics_info * dynInfo)
: m_env(&env),
  m_times(times),
  m_ics(ics),
  m_dynMain(dynInfo)
{
}

// Destructor
qoiRoutine_Data::~qoiRoutine_Data()
{
}

void
qoiRoutine(
  const QUESO::GslVector&                    paramValues,
  const QUESO::GslVector*                    paramDirection,
  const void*                                functionDataPtr,
        QUESO::GslVector&                    qoiValues,
        QUESO::DistArray<QUESO::GslVector*>* gradVectors,
        QUESO::DistArray<QUESO::GslMatrix*>* hessianMatrices,
        QUESO::DistArray<QUESO::GslVector*>* hessianEffects)
{
  const QUESO::BaseEnvironment& env = paramValues.env();

  if (paramDirection && 
      gradVectors    &&
      hessianEffects &&
      hessianMatrices) {
    // Logic just to avoid warnings from INTEL compiler
  }

  // get data from qoi data structure
  const std::vector<double>&  times
    = ((qoiRoutine_Data *) functionDataPtr)->m_times;
  std::vector<double>&  ics
    = ((qoiRoutine_Data *) functionDataPtr)->m_ics;
  dynamics_info *       dyn
    = ((qoiRoutine_Data *) functionDataPtr)->m_dynMain;

  const unsigned int n_s = dyn->N_s;          //the number of species included in the model
  const unsigned int n_times = dyn->N_times;  //the number of time points in every time series of data
  const unsigned int pf = dyn->Params_factor;    //the number of initial conditions
  const unsigned int n_params = pf * n_s;      //the number of parameters to be calibrated
  const unsigned int dim = n_s + 1;

  //set up lambda vector for loop, right now just one
  /* std::vector<double> phiPoints(n_phis,0.); */
  /* for (unsigned int j = 0; j < n_phis; j++){ */
  /*   phiPoints[j] = phis[n_times * j]; */
  /* } */

  std::vector<double> timePoints(n_times,0.);
  for (unsigned int i = 0; i < n_times; i++){
    timePoints[i] = times[i];
  }

  //return time points of all state variables + C
  std::vector<double> returnValues(n_times * dim, 0.);

  //std::cout << "hello in qoi----------------------------------\n";
  /* for (unsigned int i = 0; i < n_params; i++){  dyn->Deltas[i] = -std::exp(paramValues[i]); } */
  for (unsigned int i = 0; i < n_params; i++){  
      /* dyn->Deltas[i] = -std::abs(paramValues[i]); */
      dyn->Deltas[i] = paramValues[i];
  }
  /* for (unsigned int i = 0; i < n_params; i++){  dyn->Deltas[i] = 0.; } */

  try{
    zikaComputeModel(ics,timePoints,dyn,returnValues);
    //std::cout<< "qoi: ret val = " <<  returnValues[7 * 9 + 6] << std::endl; 
    for (unsigned int j = 0; j < returnValues.size(); j++){
      /* std::cout << "i = " << i << " and j = " << j << "\n"; */
      qoiValues[j] = returnValues[j];
//      std::cout << "returnValues[j] = " << returnValues[j] << "\n\n";
  //    std::cout << "j = " << j << "\n\n";
  }
  }catch(int exception)
  {
    for (unsigned int j = 0; j < returnValues.size(); j++){
      /* std::cout << "i = " << i << " and j = " << j << "\n"; */
        //qoiValues same as returnValues here
      qoiValues[j] = returnValues[j];}
    }

  return;
}
