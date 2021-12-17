 /*------------------------------------------------------------------
 * Brief description of this file: 
 *
 * This is the header file for src/likelihood.cpp. 
 *-----------------------------------------------------------------*/

#ifndef __ZIKA_LIKELIHOOD_H__
#define __ZIKA_LIKELIHOOD_H__

#include "dynamics_info.h"
#include <queso/GslMatrix.h>

struct likelihoodRoutine_Data // user defined class
{
  likelihoodRoutine_Data(
      const QUESO::BaseEnvironment& env,
      const std::vector<double> & times,
      std::vector<double> & ics,
      const std::vector<double> & csc,
      double & var,
      dynamics_info * dynInfo);
 ~likelihoodRoutine_Data();

  const QUESO::BaseEnvironment* m_env;
  const std::vector<double> & m_times;
  std::vector<double> & m_ics;
  const std::vector<double> & m_csc;
  double & m_var;
  dynamics_info       * m_dynMain;
};

double likelihoodRoutine( // user defined routine
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect);

#endif
