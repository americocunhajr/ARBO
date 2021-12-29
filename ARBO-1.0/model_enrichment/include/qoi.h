/*-------------------------------------------------------------------
 * ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
 *-----------------------------------------------------------------*/

/*-------------------------------------------------------------------
 * Brief description of this file: 
 *
 * This is the header file from src/qoi.cpp. 
 *-----------------------------------------------------------------*/

#ifndef __ZIKA_QOI_H__
#define __ZIKA_QOI_H__

#include "dynamics_info.h"
#include <queso/GslMatrix.h>
#include <queso/DistArray.h>

struct qoiRoutine_Data
{
  qoiRoutine_Data(
      const QUESO::BaseEnvironment& env,
      const std::vector<double> & times,
      std::vector<double> & ics,
      dynamics_info * dynInfo);
 ~qoiRoutine_Data();

  const QUESO::BaseEnvironment* m_env;
  const std::vector<double> & m_times;
  std::vector<double> & m_ics;
  dynamics_info       * m_dynMain;
};

void qoiRoutine(
  const QUESO::GslVector&                     paramValues,
  const QUESO::GslVector*                     paramDirection,
  const void*                                 functionDataPtr,
        QUESO::GslVector&                     qoiValues,
        QUESO::DistArray<QUESO::GslVector* >* gradVectors,
        QUESO::DistArray<QUESO::GslMatrix* >* hessianMatrices,
        QUESO::DistArray<QUESO::GslVector* >* hessianEffects);

#endif
