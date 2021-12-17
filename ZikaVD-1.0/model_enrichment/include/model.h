 /*------------------------------------------------------------------
 * Brief description of this file: 
 *
 * This is the header file from src/model.cpp 
 *-----------------------------------------------------------------*/

#ifndef __ZIKA_MODEL_H__
#define __ZIKA_MODEL_H__

#include "dynamics_info.h"
#include <vector>

void
zikaComputeModel(
  std::vector<double>&  initialValues,
  std::vector<double>&  timePoints,
  dynamics_info*        p_dyn,
  std::vector<double>&  returnValues);

#endif
