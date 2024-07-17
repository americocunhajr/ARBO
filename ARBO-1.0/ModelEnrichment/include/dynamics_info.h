/*-------------------------------------------------------------------
 * ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
 *-----------------------------------------------------------------*/

#ifndef __DYNAMICS_INFO_H__
#define __DYNAMICS_INFO_H__

//c++
#include <vector>

// define struct that holds all dyanamical system info, except params
struct dynamics_info { dynamics_info(
  const unsigned int & n_s,
  const unsigned int & n_times,
  const unsigned int & inad_type,
  const unsigned int & params_factor,
  std::vector<double> & deltas);
 ~dynamics_info();

  const unsigned int & N_s;
  const unsigned int & N_times;
  const unsigned int & Inad_type;
  const unsigned int & Params_factor;
  std::vector<double> & Deltas;
};
#endif
