#include "dynamics_info.h"

//Constructor
dynamics_info::dynamics_info(
    const unsigned int & n_s,
    const unsigned int & n_times,
    const unsigned int & inad_type,
    const unsigned int & params_factor,
    std::vector<double> & deltas)
:
  N_s(n_s),
  N_times(n_times),
  Inad_type(inad_type),
  Params_factor(params_factor),
  Deltas(deltas)
{
}

//Destructor
dynamics_info::~dynamics_info()
{
}
