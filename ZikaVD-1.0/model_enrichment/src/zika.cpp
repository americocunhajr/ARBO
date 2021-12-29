/*-------------------------------------------------------------------
 * ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
 *-----------------------------------------------------------------*/

/*-------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This is the inverse problem to calibrate the model discrepancy parameters
 *-----------------------------------------------------------------*/

#include <compute.h>

int main(int argc, char* argv[])
{
  // Initialize QUESO environment
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);

  // Call application
  computeParams(*env);

  // Finalize QUESO environment
  delete env;
  MPI_Finalize();

  return 0;
}
