/*-------------------------------------------------------------------
 * This file is divided in two parts:
 * - the first one handles the statistical inverse problem (SIP) for estimating
 *   the discrepancy parameters '\delta_i', i=1,...,14
 * - the second part handles the statistical forward problem (SFP) for
 *   predicting the QoI
 *-----------------------------------------------------------------*/
/* #include <fstream> */
#include "compute.h"
#include "likelihood.h"
#include "qoi.h"
#include "dynamics_info.h"
//queso
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GenericVectorRV.h>
#include <queso/GaussianVectorRV.h>
#include <queso/LogNormalVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include <sys/time.h>
#include <cmath>
#include <vector>

void computeParams(const QUESO::FullEnvironment& env) {
  struct timeval timevalNow;
  
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "\nBeginning run of 'Zika plus discrepancy' example at "
              << ctime(&timevalNow.tv_sec)
              << "\n my fullRank = "         << env.fullRank()
              << "\n my subEnvironmentId = " << env.subId()
              << "\n my subRank = "          << env.subRank()
              << "\n my interRank = "        << env.inter0Rank()
               << std::endl << std::endl;
  }

  // Just examples of possible calls
  if ((env.subDisplayFile()       ) && 
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Beginning run of 'Zika plus discrepancy' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  env.fullComm().Barrier();
  env.subComm().Barrier();  // Just an example of a possible call
  
  //================================================================
  // Statistical inverse problem (SIP): find posterior PDF for \delta's
  //================================================================
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'SIP -> all parameters estimation' at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  //------------------------------------------------------
  // SIP Step 0 of 6: Read in the data
  //------------------------------------------------------
  unsigned int n_s;  //number of species in model
  double var = 25000000;            //variance in the data
  //type of inadequacy model: right now only one type, might include more later
  unsigned int inad_type = 1;


  //number of species in SEIR-SEI model is 7:
  //S_h, E_h, I_h, R_h, S_v, E_v, I_v
  n_s = 7;
  unsigned int dim = n_s + 1;
//  std::cout << "inadequacy type = " << inad_type << "\n\n";

  //inad_type:
  //0: nothing
  //1: 2 terms per state variable, linear in state and derivative
  //STILL TODO!!!! 
  //2: same as 1 but with hyperparameters for mean and variance
  //3: quadratic 
  unsigned int params_factor;
  if( inad_type == 0 ) { params_factor = 1;}
  if( inad_type == 1 ) { params_factor = 2;}
  if( inad_type == 2 ) { params_factor = 6;}
  if( inad_type == 3 ) { params_factor = 2 * n_s;}
  unsigned int n_delta = params_factor*n_s;         //the model discrepancy terms
  unsigned int n_params = n_delta;                  //no other parameters to calibrate
  unsigned int n_weeks = 52;

  //read in data points
  double tmpPhis; //number of initial conditions (right now only one)
  double tmpWeeks;
  double tmpx;
  int numLines = 0;

  FILE *dataFile;
  dataFile = fopen("./inputs/data.txt","r");

  std::vector<double> weeks(n_weeks, 0.);
  std::vector<double> new_cases(n_weeks, 0.);
  std::vector<double> initialValues(dim, 0.);
  // MUST SET BY HAND HERE, should change to input
  double rep_factor = 1.;    // no under-reporting
  //double rep_factor = 10./9; // 10% under-reporting
  //double rep_factor = 2.0;     // 50% under-reporting
  while (fscanf(dataFile,"%lf %lf ", &tmpWeeks, &tmpx) != EOF) {
    weeks[numLines]    = tmpWeeks;
    //add 10% for under-reporting
    new_cases[numLines]     = rep_factor * tmpx;
    numLines++;
  }
  //count cumulative sum of new cases
  std::vector<double> cum_sum_cases(n_weeks, 0.);
  cum_sum_cases[0] = new_cases[0];
  for (unsigned int i = 1; i < n_weeks; i++){
      cum_sum_cases[i] = cum_sum_cases[i-1] + new_cases[i];
  }

  std::vector<double> times(n_weeks,0.); //convert time from weeks to days
  for (unsigned int i = 0; i < times.size(); i++){
    times[i] = weeks[i] * 7;
  }

  //Set initial values
  //S_h, E_h, I_h, R_h, S_v E_v, I_v, C
  double nh = 206 * pow(10,6);
  double nv = 1;
  double ci = rep_factor * 8201.0;
  double ehi = ci;
  double ihi = ci;
  double rhi = 29639.0;
  double shi = nh - ehi - ihi - rhi;
  double ivi = 0.00022;
  double evi = ivi;
  double svi = nv - evi - ivi;
  //set initial values
  initialValues[0] = shi;
  initialValues[1] = ehi;
  initialValues[2] = ihi;
  initialValues[3] = rhi;
  initialValues[4] = svi;
  initialValues[5] = evi;
  initialValues[6] = ivi;
  initialValues[7] = ci;

  fclose(dataFile);
  std::cout << "The number of data points is " << n_weeks << "\n\n";

  //create dummy vector to be filled with params inside likelihood
  std::vector<double> queso_params(n_params, 0.0); 

  //------------------------------------------------------
  // SIP Step 1 of 6: Instantiate the parameter space
  //------------------------------------------------------
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> paramSpace(env, "param_", n_params, NULL);

  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinValues(paramSpace.zeroVector());
  QUESO::GslVector paramMaxValues(paramSpace.zeroVector());
  //set upper and lower limits for parameters
  for (unsigned int i=0; i<n_params; ++i){
    paramMinValues[i] = -.3;//-INFINITY;
    paramMaxValues[i] = +.15;//INFINITY;
}
  //TODO: would need something similar if hyperparameters...
  /* //variance of xi */
  /* for (unsigned int i=n_xi; i<2*n_xi; ++i){ */
  /*   paramMinValues[i] = -INFINITY; */
  /*   paramMaxValues[i] = INFINITY;} */
  /* //xi */
  /* for (unsigned int i=2*n_xi; i<3*n_xi; ++i){ */
  /*   paramMinValues[i] = -INFINITY; */
  /*   paramMaxValues[i] = INFINITY;} */
  //mean of a_0 and b_0

  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);

  // collect information about dynamical system
  dynamics_info dynMain(n_s, n_weeks, inad_type, params_factor, queso_params);

  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function 
  // object to be used by QUESO.
  //------------------------------------------------------
  likelihoodRoutine_Data likelihoodRoutine_Data1(env, times, initialValues, cum_sum_cases, var, &dynMain);

  QUESO::GenericScalarFunction<>
    likelihoodFunctionObj(
        "like_",
			  paramDomain,
			  likelihoodRoutine,
        static_cast<void *> (&likelihoodRoutine_Data1),
			  true); // the routine computes [ln(function)]
    
  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  // UNIFORM PRIOR
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_",paramDomain);

  QUESO::GslVector diagVec(paramSpace.zeroVector());
  // NORMAL PRIOR
  /* QUESO::GslVector meanVec(paramSpace.zeroVector()); */
  /* /1* for (unsigned int i = 0; i < n_params; i++) meanVec[i] = 1.; *1/ */

  /* for (unsigned int i = 0; i < n_params; i++) diagVec[i] = 1.; */
  /* QUESO::GslMatrix covMatrix(diagVec); */
  /* QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix> */
  /*   priorRv("prior_",paramDomain,meanVec,covMatrix); */

  //------------------------------------------------------
  // SIP Step 5 of 6: Instantiate the inverse problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<>
    postTotal("post_",  // Extra prefix before the default "rv_" prefix
           paramSpace);
        
  QUESO::StatisticalInverseProblem<>
    ip("",          // No extra prefix before the default "ip_" prefix
       NULL,
       priorRv,
       likelihoodFunctionObj,
       postTotal);

  //------------------------------------------------------
  // SIP Step 6 of 6: Solve the inverse problem, that is,
  // set the 'pdf' and the 'realizer' of the posterior RV //------------------------------------------------------
  std::cout << "Solving the SIP with Multi-Level Metropolis Hastings" 
	    << std::endl << std::endl;  

  //The following is set if use ip.solveWithBayesMetropolisHastings
  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  for (unsigned int i = 0; i < n_params; i++) {paramInitials[i] = 0;}
   //
//priorRv.realizer().realization(paramInitials);

  /* QUESO::GslVector diagVec(paramSpace.zeroVector()); */
  QUESO::GslMatrix proposalCovMatrix(diagVec);
  for (unsigned int i = 0; i < n_params; i++) proposalCovMatrix(i,i) = 1.e-4;
  //proposalCovMatrix(0,0) = 1e-6;
  //proposalCovMatrix(1,1) = 1e-6;
  //proposalCovMatrix(2,2) = 1e-6;
  //proposalCovMatrix(3,3) = 1e-6;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  /* ip.seedWithMAPEstimator(); */
  /* ip.solveWithBayesMetropolisHastings(); */

  /* ip.solveWithBayesMLSampling(); */

  //================================================================
  // Statistical forward problem (SFP)
  //================================================================
  gettimeofday(&timevalNow, NULL);
  std::cout << "Beginning 'SFP -> Undecided QoI' at " 
            << ctime(&timevalNow.tv_sec)
            << std::endl;

  //------------------------------------------------------
  // SFP Step 1 of 6: Instantiate the parameter *and* qoi spaces. 
  // SFP input RV = FIP posterior RV, so SFP parameter space
  // has been already defined.
  //------------------------------------------------------
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> qoiSpace(env, "qoi_", n_weeks * dim, NULL);

  //------------------------------------------------------
  // SFP Step 2 of 6: Instantiate the parameter domain 
  //------------------------------------------------------
  
  // Not necessary because input RV of the SFP = output RV of SIP. 
  // Thus, the parameter domain has been already defined.
  
  //------------------------------------------------------ 
  // SFP Step 3 of 6: Instantiate the qoi function object 
  // to be used by QUESO.
  //------------------------------------------------------
  qoiRoutine_Data qoiRoutine_Data(env, times, initialValues, &dynMain);
  
  QUESO::GenericVectorFunction<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    qoiFunctionObj("qoi_",
                   paramDomain,
                   qoiSpace,
                   qoiRoutine,
                   (void *) &qoiRoutine_Data);
      
  //------------------------------------------------------
  // SFP Step 4 of 6: Define the input RV
  //------------------------------------------------------
  
  // Not necessary because input RV of SFP = output RV of SIP 
  // (postRv).
      
  //------------------------------------------------------
  // SFP Step 5 of 6: Instantiate the forward problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix> qoiRv("qoi_", qoiSpace);
  
  QUESO::StatisticalForwardProblem<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    fp("",
       NULL,
       postTotal,
       qoiFunctionObj,
       qoiRv);

  //------------------------------------------------------
  // SFP Step 6 of 6: Solve the forward problem
  //------------------------------------------------------
  std::cout << "Solving the SFP with Monte Carlo" 
            << std::endl << std::endl;  
  fp.solveWithMonteCarlo(NULL);

  //------------------------------------------------------
  gettimeofday(&timevalNow, NULL);
  if ((env.subDisplayFile()       ) && 
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Ending run of 'Zika plus discrepancy' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  if (env.fullRank() == 0) {
    std::cout << "Ending run of 'Zika plus discrepancy' example at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  return;
}
