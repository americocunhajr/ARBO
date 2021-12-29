/*-------------------------------------------------------------------
 * ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
 *-----------------------------------------------------------------*/

/*-------------------------------------------------------------------
 * Brief description of this file: 
 * testing the model output of the reduced model
 *-----------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <random>
#include "model.h"
#include "dynamics_info.h"

int main()
{

  //read in info file
  //FILE *infoFile;
  //infoFile = fopen("./inputs/info.txt","r");
  //if(fscanf(infoFile,"%u %u %u %u %u %lf %u", &n_S, &n_s, &n_phis_cal, &n_phis_val, &n_times, &var, &inad_type)){};
  //fclose(infoFile);

  //std::cout << "var = " << var << "\n";

  unsigned int n_s = 7;  //number of species in reduced model
  unsigned int inad_type = 0;
  unsigned int params_factor;
    if( inad_type == 0 ) { params_factor = 1;}
    if( inad_type == 1 ) { params_factor = 2;}
  unsigned int n_delta = params_factor*n_s;         //the model inadequacy terms, should be set to zero
  unsigned int n_weeks = 52;      //no hyperparameters, for now
  double finalTime = 364.0;

  //dummy vector for now
  std::vector<double> delta(n_delta,0.);
  /* dynamics_info dynMain(n_s, n_times, delta); */
  dynamics_info dynMain(n_s, n_weeks, inad_type, params_factor, delta);
  
  std::vector<double> phiPoints(1,0.);
  for (unsigned int i = 0; i < 1; i++){
  phiPoints[i] = i; //phiPoints[1] = 1.0; phiPoints[2] = 1.1;
  }

  //double timePoint = 2.0e-5;

  std::vector<double> initialValues(n_s + 1,0.); 
  //S_h, E_h, I_h, R_h, S_v E_v, I_v, C
  double nh = 206 * pow(10,6);
  double nv = 1;
  double ci = 8201.0;
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

  std::vector<double> timePoints(n_weeks,0.);
  //for (unsigned int i = 0; i < timePoints.size(); i++){
    //timePoints[i] = (i) * finalTime/n_times;
    //std::cout << "timePoints["<<i<<"] = "<<timePoints[i]<<"\n";
    /* timePoints[i] = (i) * 2.5e-5 + 3.0e-5; */
    /* timePoints[i] = (i+1) * 5.0e-6; */
  //}
  //std::cout<< "timePoints end = " << timePoints[51] << "\n\n";
  //std::cout<< "timePoints size = " << timePoints.size() << "\n\n";
  for (unsigned int i = 0; i < timePoints.size(); i++){
    timePoints[i] = (i+1) * 7;
  }

  //open file now to write data
  std::ofstream datafile;
  datafile.open ("./inputs/data-red.txt");

    //return 100 time points of all species and temp
    std::vector<double> returnValues((n_s+1)*timePoints.size(),0.);

    zikaComputeModel(initialValues,timePoints,&dynMain,returnValues);
   
    //create measurement error
    //std::normal_distribution<double> distribution(0.0,std::sqrt(var));
    //TODO increase this appropriately for temperature measurements

    for (unsigned int i = 0; i < timePoints.size(); i++){
        /* std::cout << "i = " << i << "\n"; */
      
        for (unsigned int j = 0; j < n_s+1; j++){
          double error = 0.; //distribution(generator); //don't need to add error here
          double measurement = returnValues[(n_s+1)*i + j] + error;
          if (measurement>0){datafile << phiPoints[0] << " " << timePoints[i] << " " << measurement << "\n";}
          else {datafile << phiPoints[0] << " " << timePoints[i] << " " <<  0.0 << "\n";}
      } //datafile << "\n";
    }
  return 0;
}
