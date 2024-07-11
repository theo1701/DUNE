///////////////////////////////////////////////////////////////////////////
//	This code was prepared by Pedro Pasquini in september 2018      ///
//									///
//		please cite: 1605.01670 and 1802.02133			///
//		if you use this for scientific porpouse			///
//									///
//									///
///////////////////////////////////////////////////////////////////////////
//	This is an exemplo on how to use NON_UNITARY_prob.h library	///
//									///
//		NOTICE: This follows the notation of 1605.01670 	///
//	and does not include a normalization due to near detector 	///
// 	  for a discussion to see if you need or not to worry about	///
//	        this normalization see this: 1609.08637			///
///////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   // GLoBES library 

// #include "NON_UNITARY_prob.h" // Implementation of non-unitary


/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/
// char MYFILE[]="prob_NH_dcp-90_uni.dat"; // output
char MYFILE[] = "prob_IH_dcp-90_uni_mubarmubar.dat";
int main(int argc, char *argv[])
{ 

// Initialise and define filename of output chains
  FILE *outfile = NULL; 

  outfile = fopen(MYFILE, "w");
  if (outfile == NULL)
  {
    printf("Error opening output file.\n");
    return -1;
  }
/********************************
     INITIALIZING GLOBES
******************************/

  // Initialize libglobes  
  glbInit(argv[0]);

  // Register non-standard probability engine. This has to be done
   // before any calls to glbAllocParams() or glbAllocProjections() 
  // glbRegisterProbabilityEngine(15,      // Number of parameters 
  //                              &my_probability_matrix,
  //                              &my_set_oscillation_parameters,
  //                              &my_get_oscillation_parameters,
  //                              NULL);


  char AEDLFILE1[]="2020_nova_app2.248+2.066.glb";
  glbInitExperiment(AEDLFILE1,&glb_experiment_list[0],&glb_num_of_exps) ; 
 

/********************************
     DECLARING VARIABLES
******************************/



double M_PI_DEF = 3.1415926535;

//-------------------------
// IH
double true_theta12 = asin(sqrt(0.307));
double true_theta13 = asin(sqrt(0.02222));
double true_theta23 = asin(sqrt(0.568));
double true_sdm = 7.41e-5;
double true_ldm = -2.487e-3;

// Define standard oscillation parameters (cf. https://arxiv.org/pdf/1405.7540v3.pdf) 
// double true_theta12 = asin(sqrt(0.310));
// double true_theta13 = asin(sqrt(0.02237));
// double true_theta23 = asin(sqrt(0.563));
double true_deltacp = -90*M_PI_DEF/180;

//   // Define Non-Unitary PARAMETERS

// // MODULUES
// // Notice, for standard 3-nu oscillation alphaii=1 other zero.
//   double true_alpha00 =1.0, true_alpha10 =0.0, true_alpha20 =0.0, true_alpha11 =1.0,true_alpha21 =0.0, true_alpha22 =1.0;
// // PHASES
//   double true_phi10 = 0.0*M_PI_DEF, true_phi20 =0.0, true_phi21 = 0.0;

  /* Initialize parameter and projection vector(s) */
  glb_params true_values = glbAllocParams();
  //glb_params test_values = glbAllocParams();
  //glb_params output = glbAllocParams();
  //glb_projection NSI_chi = glbAllocProjection(); 

/********************************
     DECLARING TRUE
******************************/ 

// Standard parameter 
  glbDefineParams(true_values,true_theta12,true_theta13,true_theta23,true_deltacp,true_sdm,true_ldm);
  
/*glbSetOscillationParameters(true_values);
  glbSetRates();*/

//Non-unitary parameters
  // glbSetOscParams(true_values,true_alpha00, GLB_ALPHA_00); 
  // glbSetOscParams(true_values,true_alpha10, GLB_ALPHA_10);   
  // glbSetOscParams(true_values,true_phi10, GLB_PHI_10);   
  // glbSetOscParams(true_values,true_alpha20, GLB_ALPHA_20);
  // glbSetOscParams(true_values,true_phi20, GLB_PHI_20);
  // glbSetOscParams(true_values,true_alpha11, GLB_ALPHA_11);
  // glbSetOscParams(true_values,true_alpha21, GLB_ALPHA_21);
  // glbSetOscParams(true_values,true_phi21, GLB_PHI_21);
  // glbSetOscParams(true_values,true_alpha22, GLB_ALPHA_22);
glbSetDensityParams(true_values,1.0,GLB_ALL);
//fprintf(stdout, "%g %g \n", true_alpha00,true_alpha11);
/* Calculation of probability */
glbSetOscillationParameters(true_values);
  glbSetRates();
double e;
double p;

// nu_1 is nu_e, nu_2 is nu_mu, and nu_3 is nu_tau
for(e=0;e<=5;e +=.001)
{
p=glbProfileProbability(0,2, 2, -1,e);
// q=glbProfileProbability(0,2, 1, -1,e);
// r=glbProfileProbability(0,2, 2, +1,e); 
// s=glbProfileProbability(0,2, 2, -1,e);
// fprintf(outfile, "%g %g %g %g %g \n", e,p, q, r, s);
fprintf(outfile, "%g, %g \n", e, p);
fprintf(stdout, "%g, %g \n", e, p);
}
  

/********************************
     DECLARING TEST
******************************/ 
	//glbCopyParams(true_values,test_values);


/********************************
     DECLARING ERROR
******************************/ 
  /*glb_params input_errors = glbAllocParams();
  glbDefineParams(input_errors,0,0.0026,0,0, 0,0);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
// Non-unitary parameters
  glbSetOscParams(input_errors,0, GLB_ALPHA_00); 
  glbSetOscParams(input_errors,0, GLB_ALPHA_10);   
  glbSetOscParams(input_errors,0, GLB_PHI_10);   
  glbSetOscParams(input_errors,0, GLB_ALPHA_20);
  glbSetOscParams(input_errors,0, GLB_PHI_20);
  glbSetOscParams(input_errors,0, GLB_ALPHA_11);
  glbSetOscParams(input_errors,0, GLB_ALPHA_21);
  glbSetOscParams(input_errors,0, GLB_PHI_21);
  glbSetOscParams(input_errors,0, GLB_ALPHA_22);*/



  // Destroy parameter and projection vector(s) 
  glbFreeParams(true_values);
  
  exit(0);
}


 
