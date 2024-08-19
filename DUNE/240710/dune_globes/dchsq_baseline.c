/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * Example: Find the sgn(dm31^2)-degeneracy and compute shape with glbChiSys
 * Compile with ``make example3''
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h> /* GLoBES library */
#include "myio.h"          /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE[] = "dchsq_baseline.dat";

int main(int argc, char *argv[])
{
  /* Initialize libglobes */
  glbInit(argv[0]);

  /* Initialize experiment NFstandard.glb */
  glbInitExperiment("DUNE_GLoBES.glb", &glb_experiment_list[0], &glb_num_of_exps);
  /*glbInitExperiment("2020_T2K_5nu_app.glb",&glb_experiment_list[0],&glb_num_of_exps);
 glbInitExperiment("T2K_2020_5disapp.glb",&glb_experiment_list[0],&glb_num_of_exps);*/

  /* Intitialize output */
  //  InitOutput(MYFILE,"Format: Log(10,s22th13)   deltacp   chi^2 \n");
  FILE *outfile = NULL;

  outfile = fopen(MYFILE, "w");
  if (outfile == NULL)
  {
    printf("Error opening output file.\n");
    return -1;
  }

  double my_M_PI=3.1415926535;

  /* Define standard oscillation parameters */
  // double theta12 = asin(sqrt(0.31));
  // double theta13 = asin(sqrt(0.084)) / 2;
  // double theta23 = asin(sqrt(0.5));
  // double deltacp = -170 * M_PI / 180;
  // double sdm = 7.41e-5; // used
  // double atmdm = 2.44e-3;
  // double ldm = atmdm + pow(cos(theta12), 2) * sdm - cos(deltacp) * sin(theta13) * sin(2 * theta12) * tan(theta23) * sdm;
  // double m;

  double true_theta12 = 33.5 * my_M_PI / 180;
  double true_theta13 = 8.5 * my_M_PI / 180;
  double true_theta23 = 45 * my_M_PI / 180;
  // double true_deltacp = -90 * my_M_PI / 180;
  double true_sdm = 7.5e-5;
  double true_ldm = 2.45e-3;
  double ih_ldm = -2.46e-3;

  double min_chsq, res;

  /*  */
  for (double true_deltacp = -my_M_PI; true_deltacp <= my_M_PI; true_deltacp += 0.05)
  {
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();

    glbDefineParams(true_values, true_theta12, true_theta13, true_theta23, true_deltacp, true_sdm, true_ldm);
    glbSetDensityParams(true_values, 1.0, GLB_ALL);

    glbDefineParams(test_values, true_theta12, true_theta13, true_theta23, true_deltacp, true_sdm, ih_ldm);
    glbSetDensityParams(test_values, 1.0, GLB_ALL);

    glbSetOscillationParameters(true_values);
    glbSetRates();

    glbFreeParams(true_values);
    glbFreeParams(test_values);
  }

  double min_chsq, res;

  for (double true_deltacp = -my_M_PI; true_deltacp <= my_M_PI; true_deltacp += 0.05)
  {
    min_chsq = 1e10;
    glbSetOscParams(true_values, true_deltacp, GLB_DELTA_CP);
    glbSetDensityParams(true_values, 1.0, GLB_ALL);
    glbSetOscillationParameters(true_values);
    glbSetRates();

    for (double dtest = -my_M_PI; dtest <= my_M_PI; dtest += 0.05)
    {
      glbSetOscParams(test_values, dtest, GLB_DELTA_CP);
      glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);
      res = glbChiSys(test_values, GLB_ALL, GLB_ALL);

      if (res < min_chsq)
      {
        min_chsq = res;
      }
    }
    fprintf(stdout, "%lf %lf \n", true_deltacp, min_chsq);
    fprintf(outfile, "%lf %lf \n", true_deltacp, min_chsq);
    
  }

  exit(0);
}
