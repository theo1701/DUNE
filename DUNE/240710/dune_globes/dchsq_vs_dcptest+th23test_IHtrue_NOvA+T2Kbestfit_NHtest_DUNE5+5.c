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
char MYFILE[] = "dchsq_vs_dcptest+th23test_IHtrue_NOvA+T2Kbestfit_dcp+90_NHtest_DUNE5+5.dat";

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

  // double M_PI=3.1415926535;

  /* Define standard oscillation parameters */
  // double theta12 = asin(sqrt(0.31));
  // double theta13 = asin(sqrt(0.084)) / 2;
  // double theta23 = asin(sqrt(0.5));
  // double deltacp = -170 * M_PI / 180;
  // double sdm = 7.41e-5; // used
  // double atmdm = 2.44e-3;
  // double ldm = atmdm + pow(cos(theta12), 2) * sdm - cos(deltacp) * sin(theta13) * sin(2 * theta12) * tan(theta23) * sdm;
  // double m;

  double theta12 = asin(sqrt(0.307));
  double theta13 = asin(sqrt(0.02222));
  double theta23 = asin(sqrt(0.5));
  double deltacp = +90 * M_PI / 180;
  double sdm = 7.41e-5; // used
  double atmdm = -2.487e-3;
  double ldm = atmdm + sdm;
  double m;

  /* Initialize parameter vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params central_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params deg_pos = glbAllocParams();

  /* The simulated data are computed */
  glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm);
  glbSetDensityParams(true_values, 1.0, GLB_ALL);

  glbDefineParams(deg_pos, theta12, theta13, theta23, deltacp, sdm, ldm);
  glbSetDensityParams(deg_pos, 1.0, GLB_ALL);

  glbSetOscillationParameters(true_values);
  glbSetRates();

  //   glbGetOscParams(deg_pos,GLB_THETA_13);
  // glbGetOscParams(deg_pos,GLB_DELTA_CP);

  /* If degeneracy at low enough confidence level: compute section */

  double thetheta13, x, y, z, k, l, res, theldm, thetheta23, thedeltacp, theatmdm;

  /*glbInitExperiment("2019_nova_app5+5.glb",&glb_experiment_list[0],&glb_num_of_exps);
   glbInitExperiment("2019_nova_disapp5+5.glb",&glb_experiment_list[0],&glb_num_of_exps);
   glbInitExperiment("2019_T2K_5nu_app.glb",&glb_experiment_list[0],&glb_num_of_exps);
   glbInitExperiment("2019_T2K_5anu_app.glb",&glb_experiment_list[0],&glb_num_of_exps);
   glbInitExperiment("T2K_2018_5+5disapp.glb",&glb_experiment_list[0],&glb_num_of_exps);*/

  double s = 0.00058; // (0.02397 - 0.02047) / 6
  double a = 0.01e-3;
  //
  for (y = -180; y <= 180; y = y + 10)
  {
    for (l = 0.41; l <= 0.61; l = l + 0.01)
    {
      m = 10000;
      for (x = 0.02047; x <= 0.02397; x = x + s)
      {
        // for (k = 2.32e-3 - 3 * a; k <= 2.32e-3 + 3.01 * a; k = k + a)
        for (k = 2.426e-3; k <= 2.586e-3; k = k + a)
        { /* Set vector of test values */
          thetheta23 = asin(sqrt(l));
          thedeltacp = y * M_PI / 180.0;
          thetheta13 = asin(sqrt(x));
          // theatmdm = -k;
          glbSetOscParams(deg_pos, thetheta23, GLB_THETA_23);
          glbSetOscParams(deg_pos, thedeltacp, GLB_DELTA_CP);
          glbSetOscParams(deg_pos, thetheta13, GLB_THETA_13);
          // theldm = theatmdm + pow(cos(theta12), 2) * sdm - cos(thedeltacp) * sin(thetheta13) * sin(2 * theta12) * tan(thetheta23) * sdm;
          theldm = k;
          glbSetOscParams(deg_pos, theldm, GLB_DM_31);

          /* Compute Chi^2 for all loaded experiments and all rules */
          glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);
          res = glbChiSys(deg_pos, GLB_ALL, GLB_ALL);
          double prior13 = pow(((x - 0.02224) / s), 2);

          res = res + prior13;

          if (res < m)
          {

            // marginalization;
            m = res;
          }
        }
      }

      fprintf(stdout, "%lf %lf %lf \n", y, l, m);
      fprintf(outfile, "%lf %lf %lf \n", y, l, m);
    }
    fprintf(outfile, "\n");
  }

  /* Destroy parameter vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(central_values);
  glbFreeParams(input_errors);
  glbFreeParams(deg_pos);

  exit(0);
}
