
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "3x3-C/zhetrd3.h"
#include "3x3-C/zheevq3.h"

#include <globes/globes.h> /* GLoBES library */

// Macros
#define SQR(x) ((x) * (x))                         // x^2
#define SQR_ABS(x) (SQR(creal(x)) + SQR(cimag(x))) // |x|^2

/* Output file */

#define GLB_EPS_E_E 6

// #define GLB_EPS_E_E 6
// #define GLB_EPS_E_M 7           /* Index of non-standard parameter */
// #define GLB_EPS_E_T 8            /* Index of non-standard parameter */
// #define GLB_PHI_E_M 9            /* Index of non-standard parameter */
// #define GLB_PHI_E_T 10

/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/

double th12;
double th13;
double th23;
// double eps_e_e;
// double eps_e_m;
double eps_e_e;

double dcp;

// double phi_e_m;

double d21;
double d31;

/***************************************************************************
 * Store oscillation parameters in internal data structures.               *
 * For more sophisticated probability engines, this would be the right     *
 * place to pre-compute the mixing matrix and parts of the Hamiltonian in  *
 * order to speed up the calls to the actual probability matrix function.  *
 ***************************************************************************/
int my_set_oscillation_parameters(glb_params p, void *user_data)
{
  th12 = glbGetOscParams(p, GLB_THETA_12);
  th13 = glbGetOscParams(p, GLB_THETA_13);
  th23 = glbGetOscParams(p, GLB_THETA_23);
  // eps_e_e = glbGetOscParams(p, GLB_EPS_E_E);
  // eps_e_m   = glbGetOscParams(p, GLB_EPS_E_M);
  dcp = glbGetOscParams(p, GLB_DELTA_CP);
  // phi_e_m    = glbGetOscParams(p, GLB_PHI_E_M);
  d21 = glbGetOscParams(p, GLB_DM_21); /* Convert to GeV^2 */
  d31 = glbGetOscParams(p, GLB_DM_31); /* Convert to GeV^2 */
  eps_e_e = glbGetOscParams(p, GLB_EPS_E_E);

  return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/
int my_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p, th12, GLB_THETA_12);
  glbSetOscParams(p, th13, GLB_THETA_13);
  glbSetOscParams(p, th23, GLB_THETA_23);
  // glbSetOscParams(p, eps_e_e, GLB_EPS_E_E);
  // glbSetOscParams(p, eps_e_m, GLB_EPS_E_M);
  glbSetOscParams(p, dcp, GLB_DELTA_CP);

  // glbSetOscParams(p, phi_e_m, GLB_PHI_E_M);
  glbSetOscParams(p, d21, GLB_DM_21); /* Convert to eV^2 */
  glbSetOscParams(p, d31, GLB_DM_31); /* Convert to eV^2 */
  glbSetOscParams(p, eps_e_e, GLB_EPS_E_E);

  return 0;
}

/***************************************************************************
 * Calculate oscillation probabilities.                                    *
 * Since for our setup, only P_ee is required, all other entries of P are  *
 * set to zero for simplicity. Furthermore, we neglect matter effects and  *
 * the filter feature (parameter filter_sigma).                            *
 * The formula for P_ee is Eq. (36) from hep-ph/0502147.                   *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:            The buffer where the probabilities are to be stored     *
 *   cp_sign:      +1 if probalities for neutrinos are requested, -1 for   *
 *                 anti-neutrinos.                                         *
 *   E:            The neutrino energy in GeV                              *
 *   psteps:       Number of constant density layers in the matter profile *
 *   length:       The lengths of these layers in km                       *
 *   density:      The individual densities of these layers in g/cm^3      *
 *   filter_sigma: Width of low-pass filter as given in the AEDL file      *
 ***************************************************************************/

int my_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

  double complex A[3][3];
  double complex Q[3][3];
  double w[3];
  double rho, Ve;
  double sin12, sin13, sin23, cos12, cos13, cos23;
  double sinb12, sinb13, sinb23, cosb12, cosb13, cosb23;

  int z, L, ii, jj;
  /* Set all probabilities to zero initially */
  for (ii = 0; ii < 3; ii++)
  {
    for (jj = 0; jj < 3; jj++)
    {
      P[ii][jj] = 0.0;
    }
  }

  rho = 2.8;

  /* Calculate total baseline */
  L = 0.0;

  double complex Pee[psteps];
  double complex Pem[psteps];
  double complex Pet[psteps];

  double complex Pme[psteps];
  double complex Pmm[psteps];
  double complex Pmt[psteps];

  double complex Pte[psteps];
  double complex Ptm[psteps];
  double complex Ptt[psteps];

  for (z = 0; z < psteps; z++)
  {
    L += length[z];
    // L = GLB_KM_TO_EV(L) * 1.0e9;      /* Convert to GeV^{-1} */

    // double L = 810 * 1.0e9;

    sin12 = sin(th12);
    sin13 = sin(th13);
    sin23 = sin(th23);
    cos12 = cos(th12);
    cos13 = cos(th13);
    cos23 = cos(th23);

    double complex Ue1 = cos12 * cos13;
    double complex Ue2 = sin12 * cos13;
    double complex Ue3 = sin13 * cexp(-I * cp_sign * dcp);

    double complex Um1 = -sin12 * cos23 - cos12 * sin23 * sin13 * cexp(I * cp_sign * dcp);
    double complex Um2 = cos12 * cos23 - sin12 * sin23 * sin13 * cexp(I * cp_sign * dcp);
    double complex Um3 = sin23 * cos13;

    double complex Ut1 = sin12 * sin23 - cos12 * cos23 * sin13 * cexp(I * cp_sign * dcp);
    double complex Ut2 = -cos12 * sin23 - sin12 * cos23 * sin13 * cexp(I * cp_sign * dcp);
    double complex Ut3 = cos23 * cos13;

    Ve = cp_sign * 0.5 * 0.000076 * density[z]; // sqrt(2) * G_F * N_e

    double complex eps00 = Ve * eps_e_e;
    double complex eps01 = 0;
    double complex eps02 = 0;

    double complex eps10 = 0;
    double complex eps11 = 0;
    double complex eps12 = 0;

    double complex eps20 = 0;
    double complex eps21 = 0;
    double complex eps22 = 0;

    // double complex eps00 = cp_sign * eps_e_e * 1e18;
    // double complex eps01 = cp_sign * eps_e_m * cexp(I * cp_sign * phi_e_m) * 1e18;
    // double complex eps02 = cp_sign * eps_e_t * cexp(I * cp_sign * phi_e_t) * 1e18;

    // double complex eps10 = conj(eps01);
    // b double complex eps11 = 0;
    // double complex eps12 = 0;

    // double complex eps20 = conj(eps02);
    // double complex eps21 = 0;
    // double complex eps22 = 0;

    // fprintf(stdout, "%g %g %g\n", cp_sign*0.5*0.000076*rho, cp_sign*0.5*0.000076*density[z], length[z] );

    A[0][0] = (0.5 / E) * (Ue2 * conj(Ue2) * d21 + Ue3 * conj(Ue3) * d31) + Ve + eps00;
    A[0][1] = (0.5 / E) * (Ue2 * conj(Um2) * d21 + Ue3 * conj(Um3) * d31) + eps01;
    A[0][2] = (0.5 / E) * (Ue2 * conj(Ut2) * d21 + Ue3 * conj(Ut3) * d31) + eps02;

    A[1][0] = (0.5 / E) * (Um2 * conj(Ue2) * d21 + Um3 * conj(Ue3) * d31) + eps10;
    A[1][1] = (0.5 / E) * (Um2 * conj(Um2) * d21 + Um3 * conj(Um3) * d31) + eps11;
    A[1][2] = (0.5 / E) * (Um2 * conj(Ut2) * d21 + Um3 * conj(Ut3) * d31) + eps12;

    A[2][0] = (0.5 / E) * (Ut2 * conj(Ue2) * d21 + Ut3 * conj(Ue3) * d31) + eps20;
    A[2][1] = (0.5 / E) * (Ut2 * conj(Um2) * d21 + Ut3 * conj(Um3) * d31) + eps21;
    A[2][2] = (0.5 / E) * (Ut2 * conj(Ut2) * d21 + Ut3 * conj(Ut3) * d31) + eps22;

    zheevq3(A, Q, w);

    double L1 = w[0];
    double L2 = w[1];
    double L3 = w[2];

    double complex Ue1f = Q[0][0];
    double complex Ue2f = Q[0][1];
    double complex Ue3f = Q[0][2];

    double complex Um1f = Q[1][0];
    double complex Um2f = Q[1][1];
    double complex Um3f = Q[1][2];

    double complex Ut1f = Q[2][0];
    double complex Ut2f = Q[2][1];
    double complex Ut3f = Q[2][2];

    Pee[z] = conj(Ue1f) * Ue1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Ue2f) * Ue2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Ue3f) * Ue3f * cexp(-I * 4 * 1.27 * L3 * length[z]);
    Pem[z] = conj(Ue1f) * Um1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Ue2f) * Um2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Ue3f) * Um3f * cexp(-I * 4 * 1.27 * L3 * length[z]);
    Pet[z] = conj(Ue1f) * Ut1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Ue2f) * Ut2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Ue3f) * Ut3f * cexp(-I * 4 * 1.27 * L3 * length[z]);

    Pme[z] = conj(Um1f) * Ue1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Um2f) * Ue2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Um3f) * Ue3f * cexp(-I * 4 * 1.27 * L3 * length[z]);
    Pmm[z] = conj(Um1f) * Um1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Um2f) * Um2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Um3f) * Um3f * cexp(-I * 4 * 1.27 * L3 * length[z]);
    Pmt[z] = conj(Um1f) * Ut1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Um2f) * Ut2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Um3f) * Ut3f * cexp(-I * 4 * 1.27 * L3 * length[z]);

    Pte[z] = conj(Ut1f) * Ue1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Ut2f) * Ue2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Ut3f) * Ue3f * cexp(-I * 4 * 1.27 * L3 * length[z]);
    Ptm[z] = conj(Ut1f) * Um1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Ut2f) * Um2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Ut3f) * Um3f * cexp(-I * 4 * 1.27 * L3 * length[z]);
    Ptt[z] = conj(Ut1f) * Ut1f * cexp(-I * 4 * 1.27 * L1 * length[z]) + conj(Ut2f) * Ut2f * cexp(-I * 4 * 1.27 * L2 * length[z]) + conj(Ut3f) * Ut3f * cexp(-I * 4 * 1.27 * L3 * length[z]);
  }

  P[0][0] = SQR_ABS(Pee[0]);
  P[0][1] = SQR_ABS(Pem[0]);
  P[0][2] = SQR_ABS(Pet[0]);

  P[1][0] = SQR_ABS(Pme[0]);
  P[1][1] = SQR_ABS(Pmm[0]);
  P[1][2] = SQR_ABS(Pmt[0]);

  P[2][0] = SQR_ABS(Pte[0]);
  P[2][1] = SQR_ABS(Ptm[0]);
  P[2][2] = SQR_ABS(Ptt[0]);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
inline void zhetrd3(double complex A[3][3], double complex Q[3][3],
                    double d[3], double e[2])
// ----------------------------------------------------------------------------
// Reduces a hermitian 3x3 matrix to real tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
  const int n = 3;
  double complex u[n], q[n];
  double complex omega, f;
  double K, h, g;
  int i, j;

  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (i = 0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (j = 0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  // Bring first row and column to the desired form
  h = SQR_ABS(A[0][1]) + SQR_ABS(A[0][2]);
  if (creal(A[0][1]) > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f = g * A[0][1];
  u[1] = conj(A[0][1]) - g;
  u[2] = conj(A[0][2]);

  omega = h - f;
  if (creal(omega) > 0.0)
  {
    omega = 0.5 * (1.0 + conj(omega) / omega) / creal(omega);
    K = 0.0;
    for (i = 1; i < n; i++)
    {
      f = conj(A[1][i]) * u[1] + A[i][2] * u[2];
      q[i] = omega * f;           // p
      K += creal(conj(u[i]) * f); // u* A u
    }
    K *= 0.5 * SQR_ABS(omega);

    for (i = 1; i < n; i++)
      q[i] = q[i] - K * u[i];

    d[0] = creal(A[0][0]);
    d[1] = creal(A[1][1]) - 2.0 * creal(q[1] * conj(u[1]));
    d[2] = creal(A[2][2]) - 2.0 * creal(q[2] * conj(u[2]));

    // Store inverse Householder transformation in Q
#ifndef EVALS_ONLY
    for (j = 1; j < n; j++)
    {
      f = omega * conj(u[j]);
      for (i = 1; i < n; i++)
        Q[i][j] = Q[i][j] - f * u[i];
    }
#endif

    // Calculate updated A[1][2] and store it in f
    f = A[1][2] - q[1] * conj(u[2]) - u[1] * conj(q[2]);
  }
  else
  {
    for (i = 0; i < n; i++)
      d[i] = creal(A[i][i]);
    f = A[1][2];
  }

  // Make (23) element real
  e[1] = cabs(f);
#ifndef EVALS_ONLY
  if (e[1] != 0.0)
  {
    f = conj(f) / e[1];
    for (i = 1; i < n; i++)
      Q[i][n - 1] = Q[i][n - 1] * f;
  }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
int zheevq3(double complex A[3][3], double complex Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to real tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   zhetrd3()
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double e[3];                // The third element is used only as temporary workspace
  double g, r, p, f, b, s, c; // Intermediate storage
  double complex t;
  int nIter;
  int m;
  int l, i, j, k;

  // Transform A to real tridiagonal form by the Householder method
  zhetrd3(A, Q, w, e);

  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
  for (l = 0; l < n - 1; l++)
  {
    nIter = 0;
    while (1)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m = l; m <= n - 2; m++)
      {
        g = fabs(w[m]) + fabs(w[m + 1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;

      if (nIter++ >= 30)
        return -1;

      // Calculate g = d_m - k
      g = (w[l + 1] - w[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = w[m] - w[l] + e[l] / (g + r);
      else
        g = w[m] - w[l] + e[l] / (g - r);

      s = c = 1.0;
      p = 0.0;
      for (i = m - 1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c = g / f;
          r = sqrt(SQR(c) + 1.0);
          e[i + 1] = f * r;
          c *= (s = 1.0 / r);
        }
        else
        {
          s = f / g;
          r = sqrt(SQR(s) + 1.0);
          e[i + 1] = g * r;
          s *= (c = 1.0 / r);
        }

        g = w[i + 1] - p;
        r = (w[i] - g) * s + 2.0 * c * b;
        p = s * r;
        w[i + 1] = g + p;
        g = c * r - b;

        // Form eigenvectors
#ifndef EVALS_ONLY
        for (k = 0; k < n; k++)
        {
          t = Q[k][i + 1];
          Q[k][i + 1] = s * Q[k][i] + c * t;
          Q[k][i] = c * Q[k][i] - s * t;
        }
#endif
      }
      w[l] -= p;
      e[l] = g;
      e[m] = 0.0;
    }
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/
int main(int argc, char *argv[])
{
  char MYFILE[] = "chsq_nu-fit_NH_SI_true_ee_test.dat";

  // Initialise and define filename of output chains
  FILE *outfile = NULL;

  outfile = fopen(MYFILE, "w");
  if (outfile == NULL)
  {
    printf("Error opening output file.\n");
    return -1;
  }

  /* Initialize libglobes */
  glbInit(argv[0]);
  // glbSelectMinimizer(GLB_MIN_POWELL);

  glbRegisterProbabilityEngine(7, /* Number of parameters */
                               &my_probability_matrix,
                               &my_set_oscillation_parameters,
                               &my_get_oscillation_parameters,
                               NULL);

  /* Initialize experiment */
  char AEDLFILE1[] = "DUNE_GLoBES.glb";
  glbInitExperiment(AEDLFILE1, &glb_experiment_list[0], &glb_num_of_exps);
  /* Intitialize output */

  double my_M_PI = 3.1415926535;

  // nu-fit
  double true_theta12 = asin(sqrt(0.307));
  double true_theta13 = asin(sqrt(0.02224));
  double true_theta23 = asin(sqrt(0.454));
  double true_sdm = 7.41e-5;
  double true_ldm = 2.505e-3;
  double true_deltacp = 232 * my_M_PI / 180;

  /* NSI */
  double true_eps_e_e = 0; // SI for true
  
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  
  glbDefineParams(true_values, true_theta12, true_theta13, true_theta23, true_deltacp, true_sdm, true_ldm);
  glbSetOscParams(true_values, true_eps_e_e, GLB_EPS_E_E);
  glbSetOscillationParameters(true_values);
  glbSetDensityParams(true_values, 1.0, GLB_ALL);
  glbSetRates();

  double m, res;

  double sigma_theta13 = (0.02397 - 0.02047) / 6;

  for (double test_eps_e_e = 0; test_eps_e_e <= 1.0; test_eps_e_e += 0.025)
  {
    for (double marg_theta13 = 0.02047; marg_theta13 <= 0.02397; marg_theta13 += (0.02397 - 0.02047) / 3)
    {
      for (double marg_theta23 = 0.411; marg_theta23 <= 0.606; marg_theta23 += (0.606 - 0.411) / 3)
      {
        for (double marg_theta12 = 0.275; marg_theta12 <= 0.344; marg_theta12 += (0.344 - 0.275) / 3)
        {
          for (double marg_deltacp = 139; marg_deltacp <= 350; marg_deltacp += 20)
          {
            for (double marg_sdm = 6.81e-5; marg_sdm <= 8.03e-5; marg_sdm += (8.03e-5 - 6.81e-5) / 3)
            {
              for (double marg_ldm = 2.426e-3; marg_ldm <= 2.586e-3; marg_ldm += (2.586e-3 - 2.426e-3) / 3)
              {
                double thetheta13 = asin(sqrt(marg_theta13));
                double thetheta23 = asin(sqrt(marg_theta23));
                double thetheta12 = asin(sqrt(marg_theta12));
                double thedeltacp = marg_deltacp * my_M_PI / 180;

                glbDefineParams(test_values, thetheta12, thetheta13, thetheta23, thedeltacp, marg_sdm, marg_ldm);
                glbSetDensityParams(test_values, 1.0, GLB_ALL);
                glbSetOscParams(test_values, test_eps_e_e, GLB_EPS_E_E);

                glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);
                res = glbChiSys(test_values, GLB_ALL, GLB_ALL);
                double prior13 = pow(((marg_theta13 - 0.02224) / sigma_theta13), 2);

                res += prior13;

                if (res < m)
                {
                  m = res;
                }
              }
            }
          }
        }
      }
    }
    fprintf(stdout, "%lf %lf \n", test_eps_e_e, m);
    fprintf(outfile, "%lf %lf \n", test_eps_e_e, m);
  }

  glbFreeParams(true_values);
  glbFreeParams(test_values);

  exit(0);
}
