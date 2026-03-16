
/*---------------------------------------------------------

 Compute the spectral decomposition (eigenvalues and eigenvectors) of a 2x2 metric
   M = tR * L * R
   where R is the transformation matrix the lines of which are the eigenvectors

   Met = [ 0 1 ]
         [   2 ]

   We always have : eigVal[0] < eigVal[1]

   The eigenvectors (U,V) are in the array eigVec[4] as :
   R = [ u_x = eigVec[0]  u_y = eigVec[1] ]
       [ v_x = eigVec[2]  v_y = eigVec[3] ]

   Remark: R is a rotation matrix thus R = [ cos(a)  -sin(a) ]
                                           [ sin(a)   cos(a) ]

   => eigVec[0] should be equal to eigVec[3]

---------------------------------------------------------*/

#include "mesh.h"

/*

 Compute 2d spectral decomposition of a symmetric matrix.
 By convention we have always : |eigVal[0]| < |eigVal[1]|.

*/
int Sol_Eigen2d(double met[3], double eigVal[2], double eigVec[4])
{
  int    i;
  double nrm, sign, vecNrm, dd, tmp;
  double mat[3];

  nrm = fabs(met[0]);
  if (met[0] >= 0.)
    sign = 1.;
  else
    sign = -1;
  for (int i = 1; i < 3; ++i) {
    if (fabs(met[i]) > nrm) {
      nrm = fabs(met[i]);
      if (met[i] >= 0.)
        sign = 1.;
      else
        sign = -1;
    }
  }

  //--- a null matrix
  if (nrm < 1.e-100) {
    eigVal[0] = 0.;
    eigVal[1] = 0.;
    eigVec[0] = 1.;
    eigVec[1] = 0.;
    eigVec[2] = 0.;
    eigVec[3] = 1.;
    return 2;
  }

  nrm = sign * nrm;

  //--- compute eigenvalues
  for (i = 0; i < 3; ++i)
    mat[i] = met[i] / nrm;

  dd = (mat[0] - mat[2]) * (mat[0] - mat[2]) + 4. * mat[1] * mat[1];

  //--- Filter a diagonal matrix with two identical eigenvalues
  if (fabs(dd) < 1.e-24) { // eigVal[0] = eigVal[1] = nrm (or 1 after normalization)
    eigVal[0] = nrm;
    eigVal[1] = nrm;
    eigVec[0] = 1.;
    eigVec[1] = 0.;
    eigVec[2] = 0.;
    eigVec[3] = 1.;
    return 1;
  }
  else if (dd < 0.) {
    printf("  ## ERROR: M is always diag in IR => dd = %le must be > 0", dd);
    printf("  For metric = [%le %le %le]  normalized to [%le %le %le]\n",
        met[0], met[1], met[2], mat[0], mat[1], mat[2]);
    exit(0);
  }

  dd        = sqrt(dd);
  eigVal[0] = 0.5 * (mat[0] + mat[2] - dd);
  eigVal[1] = 0.5 * (mat[0] + mat[2] + dd);

  //--- compute eigenvectors
  eigVec[2] = -mat[1];
  eigVec[3] = mat[0] - eigVal[1];
  vecNrm    = eigVec[2] * eigVec[2] + eigVec[3] * eigVec[3];

  // In fact, for diag matrix we may only check that mat[1] ~ 0
  if (vecNrm < 1.e-30) { // => diag matrix => bad sol use the other
    eigVec[2] = eigVal[1] - mat[2];
    eigVec[3] = mat[1];
    vecNrm    = eigVec[2] * eigVec[2] + eigVec[3] * eigVec[3];
    if (vecNrm < 1.e-30) { //--- in the case we have dd = 0 ie egv1 = egv2
                           //--- thus M is the Id matrix after normalisation
      if (fabs(eigVal[0] - 1.) < 1.e-12 && fabs(eigVal[1] - 1.) < 1.e-12) {
        printf("  ## INFO Eigen2d: two identical eigen values have not been filtered !!\n");
        printf("     M = [%12.8e %12.8e %12.8e]  with eig %le %le  and dd = %le\n",
            met[0], met[1], met[2], eigVal[0], eigVal[1], dd);
        printf(" vec = %le %le or %le %le with nrm = %le or %le\n", -mat[1], mat[0] - eigVal[1],
            sqrt(mat[1] * mat[1] + (mat[0] - eigVal[1]) * (mat[0] - eigVal[1])), eigVec[2], eigVec[3], vecNrm);
        eigVal[0] = nrm;
        eigVal[1] = nrm;
        eigVec[0] = 1.;
        eigVec[1] = 0.;
        eigVec[2] = 0.;
        eigVec[3] = 1.;
        return 1;
      }
      else {
        fprintf(stderr, "\n  ## ERROR Eigen2d: A nul metric is obtained !\n");
        fprintf(stderr, "     M = [%12.8e %12.8e %12.8e]\n", mat[0], mat[1], mat[2]);
        return 0;
      }
    }
  }

  vecNrm = 1. / sqrt(vecNrm);
  eigVec[2] *= vecNrm;
  eigVec[3] *= vecNrm;

  if (fabs(eigVal[0]) < fabs(eigVal[1])) {
    eigVec[0] = eigVec[3];
    eigVec[1] = -eigVec[2];
    // eigVec[0] = -eigVec[3];
    // eigVec[1] =  eigVec[2];
  }
  else {
    eigVec[0] = eigVec[2];
    eigVec[1] = eigVec[3];
    eigVec[2] = -eigVec[1];
    eigVec[3] = eigVec[0];
    tmp       = eigVal[0];
    eigVal[0] = eigVal[1];
    eigVal[1] = tmp;
  }

  eigVal[0] *= nrm;
  eigVal[1] *= nrm;

  return 1;
}
