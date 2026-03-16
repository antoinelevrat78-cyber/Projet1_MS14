#include "mesh.h"

int main(int argc, char* argv[])
{
  int    iTri, iVer;
  double to, ti;

  if (argc < 2) {
    printf(" usage : mesh file \n");
    return 0;
  }

  //--- read a mesh
  to        = clock();
  Mesh* Msh = msh_read(argv[1], 0);
  ti        = clock();

  if (!Msh)
    return 0;

  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- create neigbhors Q2 version
  to = clock();
  msh_neighborsQ2(Msh);
  ti = clock();
  printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- create neigbhors with hash table
  to = clock();
  msh_neighbors(Msh);
  ti = clock();
  printf("  time hash tab neigh.  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //Quality 1
  double* Qal1 = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));

  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    Qal1[iTri] = msh_tri_quality_Q1(Msh, iTri);
  }

  msh_write2dfield_Triangles("results/quality1.solb", Msh->NbrTri, Qal1);

  //Quality2 
  double* Qal2 = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));

  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    Qal2[iTri] = msh_tri_quality_Q2(Msh, iTri);
  }

  msh_write2dfield_Triangles("results/quality2.solb", Msh->NbrTri, Qal2);

  //--- TODO: compute metric field
  double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));

  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    Met[iVer][0] = 1.;
    Met[iVer][1] = 0.;
    Met[iVer][2] = 1.;
  }

  msh_write2dmetric("results/metric.solb", Msh->NbrVer, Met);

  //--- Free memory
  if (Qal1 != NULL) {
    free(Qal1);
    Qal1 = NULL;
  }

  if (Qal2 != NULL) {
    free(Qal2);
    Qal2 = NULL;
  }

  if (Met != NULL) {
    free(Met);
    Met = NULL;
  }

  return 0;
}
