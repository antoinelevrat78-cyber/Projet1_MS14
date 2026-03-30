#include "mesh.h"


void msh_square_insertion(int nPoints) {
    Mesh* Msh = msh_init();
    Msh->Dim = 2;
    Msh->NbrVerMax = nPoints + 4;
    Msh->NbrTriMax = 2 * nPoints + 2;

    Msh->Crd = calloc(Msh->NbrVerMax + 1, sizeof(double2d));
    Msh->Tri = calloc(Msh->NbrTriMax + 1, sizeof(int3d));
    Msh->TriRef = calloc(Msh->NbrTriMax + 1, sizeof(int1d));

    Msh->NbrVer = 4; // find the 4 vertices of the init rectangle
    Msh->Crd[1][0] = 0.0; 
    Msh->Crd[1][1] = 0.0;
    Msh->Crd[2][0] = 1.0; 
    Msh->Crd[2][1] = 0.0;
    Msh->Crd[3][0] = 1.0; 
    Msh->Crd[3][1] = 1.0;
    Msh->Crd[4][0] = 0.0; 
    Msh->Crd[4][1] = 1.0;

    Msh->NbrTri = 2; //init triangles
    Msh->Tri[1][0] = 1; 
    Msh->Tri[1][1] = 2; 
    Msh->Tri[1][2] = 3;
    Msh->Tri[2][0] = 1; 
    Msh->Tri[2][1] = 3; 
    Msh->Tri[2][2] = 4;

    srand(7); // rdm gen

    for (int i = 0; i < nPoints; i++) {
        double px = (double)rand() / RAND_MAX; // rdm gen in 01
        double py = (double)rand() / RAND_MAX;

        //find the triangle where the point is generated
        int targetTri = -1;
        
        for (int t = 1; t <= Msh->NbrTri; t++) { //loop on all (simple algo)
            if (point_in_tri(Msh, t, px, py)) {
                targetTri = t; //found !
                break;
            }
        }
        
        if (targetTri != -1) {
            msh_insert_and_split(Msh, targetTri, px, py);
            msh_neighbors(Msh);
            
            for (int e = 0; e < 3; e++) {
                msh_check_and_flip(Msh, Msh->NbrTri - 2, e);
                msh_neighbors(Msh);
                msh_check_and_flip(Msh, Msh->NbrTri - 1, e);
                msh_neighbors(Msh);
                msh_check_and_flip(Msh, Msh->NbrTri, e);
                msh_neighbors(Msh);
            }
            
          }
    }
    printf("Maillage généré : %d sommets, %d triangles.\n", Msh->NbrVer, Msh->NbrTri);
    //for  CPU time
    //msh_write(Msh, "edge_flip_carre.mesh");
    //msh_write_sol_quality(Msh, "edge_flip_carre.sol");
}


void msh_delaunay_square_insertion(int nPoints) {
    Mesh* Msh = msh_init();
    Msh->Dim = 2;
    Msh->NbrVerMax = nPoints + 4;
    Msh->NbrTriMax = 2 * nPoints + 2;

    Msh->Crd = calloc(Msh->NbrVerMax + 1, sizeof(double2d));
    Msh->Tri = calloc(Msh->NbrTriMax + 1, sizeof(int3d));
    Msh->TriRef = calloc(Msh->NbrTriMax + 1, sizeof(int1d));

    Msh->NbrVer = 4; // find the 4 vertices of the init rectangle
    Msh->Crd[1][0] = 0.0; 
    Msh->Crd[1][1] = 0.0;
    Msh->Crd[2][0] = 1.0; 
    Msh->Crd[2][1] = 0.0;
    Msh->Crd[3][0] = 1.0; 
    Msh->Crd[3][1] = 1.0;
    Msh->Crd[4][0] = 0.0; 
    Msh->Crd[4][1] = 1.0;

    Msh->NbrTri = 2; //init triangles
    Msh->Tri[1][0] = 1; 
    Msh->Tri[1][1] = 2; 
    Msh->Tri[1][2] = 3;
    Msh->Tri[2][0] = 1; 
    Msh->Tri[2][1] = 3; 
    Msh->Tri[2][2] = 4;

    msh_neighbors(Msh);

    srand(7); // rdm gen

    for (int i = 0; i < nPoints; i++) {
        double px = (double)rand() / RAND_MAX; // rdm gen in 01
        double py = (double)rand() / RAND_MAX;

        msh_global_sphere_criteria(Msh, px, py);
    }
    printf("Maillage généré : %d sommets, %d triangles.\n", Msh->NbrVer, Msh->NbrTri);
    //for  CPU time
    //msh_write(Msh, "delaunay_carre50.mesh");
    //msh_write_sol_quality(Msh,"delaunay_carre50.sol");
}





int main(int argc, char* argv[])
{/*
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

*/

int nbpoint = 3000;
double to, ti;

to = clock();
msh_square_insertion(nbpoint);
ti = clock();

printf("edge_flip has taken.  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

to = clock();
msh_delaunay_square_insertion(nbpoint);
ti = clock();

printf("Delaunay old has taken.  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

return 0;
}
