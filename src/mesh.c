#include "mesh.h"

int tri2edg[3][2] = { { 1, 2 }, { 2, 0 }, { 0, 1 } };

Mesh* msh_init()
{
  Mesh* Msh = malloc(sizeof(Mesh));
  if (!Msh) return NULL;

  Msh->Dim    = 0;
  Msh->NbrVer = 0;
  Msh->NbrTri = 0;
  Msh->NbrEfr = 0;
  Msh->NbrEdg = 0;

  Msh->NbrVerMax = 0;
  Msh->NbrTriMax = 0;
  Msh->NbrEfrMax = 0;
  Msh->NbrEdgMax = 0;

  Msh->Box[0] = 1.e30; // xmin
  Msh->Box[1] = -1.e30; // xmax
  Msh->Box[2] = 1.e30; // ymin
  Msh->Box[3] = -1.e30; // ymax

  //--- Data for the list of vertices
  Msh->Crd = NULL;

  //--- Data for the list of triangles
  Msh->Tri    = NULL;
  Msh->TriVoi = NULL;
  Msh->TriRef = NULL;
  Msh->TriMrk = NULL;

  //--- Data for the list of boundary edges
  Msh->Efr    = NULL;
  Msh->EfrVoi = NULL;
  Msh->EfrRef = NULL;

  //--- Data for the list of edges
  Msh->Edg = NULL;

  return Msh;
}

Mesh* msh_read(char* file, int readEfr)
{
  char   InpFil[1024];
  float  bufFlt[2];
  double bufDbl[2];
  int    i, bufTri[4], bufEfr[3];
  int    FilVer, ref;

  int fmsh = 0;

  if (!file) return NULL;

  Mesh* Msh = msh_init();

  //--- set file name
  strcpy(InpFil, file);
  if (strstr(InpFil, ".mesh")) {
    if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
      return NULL;
    }
  }
  else {
    strcat(InpFil, ".meshb");
    if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
      strcpy(InpFil, file);
      strcat(InpFil, ".mesh");
      if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
        return NULL;
      }
    }
  }

  printf(" File %s opened Dimension %d Version %d \n", InpFil, Msh->Dim, FilVer);

  Msh->NbrVer = GmfStatKwd(fmsh, GmfVertices);
  Msh->NbrTri = GmfStatKwd(fmsh, GmfTriangles);

  Msh->NbrVerMax = Msh->NbrVer;
  Msh->NbrTriMax = Msh->NbrTri;

  //--- allocate arrays
  Msh->Crd    = calloc((Msh->NbrVerMax + 1), sizeof(double3d));
  Msh->Tri    = calloc((Msh->NbrTriMax + 1), sizeof(int3d));
  Msh->TriRef = calloc((Msh->NbrTriMax + 1), sizeof(int1d));
  Msh->TriMrk = calloc((Msh->NbrTriMax + 1), sizeof(int1d));

  //--- read vertices
  GmfGotoKwd(fmsh, GmfVertices);
  if (Msh->Dim == 2) {
    if (FilVer == GmfFloat) { // read 32 bits float
      for (i = 1; i <= Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &ref);
        Msh->Crd[i][0] = (double)bufFlt[0];
        Msh->Crd[i][1] = (double)bufFlt[1];
      }
    }
    else { // read 64 bits float
      for (i = 1; i <= Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
        Msh->Crd[i][0] = bufDbl[0];
        Msh->Crd[i][1] = bufDbl[1];
      }
    }
  }
  else {
    fprintf(stderr, "  ## ERROR: 3D is not implemented\n");
    exit(1);
  }

  //--- read triangles
  GmfGotoKwd(fmsh, GmfTriangles);
  for (i = 1; i <= Msh->NbrTri; ++i) {
    GmfGetLin(fmsh, GmfTriangles, &bufTri[0], &bufTri[1], &bufTri[2], &bufTri[3]);
    Msh->Tri[i][0] = bufTri[0];
    Msh->Tri[i][1] = bufTri[1];
    Msh->Tri[i][2] = bufTri[2];
    Msh->TriRef[i] = bufTri[3];
  }

  //--- read boundary edges
  if (readEfr == 1) {
    Msh->NbrEfr    = GmfStatKwd(fmsh, GmfEdges);
    Msh->NbrEfrMax = Msh->NbrEfr;

    Msh->Efr    = calloc((Msh->NbrEfrMax + 1), sizeof(int2d));
    Msh->EfrRef = calloc((Msh->NbrEfrMax + 1), sizeof(int1d));

    GmfGotoKwd(fmsh, GmfEdges);
    for (i = 1; i <= Msh->NbrEfr; ++i) {
      GmfGetLin(fmsh, GmfEdges, &bufEfr[0], &bufEfr[1], &bufEfr[2]);
      Msh->Efr[i][0] = bufEfr[0];
      Msh->Efr[i][1] = bufEfr[1];
      Msh->EfrRef[i] = bufEfr[2];
    }
  }

  GmfCloseMesh(fmsh);

  return Msh;
}

double* sol_read(char* file, int mshDim, int mshNbrSol)
{
  char   InpFil[1024];
  int    FilVer, SolTyp, NbrTyp, SolSiz, TypTab[GmfMaxTyp];
  float  bufFlt;
  double bufDbl;
  int    i, dim, nbrSol;

  int fsol = 0;

  if (!file) return NULL;

  double* sol = NULL;

  //--- set file name
  strcpy(InpFil, file);
  if (strstr(InpFil, ".sol")) {
    if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
      return NULL;
    }
  }
  else {
    strcat(InpFil, ".solb");
    if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
      strcpy(InpFil, file);
      strcat(InpFil, ".sol");
      if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
        return NULL;
      }
    }
  }

  printf(" File %s opened Dimension %d Version %d \n", InpFil, dim, FilVer);

  SolTyp = GmfSolAtVertices; // read only sol at vertices
  nbrSol = GmfStatKwd(fsol, SolTyp, &NbrTyp, &SolSiz, TypTab);

  if (nbrSol == 0) {
    printf("  ## WARNING: No SolAtVertices in the solution file !\n");
    return NULL;
  }
  if (dim != mshDim) {
    printf("  ## WARNING: WRONG DIMENSION NUMBER. IGNORED\n");
    return NULL;
  }
  if (nbrSol != mshNbrSol) {
    printf("  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    return NULL;
  }
  if (NbrTyp != 1) {
    printf("  ## WARNING: WRONG FIELD NUMBER. IGNORED\n");
    return NULL;
  }
  if (TypTab[0] != GmfSca) {
    printf("  ## WARNING: WRONG FIELD TYPE. IGNORED\n");
    return NULL;
  }

  sol = (double*)calloc(nbrSol + 1, sizeof(double));

  GmfGotoKwd(fsol, SolTyp);

  for (i = 1; i <= nbrSol; ++i) {
    if (FilVer == GmfFloat) {
      GmfGetLin(fsol, SolTyp, &bufFlt);
      sol[i] = (double)bufFlt;
    }
    else {
      GmfGetLin(fsol, SolTyp, &bufDbl);
      sol[i] = bufDbl;
    }
  }

  if (!GmfCloseMesh(fsol)) {
    fprintf(stderr, "  ## ERROR: Cannot close solution file %s ! \n", InpFil);
    // myexit(1);
  }

  return sol;
}

//Done with a different approach
int msh_boundingbox(Mesh* Msh)
{
  int1d iVer;

  //--- compute bounding box
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    // TODO: Set Msh->Box
  }

  return 1;
}

int msh_write(Mesh* Msh, char* file)
{
  int iVer, iTri, iEfr;
  int FilVer = 2;

  if (!Msh) return 0;
  if (!file) return 0;

  int fmsh = GmfOpenMesh(file, GmfWrite, FilVer, Msh->Dim);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  GmfSetKwd(fmsh, GmfVertices, Msh->NbrVer);
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++)
    GmfSetLin(fmsh, GmfVertices, Msh->Crd[iVer][0], Msh->Crd[iVer][1], 0);

  GmfSetKwd(fmsh, GmfTriangles, Msh->NbrTri);
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++)
    GmfSetLin(fmsh, GmfTriangles, Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2], Msh->TriRef[iTri]);

  if (Msh->NbrEfr > 0) {
    GmfSetKwd(fmsh, GmfEdges, Msh->NbrEfr);
    for (iEfr = 1; iEfr <= Msh->NbrEfr; iEfr++)
      GmfSetLin(fmsh, GmfEdges, Msh->Efr[iEfr][0], Msh->Efr[iEfr][1], Msh->EfrRef[iEfr]);
  }

  GmfCloseMesh(fmsh);

  return 1;
}

void msh_write_sol_quality(Mesh* Msh, const char* filename) {
    FILE* f = fopen(filename, "w");
    if (!f) return;

    fprintf(f, "MeshVersionFormatted 2\n");
    fprintf(f, "Dimension 2\n\n");
    fprintf(f, "SolAtTriangles\n");
    fprintf(f, "%d\n", Msh->NbrTri);
    fprintf(f, "1 1\n"); 

    for (int t = 1; t <= Msh->NbrTri; t++) {
        double q = msh_tri_quality_Q1(Msh, t);
        fprintf(f, "%f\n", q);
    }

    fprintf(f, "\nEnd\n");
    fclose(f);
    printf("sol generted : %s\n", filename);
}
// Exercice 1:

double msh_tri_quality_Q1(Mesh* Msh, int iTri) {
    double x[3], y[3], l2[3], area;
    int i;
    const double alpha1 = sqrt(3.0) / 12.0;

    for (i = 0; i < 3; i++) {
        int iVer = Msh->Tri[iTri][i];
        x[i] = Msh->Crd[iVer][0];
        y[i] = Msh->Crd[iVer][1];
    }

    // lenght of the triangle's sides 
    l2[0] = pow(x[1]-x[0], 2) + pow(y[1]-y[0], 2);
    l2[1] = pow(x[2]-x[1], 2) + pow(y[2]-y[1], 2);
    l2[2] = pow(x[0]-x[2], 2) + pow(y[0]-y[2], 2);

    // Area of the triangle
    double dx1 = x[1] - x[0];
    double dy1 = y[1] - y[0];
    double dx2 = x[2] - x[0];
    double dy2 = y[2] - y[0];

    area = 0.5 * fabs(dx1 * dy2 - dx2 * dy1);

    return alpha1 * (l2[0] + l2[1] + l2[2]) / area; 
}

double msh_tri_quality_Q2(Mesh* Msh, int iTri) {
    double x[3], y[3], l2[3], hmax = 0, area, p, rho;
    int i;
    const double alpha2 = sqrt(3.0) / 6.0; 

    for (i = 0; i < 3; i++) {
        int iVer = Msh->Tri[iTri][i];
        x[i] = Msh->Crd[iVer][0];
        y[i] = Msh->Crd[iVer][1];
    }

    l2[0] = sqrt(pow(x[1]-x[0], 2) + pow(y[1]-y[0], 2));
    l2[1] = sqrt(pow(x[2]-x[1], 2) + pow(y[2]-y[1], 2));
    l2[2] = sqrt(pow(x[0]-x[2], 2) + pow(y[0]-y[2], 2));

    for (i = 0; i < 3; i++) if (l2[i] > hmax) hmax = l2[i];

    // Area of the triangle
    double dx1 = x[1] - x[0];
    double dy1 = y[1] - y[0];
    double dx2 = x[2] - x[0];
    double dy2 = y[2] - y[0];

    area = 0.5 * dx1 * dy2 - dx2 * dy1;
    p = (l2[0] + l2[1] + l2[2]); 
    rho = 2 * area / p; // radius of the inscribded's cercle

    return alpha2 * hmax / rho; 
}

// End of exercice 1


int msh_neighborsQ2(Mesh* Msh)
{
  int iTri, iEdg, jTri, jEdg, iVer1, iVer2, jVer1, jVer2;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  //--- Compute the neighbors using a quadratic-complexity algorithm
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      //--- find the Tri different from iTri that has iVer1, iVer2 as vertices
      for (jTri = 1; jTri <= Msh->NbrTri; jTri++) {
        if (iTri == jTri)
          continue;

        for (jEdg = 0; jEdg < 3; jEdg++) {
          jVer1 = Msh->Tri[jTri][tri2edg[jEdg][0]];
          jVer2 = Msh->Tri[jTri][tri2edg[jEdg][1]];

          // TODO: compare the 4 points
          //       set the neighbors Msh->TriVoi if both edges match
        }
      }
    }
  }

  return 1;
}

int msh_neighbors(Mesh* Msh)
{
  int iTri, iEdg, iVer1, iVer2, iObj, iVoi;
  int nbrEfr = 0, maxColl = 0, occupiedHeads = 0;
  HashTable* hsh = NULL;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  // Initialize hash table: SizHead ~ NbrTri, NbrMaxObj = 3 * NbrTri 
  hsh = hash_init(Msh->NbrTri, 3 * Msh->NbrTri);

  // Fill hash table with the triangle edges
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];
      hash_add(hsh, iVer1, iVer2, iTri);
    }
  }

  //Neighbors using the hash table
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      iObj = hash_find(hsh, iVer1, iVer2);
      if (iObj != 0) {
        // Find who is the neighbor
        iVoi = (hsh->LstObj[iObj][2] == iTri) ? hsh->LstObj[iObj][3] : hsh->LstObj[iObj][2];
        Msh->TriVoi[iTri][iEdg] = iVoi;
      }
    }
  }
  
  /*
  // Here we add a part to calculate the stats of different Keys 
  // counting the sides
  for (iObj = 1; iObj <= hsh->NbrObj; iObj++) {
    if (hsh->LstObj[iObj][3] == 0) nbrEfr++;
  }

  // Collision analysis
  for (int iPos = 0; iPos < hsh->SizHead; iPos++) {
    int count = 0;
    int currObj = hsh->Head[iPos];
    if (currObj != 0) {
      occupiedHeads++;
      while (currObj != 0) {
        count++;
        currObj = hsh->LstObj[currObj][4];
      }
      if (count > maxColl) maxColl = count;
    }
  }

  // Results
  printf("\n--- Mesh Statistics ---\n");
  printf("Total unique edges:    %d\n", hsh->NbrObj);
  printf("Boundary edges:        %d\n", nbrEfr);
  printf("Max collisions:        %d\n", maxColl);
  if (occupiedHeads > 0)
    printf("Average collisions:    %.2f\n", (float)hsh->NbrObj / occupiedHeads);
  printf("-----------------------\n");
  */
  return 1;
}



HashTable* hash_init(int SizHead, int NbrMaxObj)
{
  HashTable* hsh = NULL;

  // allocate hash table
  hsh = (HashTable*)malloc(sizeof(HashTable));

  // initialize hash table
  hsh->SizHead = SizHead;
  hsh->NbrMaxObj = NbrMaxObj;
  hsh->NbrObj = 0;

  // allocate Head, LstObj
  hsh->Head = (int*)calloc(SizHead, sizeof(int));
  hsh->LstObj = (int5d*)calloc(NbrMaxObj + 1, sizeof(int5d));

  return hsh;
}

int hash_find(HashTable* hsh, int iVer1, int iVer2)
{
  //iVer1 and iVer2 are points 
  //we order them 
  int iPos, iObj;
  int vMin = (iVer1 < iVer2) ? iVer1 : iVer2;
  int vMax = (iVer1 < iVer2) ? iVer2 : iVer1;

  // Key declaration
  //Key 1
  //iPos = (vMin + vMax) % hsh->SizHead; 
  //Key 2
  //iPos = vMin % hsh->SizHead;
  //Key 3:
  iPos = (vMin + vMax * 1000) % hsh->SizHead; 


  // Traverse the linked list for this key 
  iObj = hsh->Head[iPos]; // key at first position
  while (iObj != 0) {
    if (hsh->LstObj[iObj][0] == vMin && hsh->LstObj[iObj][1] == vMax) {
      return iObj; // Edge found
    }
    iObj = hsh->LstObj[iObj][4]; // Move to next object in collision
  }

  return 0;
}

int hash_add(HashTable* hsh, int iVer1, int iVer2, int iTri)
{
  int iPos, iObj;
  int vMin = (iVer1 < iVer2) ? iVer1 : iVer2;
  int vMax = (iVer1 < iVer2) ? iVer2 : iVer1;

  // Check if edge already exists
  iObj = hash_find(hsh, vMin, vMax);

  if (iObj != 0) {
    // Edge exists, add the second triangle that shares this edge
    hsh->LstObj[iObj][3] = iTri;
    return iObj;
  }

  // Add new edge if enough space
  if (hsh->NbrObj >= hsh->NbrMaxObj) return 0;

  hsh->NbrObj++;
  iObj = hsh->NbrObj; //This becomes the new index for the edge
  //Key declaration
  //Key 1:
  //iPos = (vMin + vMax) % hsh->SizHead;
  //Key 2:
  //iPos = vMin % hsh->SizHead;
  //Key 3:
  iPos = (vMin + vMax * 1000) % hsh->SizHead; 

  // Store edge data 
  hsh->LstObj[iObj][0] = vMin;
  hsh->LstObj[iObj][1] = vMax;
  hsh->LstObj[iObj][2] = iTri; // First triangle
  hsh->LstObj[iObj][3] = 0;    // Second triangle (unknown yet)
  
  // Update linked list (insertion at head)
  hsh->LstObj[iObj][4] = hsh->Head[iPos];
  hsh->Head[iPos] = iObj;

  return iObj;
}

int hash_suppr(HashTable* hsh, int iVer1, int iVer2, int iTri)
{

  // to be implemented
  // will be done later i guess
  // ===> suppress this entry in the hash tab

  return 0;
}

int msh_write2dfield_Vertices(char* file, int nfield, double* field)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSca;

  GmfSetKwd(fmsh, GmfSolAtVertices, nfield, 1, sizfld);

  for (iVer = 1; iVer <= nfield; iVer++)
    GmfSetLin(fmsh, GmfSolAtVertices, &field[iVer]);

  GmfCloseMesh(fmsh);

  return 1;
}

int msh_write2dfield_Triangles(char* file, int nfield, double* field)
{
  int iTri;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSca;

  GmfSetKwd(fmsh, GmfSolAtTriangles, nfield, 1, sizfld);

  for (iTri = 1; iTri <= nfield; iTri++)
    GmfSetLin(fmsh, GmfSolAtTriangles, &field[iTri]);

  GmfCloseMesh(fmsh);

  return 1;
}

int msh_write2dmetric(char* file, int nmetric, double3d* metric)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSymMat;

  GmfSetKwd(fmsh, GmfSolAtVertices, nmetric, 1, sizfld);

  for (iVer = 1; iVer <= nmetric; iVer++)
    GmfSetLin(fmsh, GmfSolAtVertices, &metric[iVer][0], &metric[iVer][1], &metric[iVer][2]);

  GmfCloseMesh(fmsh);

  return 1;
}


// Project 2: first implementation of edge flip 
// interts a point and creates 3 triangles
int msh_insert_and_split(Mesh* Msh, int iTri, double x, double y) {
    if (Msh->NbrTri + 2 > Msh->NbrTriMax || Msh->NbrVer + 1 > Msh->NbrVerMax) {
        printf("Error :mesh max capacity .\n");
        return 0;
    }

    // we create a point
    int vP = ++Msh->NbrVer;
    Msh->Crd[vP][0] = x;
    Msh->Crd[vP][1] = y;

    // 3 summits of the triangle
    int v1 = Msh->Tri[iTri][0];
    int v2 = Msh->Tri[iTri][1];
    int v3 = Msh->Tri[iTri][2];
    

    // create 2 other triangles
    int iTri2 = ++Msh->NbrTri;
    int iTri3 = ++Msh->NbrTri;

    // trigonometric
    Msh->Tri[iTri][0] = v1; Msh->Tri[iTri][1] = v2; Msh->Tri[iTri][2] = vP;
    Msh->Tri[iTri2][0] = v2; Msh->Tri[iTri2][1] = v3; Msh->Tri[iTri2][2] = vP;
    Msh->Tri[iTri3][0] = v3; Msh->Tri[iTri3][1] = v1; Msh->Tri[iTri3][2] = vP;

    //Rq : need to actualize the neibors
    return 1;
}

// checks and flips if quality Q1 is better for target edge
int msh_check_and_flip(Mesh* Msh, int iTri, int iEdg) {
    int jTri = Msh->TriVoi[iTri][iEdg]; // neibor
    if (jTri <= 0) return 0; // no neibours

    // find the vertices of the quadrilatère
    int vA = Msh->Tri[iTri][iEdg]; // opposed vertice
    int vB = Msh->Tri[iTri][(iEdg + 1) % 3];
    int vC = Msh->Tri[iTri][(iEdg + 2) % 3];

    //opposed vertice of the triangle neibor
    int vD = -1;
    for (int k = 0; k < 3; k++) {
        int v = Msh->Tri[jTri][k];
        if (v != vB && v != vC) {
            vD = v;
            break;
        }
    }

    double area1 = tri_area(Msh->Crd[vA], Msh->Crd[vB], Msh->Crd[vD]);
    double area2 = tri_area(Msh->Crd[vA], Msh->Crd[vD], Msh->Crd[vC]);
    
    if (area1 <= 1e-12 || area2 <= 1e-12) return 0;

    // find the max of the 2 qualities
    double q_old = fmax(msh_tri_quality_Q1(Msh, iTri), msh_tri_quality_Q1(Msh, jTri));

    // backup before edge flip
    int3d backup_i = {Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2]};
    int3d backup_j = {Msh->Tri[jTri][0], Msh->Tri[jTri][1], Msh->Tri[jTri][2]};

    // edge flip for test
    Msh->Tri[iTri][0] = vA; Msh->Tri[iTri][1] = vB; Msh->Tri[iTri][2] = vD;
    Msh->Tri[jTri][0] = vA; Msh->Tri[jTri][1] = vD; Msh->Tri[jTri][2] = vC;

    double q_new = fmax(msh_tri_quality_Q1(Msh, iTri), msh_tri_quality_Q1(Msh, jTri));

    if (q_new < q_old) { 
        //validated, we recalculate the neibours
        return 1; 
    } else {
        // get back to prior config
        for(int k=0; k<3; k++) {
            Msh->Tri[iTri][k] = backup_i[k];
            Msh->Tri[jTri][k] = backup_j[k];
        }
        return 0;
    }
}

// finding if point is in triangle 
double tri_area(double P[2], double P1[2], double P2[2]) {
    return (P1[0] - P[0]) * (P2[1] - P[1]) - (P1[1] - P[1]) * (P2[0] - P[0]);
}

//if 3 areas are positive, point is in triangle
int point_in_tri(Mesh* Msh, int iTri, double x, double y) {
    double P[2] = {x, y};
    double* A = Msh->Crd[Msh->Tri[iTri][0]];
    double* B = Msh->Crd[Msh->Tri[iTri][1]];
    double* C = Msh->Crd[Msh->Tri[iTri][2]];


    double s1 = tri_area(P, A, B);
    double s2 = tri_area(P, B, C);
    double s3 = tri_area(P, C, A);

    
    if (s1 >= 0 && s2 >= 0 && s3 >= 0) return 1;
    

    return 0;
}


//faster algo for finding using edges and barycentric coordinates
//with safety feature so that doesnt turn in circles
int msh_locate_point(Mesh* Msh, double x, double y, int iTriStart) {
    int iTri = iTriStart;
    double P[2] = {x, y};
    int iter = 0, maxIter = Msh->NbrTri;
    int edges_neg[3], n_neg;

    while (iTri > 0 && iter < maxIter) {
        iter++;
        n_neg = 0;
        
        for (int i = 0; i < 3; i++) { //test on the 3 vertices
            int v1 = Msh->Tri[iTri][(i + 1) % 3];
            int v2 = Msh->Tri[iTri][(i + 2) % 3];
            if (tri_area(P, Msh->Crd[v1], Msh->Crd[v2]) < 0) {
                edges_neg[n_neg++] = i; //i is the edgewe follow 
            }
        }

        if (n_neg == 0) return iTri; 
        int iEdgeChoice = edges_neg[rand() % n_neg]; // selects a rdm edge to avoid cercles
        
        //go to the neibour
        int nextTri = Msh->TriVoi[iTri][iEdgeChoice];
        if (nextTri <= 0) return iTri; // if on fronteer
        iTri = nextTri;
    }
    
    return iTri;
}


// Returns 1 if point P is in iTri else 0
int msh_in_circle(Mesh* Msh, int iTri, double px, double py){
    double* p1 = Msh->Crd[Msh->Tri[iTri][0]];
    double* p2 = Msh->Crd[Msh->Tri[iTri][1]];
    double* p3 = Msh->Crd[Msh->Tri[iTri][2]];
    double x1 = p1[0], y1 = p1[1];
    double x2 = p2[0], y2 = p2[1];
    double x3 = p3[0], y3 = p3[1];

    double D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

    double x1sq = x1*x1 + y1*y1;
    double x2sq = x2*x2 + y2*y2;
    double x3sq = x3*x3 + y3*y3;
    //ox,oy is the center of the inscribded circle
    double ox = (x1sq * (y2 - y3) + x2sq * (y3 - y1) + x3sq * (y1 - y2)) / D;
    double oy = (x1sq * (x3 - x2) + x2sq * (x1 - x3) + x3sq * (x2 - x1)) / D;

    double r2 = (ox - x1)*(ox - x1) + (oy - y1)*(oy - y1); //r2 is the squared radius of this circle
    double distP2 = (ox - px)*(ox - px) + (oy - py)*(oy - py); // distance sqrd btween 

    if (distP2 <= r2) {
        return 1; // is in circle
    }
    return 0; // isnt in circle
}

//tests all triangles with msh_in_circle
int msh_global_sphere_criteria(Mesh* Msh, double px, double py) {
    int Max_Cavity = 100; // assume the cavity is at most 100 triangles (simpler for now)
    int cavityTri[Max_Cavity]; 
    int nbCav = 0;

    /* 
    //old technique
    for (int t = 1; t <= Msh->NbrTri; t++) { // for all triangles in mesh
        if (msh_in_circle(Msh, t, px, py) == 1) {
            if (nbCav < Max_Cavity) cavityTri[nbCav++] = t; // we add the triangle
        }
    }*/
    int iTriStart = msh_locate_point(Msh, px, py, 1);
    cavityTri[nbCav++] = iTriStart;

    //other option to go faster
    for (int i = 0; i < nbCav; i++) {
        int t = cavityTri[i];
        
        for (int e = 0; e < 3; e++) {
            int neigh = Msh->TriVoi[t][e]; //get 3 neibours of startTri

            // is neibour allready in cavity?
            int alrdyin = 0;
            for (int j = 0; j < nbCav; j++) {
                if (cavityTri[j] == neigh) {
                    alrdyin = 1;
                    break;
                }
            }

            //if not we test the circle
            if (!alrdyin) {
                if (msh_in_circle(Msh, neigh, px, py) == 1) {
                    if (nbCav < 100) {
                        cavityTri[nbCav++] = neigh; //then we continue for other neibours
                    }
                }
            }
        }
    }

    //find the vertices to the fronteer (in cavity) edges
    int edges[300][2]; 
    int nbEdges = 0;

    for (int i = 0; i < nbCav; i++) {
        int t = cavityTri[i];
        for (int edg = 0; edg < 3; edg++) {
            int neigh = Msh->TriVoi[t][edg];  //we get the neibour 
            int neighInCavity = 0;
            for (int j = 0; j < nbCav; j++) { // test for all triangle in cavity if it has neibour
                if (neigh == cavityTri[j]) {
                    neighInCavity = 1;
                    break;
                }
            }

            if (neighInCavity == 0) { // if neibour is not in cavity we add to edges
                edges[nbEdges][0] = Msh->Tri[t][(edg + 1) % 3];
                edges[nbEdges][1] = Msh->Tri[t][(edg + 2) % 3];
                nbEdges++;
            }
        }
    }

    int vP = ++Msh->NbrVer; //we add the new point P
    Msh->Crd[vP][0] = px;
    Msh->Crd[vP][1] = py;
    //we try to reuse the old triangles spots and recreate new one if needed
    for (int i = 0; i < nbEdges; i++) {
        int targetTri;
        if (i < nbCav) {
            targetTri = cavityTri[i]; //old
        } else {
            targetTri = ++Msh->NbrTri; //new
        }
        
        Msh->Tri[targetTri][0] = edges[i][0];
        Msh->Tri[targetTri][1] = edges[i][1];
        Msh->Tri[targetTri][2] = vP;
        Msh->TriRef[targetTri] = 1;
    }

    msh_neighbors(Msh); 

    return 1;
}

