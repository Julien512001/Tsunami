# include "tsunami.h"


void triangleMap(int index, int *map)
{
    for (int i = 0; i < 3; i++) {
      map[i] = 3*index + i;
    }
}

void edgesMap(femTsunamiProblem *myProblem, int index, int map[2][2])
{
    int i,j,k;
    
    for (j=0; j < 2; ++j) {
        int node = myProblem->edges->edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = myProblem->edges->edges[index].elem[k];
            map[k][j] = (myProblem->mesh->nElem)*3;
            if (elem >= 0) {
                for (i=0; i < 3; i++) {
                    if (myProblem->mesh->elem[elem*3 + i] == node) {
                        map[k][j] = elem*3 + i;  }}}}}
}

femMesh *meshCreate(const char *meshFileName)
{
    femMesh *mesh = malloc(sizeof(femMesh));
    int i, j, nNode, nElem, trash;
    double dtrash;
    
    FILE* file = fopen(meshFileName, "r");
    if (file == NULL) {
    	  printf("Error : cannot open mesh file :-) \n");
        exit(0); 
    }
    fscanf(file, "Number of nodes %d \n", &nNode);
    double *X = (double *) malloc(sizeof(double)*nNode);
    double *Y = (double *) malloc(sizeof(double)*nNode);
    double *H = (double *) malloc(sizeof(double)*nNode);
    for (i = 0; i < nNode; i++) {
        fscanf(file, "%d : %le %le %le \n", &trash, &X[i], &Y[i], &H[i]);
    }

    fscanf(file, "Number of triangles %d \n", &nElem);
    int *elem = (int *) malloc(sizeof(int)*3*nElem);
    for (i = 0; i < nElem; i++) {
        fscanf(file, "%d : %d %d %d \n", &trash, &elem[3*i], &elem[3*i+1], &elem[3*i+2]);
    }

    femMesh *mesh;
    mesh->elem = elem;
    mesh->H = H;
    mesh->X = X;
    mesh->Y = Y;
    mesh->nElem = nElem;
    mesh->nLocalNode = 3;
    mesh->nNode = nNode;
}

void freeMesh(femMesh *mesh)
{
    free(mesh->X);
    free(mesh->Y);
    free(mesh->elem);
    free(mesh->H);
}


femTsunamiProblem *problemCreate(const char *meshFileName)
{
    femTsunamiProblem *theProblem = malloc(sizeof(femTsunamiProblem));

    int i,j,nNode,nElem, nEdges, trash;
    double dtrash;
    
    double BathMax = 9368;

    FILE* file = fopen(meshFileName, "r");
    if (file == NULL) {
    	  printf("Error : cannot open mesh file :-) \n");
        exit(0); 
    }
    femMesh = 

    fscanf(file, "Number of edges %d \n", &nEdges);
    
    femEdge *edge = malloc(sizeof(edge)*nEdges);
    for (i = 0; i < nEdges; i++) {
        fscanf(file, "%d : %d %d : %d %d \n", &trash, &edge->node[0], &edge->node[1], &edge->elem[0], &edge->elem[1]);
    }
    femEdges *edges = malloc(sizeof(femEdges));
    edges->edges = edge;
    edges->nEdge = nEdges;
        
    double *E  = malloc(sizeof(double)*nElem*3);
    for (i = 0; i < nElem; i++)
        for (j = 0; j < 3; j++)
            E[i*3+j] = H[elem[i*3+j]]/(10*BathMax);
    double *U = malloc(sizeof(double)*nElem*3);
    double *V = malloc(sizeof(double)*nElem*3);


    theProblem->E = E;
    theProblem->U = U;
    theProblem->V = V;
    theProblem->mesh = mesh;
    theProblem->edges = edges;
}

void freeStruct(femTsunamiProblem *myProblem)
{
    
}

void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 
    //
    // A modifer :-)
    // Ici, on se contente d'initialiser le champs d'élévation avec la bathymétrie
    // pour avoir un joli plot de départ !
    //
    // On divise la bathymétrie par un facteur (10*BathMax) pour que l'échelle de couleur fonctionne
    // bien dans la fonction d'animation fournie....
    //
    femTsunamiProblem *myProblem = problemCreate(meshFileName);

    
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
    int nElem = myProblem->mesh->nElem;
    tsunamiWriteFile(baseResultName,0,E,U,V,nElem,3);
    
    freeMesh(myProblem->mesh);

}








