# include "tsunami.h"

# define nTri 3
# define nEdg 2

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

static double phi2[nTri][nTri] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
static double phi1[nEdg][nEdg] = {{0.0, 0.0},{0.0, 0.0}};

static const double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}};

double interpolate(const double *phi, double *U, int *map, int n)
{
git    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
}

void phi2f(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void femTsunamiTriangleMap(int index, int map[3])
{
    int j;
    for (j=0; j < 3; ++j) 
        map[j] = index*3 + j; 
}

void femTsumnamiIntegralsElements(double* X, double* Y, double* E, double* U, double* FE, double* FU, double* FV,double* V, double* bath,int* Elem, int nElem){
       
    static double  xLoc[3],yLoc[3],dphidx[3],dphidy[3], phi[3]; //mettre en global pour optimiser?
    double  xsi,eta,weight,jac,theta, jacSphere;
    double  y,x,z_R,e,u,v,h,f;
    int     i,j,k,elem,mapElem[3];

    for (elem=0; elem < nElem; elem++) {
        femTsunamiTriangleMap(elem,mapElem); //mettre ca ds loop suivante?
        /*
        int j;
        for (j=0; j < 3; ++j) 
            map[j] = index*3 + j; 
        */
        int *mapCoord = &(Elem[elem*nTri]);
        for (j=0; j < nTri; ++j) {
        	  xLoc[j] = X[mapCoord[j]];
        	  yLoc[j] = Y[mapCoord[j]];
              mapElem[j] = elem*nTri+j; //celci c'est femShallowTriangleMap(elem,mapElem);
        }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        dphidx[0] = (yLoc[1] - yLoc[2])/jac;
        dphidx[1] = (yLoc[2] - yLoc[0])/jac;
        dphidx[2] = (yLoc[0] - yLoc[1])/jac;
        dphidy[0] = (xLoc[2] - xLoc[1])/jac;
        dphidy[1] = (xLoc[0] - xLoc[2])/jac;
        dphidy[2] = (xLoc[1] - xLoc[0])/jac;        
        for (k=0; k < nTri; k++) {
            xsi = gaussTriangleXsi[k];
            eta = gaussTriangleEta[k];
            weight = gaussTriangleWeight[k];

            /*
            //can optimize loops for interpolations
            y=0; x=0; h=0; e=0; u=0; v=0;
            for (i=0; i < nTri; i++){
                y += phi2[k][i]*Y[mapCoord[i]];
                x += phi2[k][i]*X[mapCoord[i]];
                h += phi2[k][i]*bath[mapCoord[i]];

                e += phi2[k][i]*E[mapElem[i]];
                u += phi2[k][i]*U[mapElem[i]];
                v += phi2[k][i]*V[mapElem[i]];
            }
            */
            phi2f(xsi, eta, phi);
            y = interpolate(phi,Y,mapCoord,nTri);
            x = interpolate(phi,X,mapCoord,nTri);
            h = interpolate(phi,bath,mapCoord,nTri);

            e = interpolate(phi,E,mapElem,nTri);
            u = interpolate(phi,U,mapElem,nTri);
            v = interpolate(phi,V,mapElem,nTri);
            /*
            y = interpolate(phi2[k],Y,mapCoord,nTri);
            x = interpolate(phi2[k],X,mapCoord,nTri);
            h = interpolate(phi2[k],bath,mapCoord,nTri);

            e = interpolate(phi2[k],E,mapElem,nTri);
            u = interpolate(phi2[k],U,mapElem,nTri);
            v = interpolate(phi2[k],V,mapElem,nTri);
            */
            z_R = (4*R*R-x*x-y*y)/(4*R*R+x*x+y*y);
            //theta = asin(z_R);
            //f = 2*Omega*sin(theta);
            f = 2*Omega*z_R;  //moins d'operations comme ca
            jacSphere = (4*R*R+x*x+y*y)/(4*R*R);

            for (i=0; i < nTri; i++) {
                //FE[mapElem[i]] += ((dphidx[i]*h*u+dphidy[i]*h*v)*jacSphere+phi2[k][i]*h*(x*u+y*v)/(R*R))*jac*weight;
                //FU[mapElem[i]] += (phi2[k][i]*(f*v-Gamma*u)+dphidx[i]*g*e*jacSphere+phi2[k][i]*g*x*e/(2*R*R))*jac*weight;
                //FV[mapElem[i]] += (phi2[k][i]*(-f*u-Gamma*v)+dphidy[i]*g*e*jacSphere+phi2[k][i]*g*y*e/(2*R*R))*jac*weight;
                FE[mapElem[i]] += ((dphidx[i]*h*u+dphidy[i]*h*v)*jacSphere+phi[i]*h*(x*u+y*v)/(R*R))*jac*weight;
                FU[mapElem[i]] += (phi[i]*(f*v-Gamma*u)+dphidx[i]*g*e*jacSphere+phi[i]*g*x*e/(2*R*R))*jac*weight;
                FV[mapElem[i]] += (phi[i]*(-f*u-Gamma*v)+dphidy[i]*g*e*jacSphere+phi[i]*g*y*e/(2*R*R))*jac*weight;
            }
        }
    }
}

void femTsunamiEdgeMap(int index, int map[2][2], femEdge* edges, int nElem, int* Elem)
{
    int i,j,k;   
    for (j=0; j < 2; ++j) {
        int node = edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = edges[index].elem[k];
            map[k][j] = nElem*3;
            if (elem >= 0) {
                for (i=0; i < 3; i++) {
                    if (Elem[elem*3 + i] == node) {
                        map[k][j] = elem*3 + i;  }}}}}
}

void phi1f(double xsi, double *phi) 
{   
    phi[0] = (1.0 - xsi)/2.0;
    phi[1] = (1.0 + xsi)/2.0;   
    
}

void femTsumnamiIntegralsEdges(double* X, double* Y, double* E, double* U, double* V,double* FE, double* FU, double* FV, double* bath, int nElem, int* Elem ,int nEdge, femEdge* edges)
{
    
    static double  xEdge[2],yEdge[2],phiEdge[2];
    double  xsi,weight,jac, jacSphere;
    double  eL,eR,uL,uR,vL,vR,unL,unR;
    double  qe,qu,qv,x,y,h;
    int     i,j,k,edge,mapEdge[2][2];
    int     sizeGlo = nElem * nTri + 1; 

    for (edge=0; edge < nEdge; edge++) {
        femTsunamiEdgeMap(edge,mapEdge,edges,nElem,Elem);
        for (j=0; j < 2; ++j) {
        	  int node = edges[edge].node[j];
        	  xEdge[j] = X[node];
        	  yEdge[j] = Y[node]; }
        
        int boundary = (mapEdge[1][0] == sizeGlo-1);
        
        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        //printf("%f\n", norm);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        jac = norm / 2.0;
        for (k=0; k < nEdg; k++) {
            xsi = gaussEdgeXsi[k];
            weight = gaussEdgeWeight[k]; 
            phi1f(xsi, phiEdge);    
            //femDiscretePhi1(theSpace,xsi,phiEdge);           
            eL = interpolate(phiEdge,E,mapEdge[0],2);
            eR = boundary ? eL : interpolate(phiEdge,E,mapEdge[1],2);
            uL = interpolate(phiEdge,U,mapEdge[0],2);
            uR = interpolate(phiEdge,U,mapEdge[1],2);
            vL = interpolate(phiEdge,V,mapEdge[0],2);
            vR = interpolate(phiEdge,V,mapEdge[1],2);
            x = interpolate(phiEdge,X,edges[edge].node,2);
            y = interpolate(phiEdge,Y,edges[edge].node,2);
            h = interpolate(phiEdge,bath,edges[edge].node,2);
            //h = 100;
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            //qe =  0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) );
            //qu =  0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );
            //qv =  0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );
            qe = (eL + eR)/2.0 + sqrt(h/g)*(unL - unR)/2.0;
            qu = (unL + unR)/2.0 + sqrt(g/h)*(eL - eR)/2.0;

            jacSphere = (4*R*R+x*x+y*y)/(4*R*R);
            //jacSphere = (4*R*R + xEdge[k]*xEdge[k] + yEdge[k]*yEdge[k])/(4*R*R);
            
            for (i=0; i < 2; i++) {
                FE[mapEdge[0][i]] -= qu*h*phiEdge[i]*jacSphere*jac*weight; 
                FU[mapEdge[0][i]] -= qe*nx*g*phiEdge[i]*jacSphere*jac*weight;  
                FV[mapEdge[0][i]] -= qe*ny*g*phiEdge[i]*jacSphere*jac*weight;  
                FE[mapEdge[1][i]] += qu*h*phiEdge[i]*jacSphere*jac*weight; 
                FU[mapEdge[1][i]] += qe*nx*g*phiEdge[i]*jacSphere*jac*weight;  
                FV[mapEdge[1][i]] += qe*ny*g*phiEdge[i]*jacSphere*jac*weight;  
                }
            }}
}

void femTsunamiMultiplyInverseMatrix(double* FE, double* FU, double* FV, double* X, double* Y,int* Elem,int nElem)
{        
    double BEloc[3],BUloc[3],BVloc[3];
 
    double xLoc[3],yLoc[3],jac;
    int    elem,i,j,mapElem[3];
    
    for (elem=0; elem < nElem; elem++) {
        //femShallowTriangleMap(myProblem,elem,mapElem);
        int *mapCoord = &(Elem[elem*3]);
        for (j=0; j < 3; ++j) {
        	  xLoc[j] = X[mapCoord[j]];
        	  yLoc[j] = Y[mapCoord[j]]; 
              mapElem[j] = elem*nTri+j; //celci c'est femShallowTriangleMap(elem,mapElem);
              }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        for (i=0; i < 3; i++) {
            BEloc[i] = FE[mapElem[i]];
            BUloc[i] = FU[mapElem[i]];
            BVloc[i] = FV[mapElem[i]];
            FE[mapElem[i]] = 0.0; 
            FU[mapElem[i]] = 0.0; 
            FV[mapElem[i]] = 0.0; }
        for (i=0; i < 3; i++) { 
            for (j=0; j < 3; j++) {
                FE[mapElem[i]] += invA[i][j] * BEloc[j] / jac; 
                FU[mapElem[i]] += invA[i][j] * BUloc[j] / jac; 
                FV[mapElem[i]] += invA[i][j] * BVloc[j] / jac; }}}
}

void femCompute(int size, int* Elem, double dt, double* bath,double* E, double* U, double* V, double* FE, double* FU, double* FV, double* X, double* Y, int nElem, int nEdge, femEdge* edges)
{ 
    int i;
    for (i=0; i < size+1; i++) {
        FE[i] = 0.0;
        FU[i] = 0.0;
        FV[i] = 0.0; }

    //compute vecteur b       
    femTsumnamiIntegralsElements(X, Y, E, U, FE, FU, FV, V, bath, Elem, nElem);
    femTsumnamiIntegralsEdges(X, Y, E, U, V, FE, FU, FV, bath, nElem, Elem, nEdge, edges);
    //resoudre system
    femTsunamiMultiplyInverseMatrix(FE, FU, FV, X, Y, Elem, nElem);  

    //euler explicit
    for (i=0; i < size+1; i++) {
        E[i] += dt * FE[i];
        U[i] += dt * FU[i];
        V[i] += dt * FV[i]; } 

}

void legendreUpdateRungeKutta(int size, double* E, double* U, double* V,double dt)
{
    int nStage = 2;
    double *K = calloc(size ,sizeof(double));
    double *Uold = malloc(sizeof(double) * size); memcpy(Uold, U, sizeof(double)*size);
    double *Unext = malloc(sizeof(double) * size);

    //mettre en global comme static
    const double B[2] = {0.0, 1.0};
    const double G[2] = {1.0/2.0, 1.0/2.0};

    for (int k = 0; k < nStage; k++)
    {
        for (int i = 0; i < size; i++)
        {
            Unext[i] = Uold[i] + dt*B[k]*K[i];
        }
        //legendreComputeRightHandSide(myProblem, Unext, K);
        
        for (int i = 0; i < size; i++)
        {
            U[i] += dt*G[k]*K[i];
        }
    }
    free(Uold);
    free(Unext);
    free(K);
}

/*! \brief Compute the CFD simulation of the Tsunami
 * 
 * \param[in] dt    time step
 * \param[in] nmax  maximal iteration number
 * \param[in] sub   iteration step between 2 backup
 * \param[out] void
 */

void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 
    int i,j,nNode,nElem,trash,nEdge,nBound;
    double BathMax = 9368;
    
    // Init : batimetry, X and Y by reading "meshFileName"
    FILE* file = fopen(meshFileName,"r");
    fscanf(file, "Number of nodes %d \n",&nNode);   
    double *bath = malloc(sizeof(double)*nNode);
    double *X = malloc(sizeof(double)*nNode);
    double *Y = malloc(sizeof(double)*nNode);
    for (i = 0; i < nNode; i++){
        fscanf(file,"%d : %le %le %le\n",&trash,&X[i],&Y[i],&bath[i]);
        bath[i] = (bath[i] > BathMax) ? BathMax : bath[i];
    }
    
    // Init : Nodes of each triangle
    fscanf(file, "Number of triangles %d \n",&nElem); 
    int *elem = malloc(sizeof(int)*(3*nElem+1));
    for (i = 0; i < nElem; i++){
        fscanf(file,"%d : %d %d %d \n",&trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]);
    }
    
    // Init : Edges of each triangle
    fscanf(file, "Number of edges %d \n",&nEdge);
    femEdge *edges = (femEdge*) malloc(sizeof(femEdge)*nEdge);
    nBound = 0;
    for (i = 0; i < nEdge; i++){
        fscanf(file,"%d : %d %d : %d %d \n",&trash,&(edges[i].node[0]),&(edges[i].node[1]),&(edges[i].elem[0]), &(edges[i].elem[1]) );
        if (edges[i].elem[1] == -1) nBound++;
    }
    fclose(file);

    int size = nElem*3;
    double *E  = calloc(size+1,sizeof(double));
    double *U  = calloc(size+1,sizeof(double));
    double *V  = calloc(size+1,sizeof(double));
    
    // Init : Initial condition
    for (i = 0; i < nElem*nTri; i++){
        E[i] = tsunamiInitialConditionOkada(X[elem[i]],Y[elem[i]]);
        U[i] = 0;
        V[i] = 0;
    }

    for (i = 0; i < nTri; i++){
        phi2f(gaussTriangleXsi[i], gaussTriangleEta[i], phi2[i]);
    }
    
    for (i = 0; i < nEdg; i++){
        phi1f(gaussEdgeXsi[i],phi1[i]);
    }
    
    double *FE  = calloc(size+1,sizeof(double));
    double *FU  = calloc(size+1,sizeof(double));
    double *FV  = calloc(size+1,sizeof(double));

    // Computation of the simulation
    tsunamiWriteFile(baseResultName,0,U,V,E,nElem,3);
    for (i = 1; i < nmax; i++)
    {   
        femCompute(size, elem, dt, bath, E, U, V, FE, FU, FV, X, Y, nElem, nEdge, edges);
        if (!(i%sub)){
            tsunamiWriteFile(baseResultName,i,U,V,E,nElem,3);
        }
    }

    free(bath);
    free(E);
    free(V);
    free(U);
    free(FU);
    free(FE);
    free(FV);
    free(elem);
    free(X);
    free(Y);
    free(edges);
}