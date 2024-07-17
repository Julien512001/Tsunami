# include "tsunami.h"
// Ce projet est le fruit du travail de Emma Stiennon et Julien Moulis


typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

double interpolate(double *phi, double *U, int *map, int n)
{
    double res = 0.0;
    for (int i = 0; i < n; i++) res += U[map[i]]*phi[i];
    return res;
}

void triangleMap(int index, int map[3])
{
    for (int i = 0; i < 3; i ++) map[i] = index*3 + i;
}

void edgesMap(femEdges *edges, int nElem, int *elem, int index, int map[2][2])
{
    for (int j = 0; j < 2; j++) {
        int node = edges->edges[index].node[j];
        for (int k = 0; k < 2; k++) {
            int ielem = edges->edges[index].elem[k];
            map[k][j] = (nElem)*3;
            if (ielem >= 0) {
                for (int i = 0; i < 3; i++) {
                    if (elem[ielem*3 + i] == node) map[k][j] = ielem*3 + i;  }}}}
}

void integralsElements(int nElem, int *elem, double *X, double *Y, double *H, double *E, double *U, double *V, double *BE, double *BU, double *BV, int nEdge, femEdges *edges,int *check)
{
    int nInteg = 3, mapElem[3];
    double dphidx[3], dphidy[3],dphidxsi[3],dphideta[3],phiElem[3],xLoc[3],yLoc[3],theta[3],xsi,eta,weight,Je,dxdxsi,dxdeta,dydxsi,dydeta,u,v,e,f,x,y,h;
    for (int iElem = 0; iElem < nElem; iElem++) {
        if (check[iElem] == 0) {
            triangleMap(iElem, mapElem);
            int *mapCoord = &(elem[3*iElem]);
            for (int j=0; j < 3; ++j) {
                xLoc[j] = X[mapCoord[j]];
                yLoc[j] = Y[mapCoord[j]]; 
            }
            Je = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
            dphidx[0] = (yLoc[1] - yLoc[2])/Je;
            dphidx[1] = (yLoc[2] - yLoc[0])/Je;
            dphidx[2] = (yLoc[0] - yLoc[1])/Je;
            dphidy[0] = (xLoc[2] - xLoc[1])/Je;
            dphidy[1] = (xLoc[0] - xLoc[2])/Je;
            dphidy[2] = (xLoc[1] - xLoc[0])/Je;  
            for (int iInteg = 0; iInteg < nInteg; iInteg++) {
                xsi = gaussTriangleXsi[iInteg];
                eta = gaussTriangleEta[iInteg];
                weight = gaussTriangleWeight[iInteg];
                phiElem[0] = 1-xsi-eta;
                phiElem[1] = xsi;
                phiElem[2] = eta;
                x = interpolate(phiElem, X, mapCoord, 3);
                y = interpolate(phiElem, Y, mapCoord, 3);
                h = interpolate(phiElem, H, mapCoord, 3);
                e = interpolate(phiElem, E, mapElem, 3);
                u = interpolate(phiElem, U, mapElem, 3);
                v = interpolate(phiElem, V, mapElem, 3);
                f = 2*Omega*((4*R*R - x*x - y*y) / (4*R*R + x*x + y*y));
                for (int i = 0; i < 3; i++) {
                    BE[mapElem[i]] += weight*Je*h*( (dphidx[i]*u + dphidy[i]*v)*(4*R*R + x*x + y*y)/(4*R*R) + phiElem[i]*(x*u + y*v)/(R*R) );
                    BU[mapElem[i]] += weight*Je* ( phiElem[i]*(f*v - Gamma*u) + dphidx[i]*g*e*(4*R*R + x*x + y*y)/(4*R*R) + phiElem[i]*(g*x*e)/(2*R*R) );
                    BV[mapElem[i]] += weight*Je*( phiElem[i]*(-f*u - Gamma*v) + dphidy[i]*g*e*(4*R*R + x*x + y*y)/(4*R*R) + phiElem[i]*(g*y*e)/(2*R*R) );
                }
            }
        }
    }
        
    nInteg = 2;
    double phiEdge[2],nx,ny,norm,eL,eR,uL,uR,vL,vR,estar,unstar,unL,unR,xsiEdge,weightEdge,JeEdge,Jx,xEdge,yEdge,hEdge;
    int sizeLoc = 3;
    int size1 = nElem * sizeLoc + 1; 
    int mapEdge[2][2];
    for (int iEdge = 0; iEdge < nEdge; iEdge++) {
        edgesMap(edges, nElem, elem, iEdge, mapEdge);
        int boundary = (mapEdge[1][0] == size1-1);
        double dxdxsi = (X[edges->edges[iEdge].node[1]] - X[edges->edges[iEdge].node[0]]);
        double dydxsi = (Y[edges->edges[iEdge].node[1]] - Y[edges->edges[iEdge].node[0]]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        JeEdge = norm/2.0;
        for (int iInteg = 0; iInteg < nInteg; iInteg++) {
            xsiEdge = gaussEdgeXsi[iInteg];
            weightEdge = gaussEdgeWeight[iInteg];
            phiEdge[0] = (1-xsiEdge)/2.0;
            phiEdge[1] = (1+xsiEdge)/2.0;
            xEdge = interpolate(phiEdge, X, edges->edges[iEdge].node, 2);
            yEdge = interpolate(phiEdge, Y, edges->edges[iEdge].node, 2);
            hEdge = interpolate(phiEdge, H, edges->edges[iEdge].node, 2);
            eL = interpolate(phiEdge, E, mapEdge[0], 2);
            eR = boundary ? eL : interpolate(phiEdge, E, mapEdge[1], 2);
            uL = interpolate(phiEdge, U, mapEdge[0], 2);
            uR = interpolate(phiEdge, U, mapEdge[1], 2);
            vL = interpolate(phiEdge, V, mapEdge[0], 2);
            vR = interpolate(phiEdge, V, mapEdge[1], 2);
            unL = uL*nx + vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            unstar = (unL+unR)/2.0 + sqrt(g/hEdge)*(eL - eR)/2.0;
            estar = (eL + eR)/2.0 + sqrt(hEdge/g)*(unL - unR)/2.0;
            Jx = (1+xEdge*xEdge/(4*R*R)+yEdge*yEdge/(4*R*R));
            for (int i = 0; i < 2; i++) {
                BE[mapEdge[0][i]] -= JeEdge*weightEdge*phiEdge[i]*hEdge*unstar*Jx;
                BU[mapEdge[0][i]] -= JeEdge*weightEdge*phiEdge[i]*nx*g*estar*Jx;
                BV[mapEdge[0][i]] -= JeEdge*weightEdge*phiEdge[i]*ny*g*estar*Jx;
                BE[mapEdge[1][i]] += JeEdge*weightEdge*phiEdge[i]*hEdge*unstar*Jx;
                BU[mapEdge[1][i]] += JeEdge*weightEdge*phiEdge[i]*nx*g*estar*Jx;
                BV[mapEdge[1][i]] += JeEdge*weightEdge*phiEdge[i]*ny*g*estar*Jx;
            }
        }
    }
}

void multiplyMatrix(int nElem, int *elem, double *X, double *Y, double *BE, double *BU, double *BV)
{
    double AInv[3][3] = {   {18.0, -6.0, -6.0},
                            {-6.0, 18.0, -6.0},
                            {-6.0, -6.0, 18.0}
                        };
    double xLocMat[3], yLocMat[3], JeMat, BEloc[3], BUloc[3], BVloc[3];
    int map[3];
    for (int iElem = 0; iElem < nElem; iElem++) {
        triangleMap(iElem, map);
        int *mapCoord = &(elem[iElem*3]);
        JeMat = fabs((X[mapCoord[1]] - X[mapCoord[0]])*(Y[mapCoord[2]] - Y[mapCoord[0]]) - (Y[mapCoord[1]] - Y[mapCoord[0]])*(X[mapCoord[2]] - X[mapCoord[0]]));
        for (int i = 0; i < 3; i++) {
            BEloc[i] = BE[map[i]];
            BUloc[i] = BU[map[i]];
            BVloc[i] = BV[map[i]];
            BE[map[i]] = 0.0;
            BU[map[i]] = 0.0;
            BV[map[i]] = 0.0;
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                BE[map[i]] += AInv[i][j]*BEloc[j]/JeMat;
                BU[map[i]] += AInv[i][j]*BUloc[j]/JeMat;
                BV[map[i]] += AInv[i][j]*BVloc[j]/JeMat;
            }
        }
    }
}

void rungeKutta(double dt,int size,double *E ,double *U,double *V,int nElem,int *elem,double *X,double *Y,double *H,double *BE,double *BU,double *BV,int nEdge,femEdges *edges,int *check) 
{

    double *EUVold = malloc(sizeof(double)*3*size);
    double *Eold = EUVold;
    double *Uold = EUVold+size;
    double *Vold = EUVold+2*size;
    double *EUVpred = malloc(sizeof(double)*3*size);
    double *Epred = EUVpred;
    double *Upred = EUVpred+size;
    double *Vpred = EUVpred+2*size;
    for (int i=0; i < size; i++) {
        BE[i] = 0.0;
        BU[i] = 0.0;
        BV[i] = 0.0;
        Eold[i] = E[i];
        Uold[i] = U[i];  
        Vold[i] = V[i];
    }
    int n = 2;
    double beta[2] = {1.0,1.0}; 
    double w[2] = {1.0/2.0,1.0/2.0};
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < size; i++) {
            Epred[i] = Eold[i] + dt*beta[j]*BE[i];
        	Upred[i] = Uold[i] + dt*beta[j]*BU[i];
            Vpred[i] = Vold[i] + dt*beta[j]*BV[i];
        }
    	integralsElements(nElem, elem, X,Y,H,Epred,Upred,Vpred,BE,BU,BV, nEdge, edges, check);
        multiplyMatrix(nElem,elem,X,Y,BE,BU,BV);
        for (int i=0; i < size; i++) {
            E[i] += dt*w[j]*BE[i];
        	U[i] += dt*w[j]*BU[i]; 
            V[i] += dt*w[j]*BV[i];
        }
    }
    free(EUVold);
    free(EUVpred);
}

void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 

    int nNodes, nElem, trash, i, nEdge;
    double *X, *Y, *H;

    FILE *file = fopen(meshFileName, "r");
    fscanf(file, "Number of nodes %d \n", &nNodes);
    X = malloc(sizeof(double)*nNodes);
    Y = malloc(sizeof(double)*nNodes);
    H = malloc(sizeof(double)*nNodes);
    for (i = 0; i < nNodes; i++) fscanf(file, "%d : %le %le %le \n", &trash, &X[i], &Y[i], &H[i]);
    

    fscanf(file, "Number of triangles %d \n", &nElem);
    int *elem = (int*) malloc(3*nElem*sizeof(int));
    for (i = 0; i < nElem; i++) fscanf(file, "%d : %d %d %d \n", &trash, &elem[3*i], &elem[3*i+1], &elem[3*i+2]);
    

    fscanf(file, "Number of edges %d \n", &nEdge);
    femEdge *edge = malloc(sizeof(femEdge)*nEdge);
    femEdges *edges = malloc(sizeof(femEdges));
    edges->edges = edge;
    edges->nEdge = nEdge;
    for (i = 0; i < nEdge; i++) fscanf(file, "%d : %d %d : %d %d \n", &trash, &edges->edges[i].node[0], &edges->edges[i].node[1], &edges->edges[i].elem[0], &edges->edges[i].elem[1]);
    
    fclose(file);

    int size = nElem*3+1;
    double *EUV = calloc(3*size,sizeof(double));
    double *E = EUV;
    double *U = EUV+size;
    double *V = EUV+2*size;

    double *BEBUBV = calloc(3*size,sizeof(double));
    double *BE = BEBUBV;
    double *BU = BEBUBV+size;
    double *BV = BEBUBV+2*size;
    for (int i = 0; i < size; i++) 
        E[i] = tsunamiInitialConditionOkada(X[elem[i]],Y[elem[i]]);
    
    int *check = (int*) calloc(size, sizeof(int));

    for (int k = 1; k < nmax+1; k++)
    {   
        rungeKutta(dt,size,E,U,V,nElem,elem,X,Y,H,BE,BU,BV,nEdge,edges,check);
        for (int j = 0; j < size; j++) {
            if (H[j] < 1e6) check[j] = 1;
        }
        if (!(k%sub)) tsunamiWriteFile(baseResultName,k,U,V,E,nElem,3);
        
    }
    free(X);
    free(Y);
    free(H);
    free(elem);
    free(edges->edges);
    free(edges);
    free(EUV);
    free(BEBUBV);
}
