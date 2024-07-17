#include <stdio.h>
#include <stdlib.h>

void main (void) {


    int nNode, nElem;
    int trash;

    FILE* file = fopen("test.txt", "r");
    if (file == NULL) {
    	  printf("Error : cannot open mesh file :-) \n");
        exit(0); 
    }
    fscanf(file, "Number of nodes %d \n", &nNode);
    double *X = (double *) malloc(sizeof(double)*nNode);
    double *Y = (double *) malloc(sizeof(double)*nNode);
    for (int i = 0; i < nNode; i++) {
        fscanf(file, "%d : %le %le \n", &trash, &X[i], &Y[i]);
    }

    fscanf(file, "Number of triangles %d \n", &nElem);
    int *elem = (int *) malloc(sizeof(int)*3*nElem);
    for (int i = 0; i < nElem; i++) {
        fscanf(file, "%d : %d %d %d \n", &trash, &elem[3*i], &elem[3*i+1], &elem[3*i+2]);
    }
    for (int i = 0; i < nElem; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%d ", elem[3*i + j ]);
        }
    }

}