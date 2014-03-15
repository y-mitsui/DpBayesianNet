#ifndef DP_BAYESIAN_NETWORK_H_
#define DP_BAYESIAN_NETWORK_H_

#include "laa.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define Calloc(type,n) (type *)calloc(1,(n)*sizeof(type))

typedef struct{
	int numNode;
	int *edge;
	int *numValuePattern;
	Associate ***CPT;
}bayesianNetwork;



int combination(int n,int m);
void bayesianNetPrintSquereMatrix(int *matrix,int numMatrix);
int bayesianNetPredict(bayesianNetwork *model,int *values,int targetNode);
bayesianNetwork* bayesianNetTrain(int *data,int numSample,int dimention,int *numValuePattern);
void bayesianNetFree(bayesianNetwork *model);

#endif
