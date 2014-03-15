#ifndef DP_BAYESIAN_NETWORK_H_
#define DP_BAYESIAN_NETWORK_H_

#include "laa.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define Calloc(type,n) (type *)calloc(1,(n)*sizeof(type))

typedef struct _baysianNode{
	const char *name;
	int value;
	int numVariablePattern;
	double *probability;
}baysianNode;

typedef struct {
	int *edge;
	double score;
}bayesianNetScore;
typedef struct{
	bayesianNetScore **models;
	baysianNode *nodes;
	int numModels;
	int numNode;
	Associate ***CPT;
}bayesianNetwork;
typedef struct{
	baysianNode *node;
	int numVariable;
}bayesianNetworkOption;


int combination(int n,int m);
void printSquereMatrix(int *matrix,int numMatrix);
void *bayesianNetPredict(bayesianNetwork *model,int *values,int targetNode);
bayesianNetwork* bayesianNetTrain(int *data,int numSample,int dimention,bayesianNetworkOption *opt);

#endif
