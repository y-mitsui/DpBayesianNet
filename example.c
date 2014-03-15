#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dpBayesianNetwork.h"

/*
真の構造が B→A となるように、4値をサンプリング
*/
void setSample(int *sample,int numSample,int numVariable){
	int i;

	for(i=0;i<numSample;i++){
		sample[i*numVariable+1]=rand()%4;
		if(sample[i*numVariable+1]<2) sample[i*numVariable]=(rand()%100 < 80) ? sample[i*numVariable+1] : rand()%4;
		else sample[i*numVariable]=rand()%4;
	}
}
#define NUM_NODE 2
#define NUM_SAMPLE 120000
int main(void){
	int numVariablePattern[NUM_NODE]={4,4};
	int *sample=Malloc(int,NUM_SAMPLE*NUM_NODE);

	setSample(sample,NUM_SAMPLE,NUM_NODE);

	bayesianNetwork *model=bayesianNetTrain(sample,NUM_SAMPLE,NUM_NODE,numVariablePattern);

	bayesianNetPrintSquereMatrix(model->edge,NUM_NODE*NUM_NODE);

	int predictValue=bayesianNetPredict(model,&sample[0],0);
	printf("predict:%d\n",predictValue);

	bayesianNetFree(model);

	free(sample);
	free(model);
	return 0;
}
