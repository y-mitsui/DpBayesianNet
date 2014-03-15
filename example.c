#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dpBayesianNetwork.h"

/*
真の構造が B→A となるように、サンプリング
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
	baysianNode node[NUM_NODE];
	bayesianNetworkOption opt;
	int *sample=Malloc(int,NUM_SAMPLE*NUM_NODE);
	setSample(sample,NUM_SAMPLE,NUM_NODE);
	memset(&node,0,sizeof(node[0])*NUM_NODE);
	node[0].name="A";
	node[0].numVariablePattern=4;
	node[1].name="B";
	node[1].numVariablePattern=4;


	opt.numVariable=NUM_NODE;
	opt.node=node;
	bayesianNetwork *model=bayesianNetTrain(sample,NUM_SAMPLE,NUM_NODE,&opt);
	
	int *value=bayesianNetPredict(model,&sample[0],0);
	printf("value:%d\n",*value);
	return 0;
}
