#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dpBayesianNetwork.h"
#include "laa.h"

/*
連想配列
探索効率は線形
*/

Associate* makeLAA(int arraySize,int patternSize){
	Associate *result=Calloc(Associate,1);
	result->keys=Calloc(int,arraySize*patternSize);
	result->patternSize=patternSize;
	result->array=Calloc(void*,arraySize);
	return result;
}
int searchLAA(Associate* laa,int *pattern){
	int i,j;
	for(i=0;i<laa->numKeys;i++){
		for(j=0;j<laa->patternSize;j++){
			if(laa->keys[i*laa->patternSize+j]!=pattern[j]) break;
		}
		if(j==laa->patternSize){
			return i;
		}
	}
	return -1;
}
void* getLAA(Associate* laa,int *pattern){
	int i,j;
	for(i=0;i<laa->numKeys;i++){
		for(j=0;j<laa->patternSize;j++){
			if(laa->keys[i*laa->patternSize+j]!=pattern[j]) break;
		}
		if(j==laa->patternSize){
			return laa->array[i];
		}
	}
	return NULL;
}
void setLAA(Associate* laa,int *pattern,void* data){
	int i,j;
	for(i=0;i<laa->numKeys;i++){
		for(j=0;j<laa->patternSize;j++){
			if(laa->keys[i*laa->patternSize+j]!=pattern[j]) break;
		}
		if(j==laa->patternSize){
			break;
		}
	}
	laa->array[i]=data;
	if(i==laa->numKeys){
		memcpy(&laa->keys[i*laa->patternSize],pattern,sizeof(int)*laa->patternSize);
		laa->numKeys++;
	}
}
void freeLAA(Associate* laa){
	free(laa->keys);
	free(laa->array);
}
