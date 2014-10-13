#ifndef LAA_H_
#define LAA_H_

typedef struct{
	int *keys;
	int numKeys;
	void **array;
	int patternSize;
}Associate;

void* getLAA(Associate* laa,int *pattern);
void setLAA(Associate* laa,int *pattern,void* data);
Associate* makeLAA(int arraySize,int patternSize);
void freeLAA(Associate* laa);

#endif
