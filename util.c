#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void bayesianNetPrintSquereMatrix(int *matrix,int numMatrix){
	int width=sqrt(numMatrix);
	int i,j;
	puts("--------------");
	printf("  ");
	for(i=0;i<width;i++){
		printf("%c ",'A'+i);
	}
	puts("");
	for(i=0;i<width;i++){
		printf("%c ",'A'+i);
		for(j=0;j<width;j++){
			printf("%d ",matrix[i*width+j]);
		}
		puts("");
	}
	puts("--------------");
}
