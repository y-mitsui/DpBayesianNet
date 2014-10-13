#ifndef DP_BAYESIAN_NETWORK_H_
#define DP_BAYESIAN_NETWORK_H_

#include "laa.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define Calloc(type,n) (type *)calloc(1,(n)*sizeof(type))

typedef struct{
	int numNode;		//確率変数の数 number of random variable
	int *edge;		//ネットワーク構造  structure of network
	int *numValuePattern;	//各ノードごとの取りうる値の数 number of pattern of value for every node
	Associate ***CPT;	//条件付確率表. conditional probability table
}bayesianNetwork;


/* bayesianNetTrain
動的計画法による学習    Train using dynamic programming 
data:サンプリングデータ。要素数はnumSample*dimention    Sampling data. Number of elements is numSample*dimention.
numSample:サンプリングデータの数    Number of sampling data.
dimention:サンプリングデータの次元(確率変数の数)   Dimention of sampling data.(Number of random variables)
numValuePattern:各ノードごとの取りうる値の数   Number of pattern of value for each node
戻り値:学習結果     Return value: Result trained.
*/
bayesianNetwork* bayesianNetTrain(int *data,int numSample,int dimention,int *numValuePattern);

/* bayesianNetPredict
クラス分類を行う     Classification
model:bayesianNetTrainによって得られたコンテキスト   Context data which get by bayesianNetTrain
targetNode:目的変数   Target variable.
values:説明変数の値。(エビデンス)   Evidence
戻り値:確率が最大となるtargetNodeが取る値      Return value:Value which take the biggest probability.
*/
int bayesianNetPredict(bayesianNetwork *model,int *values,int targetNode);

/* bayesianNetFree
学習結果のリリース */
void bayesianNetFree(bayesianNetwork *model);

int combination(int n,int m);
void bayesianNetPrintSquereMatrix(int *matrix,int numMatrix);



#endif
