#ifndef DP_BAYESIAN_NETWORK_H_
#define DP_BAYESIAN_NETWORK_H_

#include "laa.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define Calloc(type,n) (type *)calloc(1,(n)*sizeof(type))

typedef struct{
	int numNode;		//確率変数の数
	int *edge;		//ネットワーク構造
	int *numValuePattern;	//各ノードごとの取りうる値の数
	Associate ***CPT;	//条件付き確率表.
}bayesianNetwork;


/* bayesianNetTrain
動的計画法による学習
data:サンプリングデータ  要素数はnumSample*dimention
numSample:サンプリングデータの数
dimention:サンプリングデータの次元(確率変数の数)
numValuePattern:各ノードごとの取りうる値の数
戻り値:学習結果
*/
bayesianNetwork* bayesianNetTrain(int *data,int numSample,int dimention,int *numValuePattern);

/* bayesianNetPredict
クラス分類を行う
model:bayesianNetTrainによって得られたコンテキスト
targetNode:目的変数
values:説明変数の値。(エビデンス)
戻り値:確率が最大となるtargetNodeが取る値
*/
int bayesianNetPredict(bayesianNetwork *model,int *values,int targetNode);

/* bayesianNetFree
学習結果のリリース */
void bayesianNetFree(bayesianNetwork *model);

int combination(int n,int m);
void bayesianNetPrintSquereMatrix(int *matrix,int numMatrix);



#endif
