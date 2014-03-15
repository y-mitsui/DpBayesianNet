#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "laa.h"
#include "dpBayesianNetwork.h"



/*
	各ノードごとに取りうる全ての親パターンで学習スコアを計算する
*/
static void getLocalScore(baysianNode *node,int numNode,int currentNode,int numLine,Associate *cHist,unsigned short parentPattern,double *scores,Associate **CPT,int numSample,int rank,int start){
	int *pattern;
	int *nums;
	Associate *tmp;
	double like;
	int i,j,patternCt,total;

	if(numLine==rank){
		/* 子ノードcurrentNodeにおける親パターンedgeについて、条件付き確率を計算する */
		tmp=makeLAA(cHist->numKeys,numLine);	//親値パターンと子ノード値ごとの総数
		pattern=Calloc(int,numNode);
		for(i=0;i<cHist->numKeys;i++){	
			for(patternCt=0,j=0;j<numNode;j++){
				if((parentPattern&((unsigned short)1<<j))!=0){
					pattern[patternCt++]=cHist->keys[i*cHist->patternSize+j];
				}
			}
			if(!(nums=(int*)getLAA(tmp,pattern))) nums=Calloc(int,node->numVariablePattern);
			for(j=0;j<node->numVariablePattern;j++){
				nums[j]+=((int*)cHist->array[i])[j];
			}
			setLAA(tmp,pattern,nums);
		}
		/* 条件付き確立を計算*/
		CPT[parentPattern]=makeLAA(tmp->numKeys,numLine);
		
		for(i=0;i<tmp->numKeys;i++){
			double *prop=Malloc(double,node->numVariablePattern);
			for(total=0,j=0;j<node->numVariablePattern;j++){
				total+=((int*)tmp->array[i])[j];
			}
			for(j=0;j<node->numVariablePattern;j++){
				prop[j]=(double)((int*)tmp->array[i])[j]/total;
				if(prop[j]==0.0) prop[j]=0.00000001;
			}			
			setLAA(CPT[parentPattern],&tmp->keys[i*tmp->patternSize],prop);
		}
		
		/* 尤度を計算 */
		like=0.0;
		total=0;
		for(i=0;i<tmp->numKeys;i++){
			double* prop=(double*)getLAA(CPT[parentPattern],&tmp->keys[i*tmp->patternSize]);
			for(j=0;j<node->numVariablePattern;j++){
				like+=log(prop[j])*((int*)tmp->array[i])[j];
				total+=((int*)tmp->array[i])[j];
			}
			free(tmp->array[i]);
		}
		if(like==0.0) like=-0.00000000001;
		scores[parentPattern]=1.0/(-2*like+(numLine*(node->numVariablePattern-1))*log(numSample));//BICの逆数
		free(pattern);
		freeLAA(tmp);
		free(tmp);
	}else{
		for(i=start;i<numNode;i++){
			if(i!=currentNode){
				getLocalScore(node,numNode,currentNode,numLine,cHist,parentPattern | (1<<i),scores,CPT,numSample,rank+1,i+1);
			}
		}
	}
}
static void getBestParents(int currentNode,int numNode,double *localScore,unsigned short *bestParents,unsigned short parentPattern,int num,int start,int rank){
	double maxScore;
	unsigned short maxScorePattern=0;
	int i;

	if(num==rank){
		unsigned short oneLeaveParentPattern=parentPattern;
		for(maxScore=0,i=0;i<numNode;i++){
			if(i!=currentNode){
				oneLeaveParentPattern=oneLeaveParentPattern & ~(1<<i);
				if(maxScore < localScore[oneLeaveParentPattern]){
					maxScore=localScore[oneLeaveParentPattern];
					maxScorePattern=oneLeaveParentPattern;
				}
			}
		}
		bestParents[parentPattern]=(localScore[parentPattern] > localScore[maxScorePattern]) ? parentPattern : maxScorePattern;
		
	}else{
		for(i=start;i<numNode;i++){
			if(i!=currentNode){
				getBestParents(currentNode,numNode,localScore,bestParents,parentPattern | (1<<i),num,i+1,rank+1);
			}
		}
	}
	
}
static void getBestLeaf(double *scores,int *leafs,unsigned short **bestParents,int *nodes,unsigned short nodePattern,double **localScore,int numNode,int num,int start,int rank){
	int i;
	if(num==rank){
		scores[nodePattern]=0;
		leafs[nodePattern]=0;
		for(i=0;i<num;i++){
			unsigned short oneLeavePattern= nodePattern & ~(1<<nodes[i]);
			double scoreTmp=scores[oneLeavePattern] + localScore[nodes[i]][bestParents[nodes[i]][oneLeavePattern]];
			if(scoreTmp > scores[nodePattern]){
				scores[nodePattern]=scoreTmp;
				leafs[nodePattern]=nodes[i];
			}
		}
	}else{
		for(i=start;i<numNode;i++){
			nodes[rank]=i;
			getBestLeaf(scores,leafs,bestParents,nodes,nodePattern | (1<<i),localScore,numNode,num,i+1,rank+1);
		}
	}
}
static void leafs2Order(int *order,int *leafs,int numNode){
	int i;
	unsigned short nodePattern=(1<<numNode)-1;
	for(i=numNode-1;i>-1;i--){
		order[i]=leafs[nodePattern];
		nodePattern=nodePattern & ~(1<<order[i]);
	}
	
}
static void Ord2net(unsigned short *parents,int *order,unsigned short **bestParents,int numNode){
	unsigned short predecs=0;
	int i;
	for(i=0;i<numNode;i++){
		parents[order[i]]=bestParents[order[i]][predecs];
		predecs|=1<<order[i];
	}
}
static void net2Matrix(unsigned short *parents,int numNode,int *edge){
	int i,j;
	for(i=0;i<numNode;i++){
		for(j=0;j<numNode;j++){
			if((parents[i]&(1<<j))!=0){
				edge[i*numNode+j]=1;
			}
		}
	}
}
int *int2Number(int num){
	int *r=Malloc(int,1);
	*r=num;
	return r;
}
int combination(int n,int m){
	double sum=1.0,sum2=1.0;
	int i;
	for(i=0;i<m;i++){
		sum*=n-i;
	}
	for(i=0;i<m;i++){
		sum2*=(m-i);
	}
	return (int)(sum/sum2);
}
bayesianNetwork* bayesianNetTrain(int *data,int numSample,int dimention,bayesianNetworkOption *opt){
	Associate *frequency;
	int i,j,k,num;
	int *pattern;
	bayesianNetwork *model;

	/* 分割表を作成 */
	frequency=makeLAA(numSample,opt->numVariable);
	pattern=Malloc(int,opt->numVariable);
	for(i=0;i<numSample;i++){
		for(j=0;j<dimention;j++){
			pattern[j]=data[i*dimention+j];
		}

		int* cur=(int*)getLAA(frequency,pattern);
		if(cur==NULL){
			setLAA(frequency,pattern,int2Number(1));
		}else{
			(*cur)++;
		}
	}
	/* 条件付き頻度表を作成*/	
	Associate **cHist=Malloc(Associate*,opt->numVariable);
	for(i=0;i<opt->numVariable;i++){
		cHist[i]=makeLAA(numSample,opt->numVariable);
		for(j=0;j<frequency->numKeys;j++){
			for(k=0;k<frequency->patternSize;k++){
				pattern[k]=(i==k) ? 0 : frequency->keys[j*frequency->patternSize+k];
			}
			int val=frequency->keys[j*frequency->patternSize+i];
			int *cur=(int*)getLAA(cHist[i],pattern);
			if(cur) cur[val]+=*((int*)frequency->array[j]);
			else{
				cur=Calloc(int,opt->node[i].numVariablePattern);
				cur[val]=*((int*)frequency->array[j]);
			}
			setLAA(cHist[i],pattern,cur);
		}
	}	
	/* ローカルスコアを計算*/
	Associate ***CPT=Malloc(Associate **,opt->numVariable);
	double **localScore=Malloc(double *,opt->numVariable);
	for(i=0;i<opt->numVariable;i++){
		localScore[i]=Malloc(double,pow(2,opt->numVariable));
		CPT[i]=Malloc(Associate *,pow(2,opt->numVariable));
		for(j=0;j<opt->numVariable;j++){
			getLocalScore(&opt->node[i],opt->numVariable,i,j,cHist[i],0,localScore[i],CPT[i],numSample,0,0);
		}
	}
	model=Malloc(bayesianNetwork,1);
	model->CPT=CPT;
	model->nodes=opt->node;
	model->numNode=opt->numVariable;
	unsigned short **bestParents=Malloc(unsigned short *,opt->numVariable);
	for(i=0;i<opt->numVariable;i++){
		bestParents[i]=Malloc(unsigned short,pow(2,opt->numVariable));
		for(j=0;j<opt->numVariable;j++){
			getBestParents(i,opt->numVariable,localScore[i],bestParents[i],0,j,0,0);
		}
	}

	
	for(num=0,i=0;i<=opt->numVariable;i++) num+=combination(opt->numVariable,i);
	double *scores=Malloc(double,pow(2,opt->numVariable));
	int *leafs=Malloc(int,pow(2,opt->numVariable));
	int *nodes=Malloc(int,opt->numVariable);
	for(i=0;i<=opt->numVariable;i++){
		getBestLeaf(scores,leafs,bestParents,nodes,0,localScore,opt->numVariable,i,0,0);
	}

	int* order=Malloc(int,opt->numVariable);
	leafs2Order(order,leafs,opt->numVariable);
	unsigned short *parents=Malloc(unsigned short,opt->numVariable);
	Ord2net(parents,order,bestParents,opt->numVariable);
	int *edge=Calloc(int,opt->numVariable*opt->numVariable);
	net2Matrix(parents,opt->numVariable,edge);

	//printSquereMatrix(edge,opt->numVariable*opt->numVariable);

	for(i=0;i<opt->numVariable;i++){
		free(localScore[i]);
		for(j=0;j<cHist[i]->numKeys;j++){
			free(cHist[i]->array[j]);
		}
		freeLAA(cHist[i]);
		free(cHist[i]);
	}
	free(localScore);
	free(cHist);
	for(i=0;i<frequency->numKeys;i++)
		free(frequency->array[i]);
	freeLAA(frequency);
	free(frequency);
	free(pattern);
	
	return model;

}

static double getJointProbability(bayesianNetwork *model,int targetNode,int targetValue,int *values,int *edge,int *nodes,int numNode,int *valuePattern,int maxNode,int rank){
	double result,prob;
	int *vals,numParents,i,j;
	unsigned short parentsPattern;
	if(maxNode==rank){
		prob=1.0;
		vals=Malloc(int,numNode);
		for(i=0;i<maxNode;i++){
			memset(vals,0,sizeof(int)*numNode);
			if(targetNode!=i && values[i]!=-1 && valuePattern[i]!=values[i]){ prob=0.0; goto FINISH; }
			for(parentsPattern=numParents=0,j=0;j<maxNode;j++){			
				if(edge[i*maxNode+j]==1){
					parentsPattern|=1<<j;
					vals[numParents++]=(j==targetNode) ? targetValue : valuePattern[j];
					if(targetNode!=i && values[j]!=-1 && valuePattern[j]!=values[j]){ prob=0.0; goto FINISH; }
				}
			}
			if(i!=targetNode){
				prob*=((double*)getLAA(model->CPT[i][parentsPattern],vals))[valuePattern[i]];
			}else{
				prob*=((double*)getLAA(model->CPT[targetNode][parentsPattern],vals))[targetValue];
			}
		}
FINISH:
		free(vals);
		return prob;
	}
	if(rank==targetNode){
		return getJointProbability(model,targetNode,targetValue,values,edge,nodes,numNode+1,valuePattern,maxNode,rank+1);
	}
	result=0.0;
	for(i=0;i<model->nodes[rank].numVariablePattern;i++){
		valuePattern[rank]=i;
		result+=getJointProbability(model,targetNode,targetValue,values,edge,nodes,numNode+1,valuePattern,maxNode,rank+1);

	}
	return result;
}
/* 
列挙法による確率推論
*/
static double bayesianNetworkGetProbability(bayesianNetwork *model,int targetNode,int targetValue,int *values){
	int *nodes,*valuePattern;
	nodes=Calloc(int,model->numNode);
	valuePattern=Calloc(int,model->numNode);
	double p1=getJointProbability(model,targetNode,targetValue,values,model->models[0]->edge,nodes,0,valuePattern,model->numNode,0);
	memset(valuePattern,0,sizeof(int)*model->numNode);
	printf("p1:%.15lf\n",p1);
	double p2=getJointProbability(model,-1,targetValue,values,model->models[0]->edge,nodes,0,valuePattern,model->numNode,0);
	printf("p2:%.15lf p1/p2:%.15lf\n",p2,p1/p2);
	return p1/p2;
}

void *bayesianNetPredict(bayesianNetwork *model,int *values,int targetNode){
	double maxProb=0.0;
	int maxIdx,i;

	values[targetNode]=-1;
	for(i=0;i<model->nodes[targetNode].numVariablePattern;i++){
		double prob=bayesianNetworkGetProbability(model,targetNode,i,values);
		printf("final prob:%.15lf\n",prob);
		if(prob > maxProb){
			maxProb=prob;
			maxIdx=i;
		}
	}
	return int2Number(maxIdx);
}
